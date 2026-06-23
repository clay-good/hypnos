"""Validate the dataset against its JSON Schema and check cross-references.

Two layers:

* **Schema validation** — every model record conforms to ``model.schema.json``.
* **Integrity checks** — referential and semantic invariants the schema cannot
  express on its own (citations resolve, kernel bindings exist, record tier is
  the worst contributing parameter tier, envelopes are well-ordered).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

from .load import Dataset, find_dataset_dir, load
from .models import TIER_RANK, worst_tier


class ValidationError(Exception):
    pass


def _load_schema(root: Path, name: str = "model.schema.json") -> dict:
    schema_path = root / "schema" / name
    with schema_path.open(encoding="utf-8") as fh:
        return json.load(fh)


def validate_dataset(ds: Optional[Dataset] = None) -> List[str]:
    """Return a list of problem strings. Empty list == dataset is valid."""
    if ds is None:
        ds = load()
    root = find_dataset_dir()
    problems: List[str] = []

    # --- schema layer -----------------------------------------------------
    try:
        import jsonschema

        schema = _load_schema(root)
        validator = jsonschema.Draft7Validator(schema)
        for m in ds:
            for err in sorted(validator.iter_errors(m.raw), key=lambda e: e.path):
                loc = "/".join(str(p) for p in err.path)
                problems.append(f"[schema] {m.id}: {loc}: {err.message}")
        # covariate-equation library records (v0.7 §3.2) against their own schema
        eq_schema_path = root / "schema" / "covariate_equation.schema.json"
        if eq_schema_path.is_file():
            eq_validator = jsonschema.Draft7Validator(_load_schema(root, "covariate_equation.schema.json"))
            for eid, rec in getattr(ds, "covariate_equations", {}).items():
                for err in sorted(eq_validator.iter_errors(rec), key=lambda e: e.path):
                    loc = "/".join(str(p) for p in err.path)
                    problems.append(f"[schema] covariate_equation {eid}: {loc}: {err.message}")
    except ImportError:  # pragma: no cover - jsonschema is a hard dependency
        problems.append("[schema] jsonschema not installed; skipped schema validation")

    # --- integrity layer --------------------------------------------------
    from .export.registry import (  # local import to avoid cycle
        INTERACTION_KERNELS,
        KERNELS,
        LA_KERNELS,
        PD_KERNELS,
        VOLATILE_KERNELS,
    )

    known_kernels = (set(KERNELS) | set(PD_KERNELS) | set(INTERACTION_KERNELS)
                     | set(VOLATILE_KERNELS) | set(LA_KERNELS))

    known_citations = set(ds.citations.keys())
    for m in ds:
        # record-level citation resolves
        if m.primary_citation not in known_citations:
            problems.append(f"[cite] {m.id}: primary_citation '{m.primary_citation}' not in citations/")
        # per-parameter citations resolve when present
        for p in m.parameters:
            if p.primary_citation and p.primary_citation not in known_citations:
                problems.append(
                    f"[cite] {m.id}: parameter {p.symbol} cites unknown '{p.primary_citation}'"
                )
        # failure-mode citations resolve when present
        for fm in m.known_failure_modes:
            if fm.citation and fm.citation not in known_citations:
                problems.append(f"[cite] {m.id}: failure-mode cites unknown '{fm.citation}'")
        # predictive-performance citations resolve when present (a performance
        # number must be traceable to a source, never asserted bare — spec §5/§9)
        for pp in m.predictive_performance:
            cid = pp.get("citation")
            if cid and cid not in known_citations:
                problems.append(
                    f"[cite] {m.id}: predictive_performance ({pp.get('metric')}) cites unknown '{cid}'"
                )

        # record tier == worst contributing parameter tier (the "worst input wins" invariant)
        param_tiers = [p.tier for p in m.parameters]
        if param_tiers:
            expected = worst_tier(param_tiers)
            if TIER_RANK[m.tier] < TIER_RANK[expected]:
                problems.append(
                    f"[tier] {m.id}: record tier {m.tier} is better than worst "
                    f"parameter tier {expected} (worst input must win)"
                )

        # kernel binding resolves to a registered function
        if m.kernel_implemented:
            fn = m.kernel_function
            if fn not in known_kernels:
                problems.append(f"[kernel] {m.id}: kernel.function '{fn}' not registered")

        # envelope ranges well-ordered (min <= max)
        env = m.applicability_envelope
        for name in ("age_years", "weight_kg", "height_cm", "bmi_kg_m2",
                     "crcl_ml_min", "albumin_g_dl", "ejection_fraction_pct",
                     "pma_weeks", "postnatal_age_days", "gestational_age_weeks"):
            rng = getattr(env, name)
            if rng.min is not None and rng.max is not None and rng.min > rng.max:
                problems.append(f"[envelope] {m.id}: {name} min {rng.min} > max {rng.max}")

        # organ-tolerance citations resolve (a cited standing claim must be traceable —
        # the same never-assert-bare rule the rest of the dataset follows; v0.5 §C)
        for ot in env.organ_tolerance:
            cid = ot.get("citation")
            if cid and cid not in known_citations:
                problems.append(
                    f"[cite] {m.id}: organ_tolerance ({ot.get('axis')}) cites unknown '{cid}'"
                )

        # local-anesthetic absorption citations resolve, and site ranks are unique
        # (the rank is the robust, curated direction — duplicates would be a slip)
        absn = m.absorption
        if absn:
            ac = absn.get("primary_citation")
            if ac and ac not in known_citations:
                problems.append(f"[cite] {m.id}: absorption cites unknown '{ac}'")
            ranks = [s.get("rank") for s in absn.get("site_rates", [])]
            for s in absn.get("site_rates", []):
                sc = s.get("citation")
                if sc and sc not in known_citations:
                    problems.append(
                        f"[cite] {m.id}: absorption site '{s.get('site')}' cites unknown '{sc}'"
                    )
            if len(ranks) != len(set(ranks)):
                problems.append(f"[absorption] {m.id}: duplicate site_rates rank(s) {ranks}")

        # --- local-anesthetic toxicity thresholds (v0.6 LA1 spec §4) ------
        problems.extend(_check_toxicity_thresholds(m, ds, known_citations))

        # id prefix matches declared subsystem
        if m.id.split(".")[0] != m.subsystem:
            problems.append(f"[id] {m.id}: id prefix does not match subsystem '{m.subsystem}'")

        # review-state integrity (the pending_human_review governance line)
        problems.extend(_check_review_state(m))

        # --- variability layer (v0.2 spec §9) -----------------------------
        problems.extend(_check_variability(m, known_citations))

        # --- external-validation layer (v0.4 spec §4) ---------------------
        problems.extend(_check_external_validation(m))

        # --- tier falsification (v0.4 spec §5) — advisory, humans still gate ---
        problems.extend(_check_tier_falsification(m))

        # --- estimation-uncertainty layer (v0.3 spec §9) ------------------
        problems.extend(_check_estimation(m, known_citations))

        # --- covariate-model layer (v0.7 spec §9) -------------------------
        problems.extend(_check_covariate_model(m, ds, known_citations))

        # --- developmental layer (v0.8 spec §8) ---------------------------
        problems.extend(_check_developmental(m, known_citations))

        # --- pharmacogenomic layer (v0.9 spec §8) -------------------------
        problems.extend(_check_pharmacogenomics(m, known_citations))

    # --- covariate-equation library (v0.7 spec §9) ------------------------
    problems.extend(_check_covariate_equations(ds, known_citations))

    # drug-level protein-binding citations resolve (v0.5 §B3 binding failure mode —
    # a binding-sensitivity claim must be traceable, like every other curated claim)
    for d in getattr(ds, "drugs", {}).values():
        pb = d.get("protein_binding") or {}
        cid = pb.get("citation")
        if cid and cid not in known_citations:
            problems.append(
                f"[cite] drug {d.get('name')}: protein_binding cites unknown '{cid}'"
            )
        # drug-level free_fraction_model (v0.6 LA3): the non-linear saturable-binding
        # model. Must cite, use a known type, carry a positive capacity, and only sit on
        # a drug whose binding is actually saturable (a linear drug has no saturation).
        ffm = pb.get("free_fraction_model") or {}
        if ffm:
            fcid = ffm.get("citation")
            if fcid and fcid not in known_citations:
                problems.append(
                    f"[cite] drug {d.get('name')}: free_fraction_model cites unknown '{fcid}'"
                )
            if ffm.get("type") not in ("capacity_limited",):
                problems.append(
                    f"[free-fraction] drug {d.get('name')}: free_fraction_model type "
                    f"'{ffm.get('type')}' not supported (expected 'capacity_limited')"
                )
            cap = ffm.get("binding_capacity_ug_ml")
            if cap is None or cap <= 0:
                problems.append(
                    f"[free-fraction] drug {d.get('name')}: free_fraction_model needs a "
                    "positive binding_capacity_ug_ml"
                )
            if not pb.get("saturable"):
                problems.append(
                    f"[free-fraction] drug {d.get('name')}: free_fraction_model present but "
                    "binding is not marked saturable (a non-saturable drug has no saturation model)"
                )
        # drug-level cardiotoxicity_class (v0.6 LA2): a stereochemistry-driven claim
        # that must be traceable + use the controlled rank/margin vocabulary.
        cc = d.get("cardiotoxicity_class") or {}
        if cc:
            ccid = cc.get("citation")
            if ccid and ccid not in known_citations:
                problems.append(
                    f"[cite] drug {d.get('name')}: cardiotoxicity_class cites unknown '{ccid}'"
                )
            if cc.get("rank") not in (None, "high", "intermediate", "low"):
                problems.append(
                    f"[cardiotoxicity] drug {d.get('name')}: invalid rank '{cc.get('rank')}'"
                )
            if cc.get("cns_to_cvs_margin") not in (None, "narrow", "moderate", "wide"):
                problems.append(
                    f"[cardiotoxicity] drug {d.get('name')}: invalid cns_to_cvs_margin "
                    f"'{cc.get('cns_to_cvs_margin')}'"
                )

    return problems


def _check_toxicity_thresholds(m, ds, known_citations) -> List[str]:
    """Local-anesthetic toxicity-threshold checks (v0.6 LA1 spec §4) — the
    safety-critical invariants the schema cannot express on its own:

    * a threshold is a RANGE with ``low < high`` (no false-precision lines);
    * every threshold declares its ``basis`` (total vs free) — enforced by schema,
      re-checked here for defense in depth;
    * a ``total_plasma`` basis on a ``saturable`` drug MUST carry the free-fraction
      ``saturation_caveat`` (total under-predicts free risk exactly when risk is
      highest — the documented failure mode a naive view would hide);
    * only the ``local_anesthetics`` subsystem may carry thresholds;
    * citations resolve.
    """
    problems: List[str] = []
    thresholds = m.toxicity_thresholds
    if not thresholds:
        return problems
    if m.subsystem != "local_anesthetics":
        problems.append(
            f"[toxicity] {m.id}: toxicity_thresholds present on a non-local_anesthetics "
            f"subsystem '{m.subsystem}'")
    drug = ds.drug(m.drug_name) if hasattr(ds, "drug") else (
        getattr(ds, "drugs", {}) or {}).get(m.drug_name)
    saturable = bool(((drug or {}).get("protein_binding") or {}).get("saturable"))
    for t in thresholds:
        if t.low is None or t.high is None:
            problems.append(
                f"[toxicity] {m.id}: {t.endpoint} threshold is not a range "
                "(low/high required — a single-value threshold is forbidden)")
        elif t.low >= t.high:
            problems.append(
                f"[toxicity] {m.id}: {t.endpoint} threshold low {t.low} >= high {t.high} "
                "(must be a range with low < high — no false-precision line)")
        if t.basis not in ("total_plasma", "free_plasma"):
            problems.append(
                f"[toxicity] {m.id}: {t.endpoint} threshold has invalid basis {t.basis!r}")
        # the load-bearing free-fraction guard (v0.6 §3.2/§4)
        if saturable and t.basis == "total_plasma" and not (t.saturation_caveat or "").strip():
            problems.append(
                f"[toxicity] {m.id}: {t.endpoint} threshold is total_plasma on a saturable-"
                "binding drug but carries no saturation_caveat (the free fraction rises "
                "non-linearly when total is highest — total under-predicts risk)")
        if t.primary_citation and t.primary_citation not in known_citations:
            problems.append(
                f"[cite] {m.id}: {t.endpoint} toxicity threshold cites unknown "
                f"'{t.primary_citation}'")
    return problems


# cv_percent vs omega2 consistency tolerance: the stored convenience value is
# rounded for human display, so allow ~1 CV-percentage-point of slack.
_CV_TOL = 1.0


def _check_variability(m, known_citations) -> List[str]:
    """Variability-layer consistency checks (spec §9): cv<->omega2, citations
    resolve, and variability_status matches the curated contents."""
    problems: List[str] = []

    has_bsv = False
    for p in m.parameters:
        v = p.variability
        if v is None:
            continue
        if v.omega2 is not None:
            has_bsv = True
        # (1) cv_percent recomputes from omega2 within tolerance (Trap 1)
        if v.cv_percent is not None and v.cv_from_omega2 is not None:
            if abs(v.cv_percent - v.cv_from_omega2) > _CV_TOL:
                problems.append(
                    f"[variability] {m.id}: parameter {p.symbol} cv_percent "
                    f"{v.cv_percent:g} disagrees with omega2-derived "
                    f"{v.cv_from_omega2:.1f} (>{_CV_TOL} pp — variance/SD/CV confusion?)"
                )
        # (2) variability citations resolve
        if v.primary_citation and v.primary_citation not in known_citations:
            problems.append(
                f"[cite] {m.id}: parameter {p.symbol} variability cites unknown "
                f"'{v.primary_citation}'"
            )

    re_ = m.residual_error
    if re_ is not None and re_.primary_citation and re_.primary_citation not in known_citations:
        problems.append(f"[cite] {m.id}: residual_error cites unknown '{re_.primary_citation}'")

    ob = m.omega_block
    if ob is not None and ob.primary_citation and ob.primary_citation not in known_citations:
        problems.append(f"[cite] {m.id}: omega_block cites unknown '{ob.primary_citation}'")

    # (3) variability_status matches the actual contents
    status = m.variability_status
    if status == "none":
        if has_bsv or re_ is not None or ob is not None:
            problems.append(
                f"[variability] {m.id}: variability_status 'none' but curated "
                "BSV/residual/omega_block is present"
            )
    elif status == "full":
        if ob is None:
            problems.append(
                f"[variability] {m.id}: variability_status 'full' requires an omega_block"
            )
        if not has_bsv:
            problems.append(f"[variability] {m.id}: variability_status 'full' requires BSV")
    elif status in ("partial", "diagonal"):
        if not has_bsv:
            problems.append(
                f"[variability] {m.id}: variability_status '{status}' requires at least "
                "one parameter carrying BSV"
            )
        if status == "diagonal" and ob is not None and ob.complete:
            problems.append(
                f"[variability] {m.id}: variability_status 'diagonal' but a complete "
                "omega_block is present (should be 'full')"
            )

    return problems


# RSE recompute tolerance (percentage points) and CI<->SE relative tolerance.
_RSE_TOL = 1.0
_CI_REL_TOL = 0.05


def _check_estimation(m, known_citations) -> List[str]:
    """Estimation-uncertainty consistency checks (v0.3 spec §9): the SE/RSE/CI/scale
    traps (§4), citation resolution, and uncertainty_status ↔ contents.

    The RSE-vs-CV separation (Trap 1) is enforced *structurally* — estimation lives in
    its own block beside `variability`, so an RSE can never be silently filed as a BSV
    CV — not by a numeric check (both are plausible magnitudes; that is the human
    line item). These checks guard the numeric traps a machine *can* catch."""
    problems: List[str] = []
    has_estimation = False

    for p in m.parameters:
        e = p.estimation
        if e is None:
            continue
        if e.se is not None:
            has_estimation = True
            # (Trap 2) scale is mandatory when an SE is present
            if e.scale not in ("natural", "log"):
                problems.append(
                    f"[estimation] {m.id}: parameter {p.symbol} has an SE but no/invalid "
                    f"`scale` (got {e.scale!r}) — a log-scale SE read as natural is silently wrong"
                )
            # (Trap 2) rse_percent recomputes from se & central on the declared scale
            ref = e.rse_from_se(p.central)
            if e.rse_percent is not None and ref is not None:
                if abs(e.rse_percent - ref) > _RSE_TOL:
                    problems.append(
                        f"[estimation] {m.id}: parameter {p.symbol} rse_percent "
                        f"{e.rse_percent:g} disagrees with se-derived {ref:.1f} "
                        f"(>{_RSE_TOL} pp on scale '{e.scale}' — SE/RSE/scale confusion?)"
                    )
            # (Trap 3) ci95 consistent with se for a symmetric asymptotic interval
            if (e.method == "asymptotic_covariance" and e.scale == "natural"
                    and e.ci95 and e.ci95.get("low") is not None
                    and e.ci95.get("high") is not None):
                width = e.ci95["high"] - e.ci95["low"]
                expected = 2 * 1.96 * e.se
                if expected > 0 and abs(width - expected) > _CI_REL_TOL * expected:
                    problems.append(
                        f"[estimation] {m.id}: parameter {p.symbol} ci95 width {width:g} "
                        f"inconsistent with 2·1.96·se={expected:g} (asymptotic interval — Trap 3)"
                    )
        # estimation citation resolves
        if e.primary_citation and e.primary_citation not in known_citations:
            problems.append(
                f"[cite] {m.id}: parameter {p.symbol} estimation_uncertainty cites unknown "
                f"'{e.primary_citation}'"
            )

    ec = m.estimate_covariance
    if ec is not None and ec.primary_citation and ec.primary_citation not in known_citations:
        problems.append(f"[cite] {m.id}: estimate_covariance cites unknown '{ec.primary_citation}'")

    # uncertainty_status matches the actual contents (v0.3 §5/§9 item 4)
    status = m.uncertainty_status
    if status == "none":
        if has_estimation or ec is not None:
            problems.append(
                f"[estimation] {m.id}: uncertainty_status 'none' but curated estimation "
                "uncertainty / estimate_covariance is present"
            )
    elif status == "correlated":
        if ec is None:
            problems.append(
                f"[estimation] {m.id}: uncertainty_status 'correlated' requires an estimate_covariance"
            )
    elif status == "marginal":
        if not has_estimation:
            problems.append(
                f"[estimation] {m.id}: uncertainty_status 'marginal' requires at least one "
                "parameter carrying an estimation SE"
            )

    return problems


# Equation input variables that are universally available from a standard patient
# (so a model need not list them in covariates.required to bind an equation needing them).
_STANDARD_COVARIATES = {"weight", "height", "height_cm", "age", "sex", "bmi"}
# Map an equation input variable to the model-covariate name that supplies it.
_INPUT_TO_COVARIATE = {"height_cm": "height"}


def _check_covariate_model(m, ds, known_citations) -> List[str]:
    """Covariate-model layer checks (v0.7 spec §9):

    1. every ``derived_inputs[].equation`` resolves to a library record;
    2. every ``used_for`` symbol resolves to a real parameter;
    3. the equation's input covariates are available to the model (Trap 2 — the
       dimensional/availability half a machine can check; the cm-vs-m unit check is
       the human line item);
    4. ``covariate_sensitivity_status`` matches the curated contents;
    5. every covariate-layer ``primary_citation`` resolves.
    """
    problems: List[str] = []
    cm = m.covariate_model
    status = m.covariate_sensitivity_status

    # (4a) 'computed' is caller-side and must never appear in the dataset
    if status == "computed":
        problems.append(
            f"[covariate] {m.id}: covariate_sensitivity_status 'computed' is caller-side "
            "(a supplied covariate-value distribution) and must not appear in the dataset"
        )
    # (4b) declared <-> covariate_model presence
    if cm is None:
        if status == "declared":
            problems.append(
                f"[covariate] {m.id}: covariate_sensitivity_status 'declared' but no "
                "covariate_model block is present"
            )
        return problems
    if status != "declared":
        problems.append(
            f"[covariate] {m.id}: covariate_model present but covariate_sensitivity_status "
            f"is '{status}' (expected 'declared')"
        )

    param_symbols = {p.symbol for p in m.parameters}
    declared_covs = set(m.covariates.get("required", [])) | set(m.covariates.get("optional", []))
    eq_lib = getattr(ds, "covariate_equations", {})

    for d in cm.derived_inputs:
        # (1) equation resolves to a library record
        rec = eq_lib.get(d.equation)
        if rec is None:
            problems.append(
                f"[covariate] {m.id}: derived input '{d.quantity}' binds equation "
                f"'{d.equation}' which is not in covariate_equations/"
            )
        else:
            # (1b) the equation's quantity matches the binding's declared quantity
            if rec.get("quantity") and rec.get("quantity") != d.quantity:
                problems.append(
                    f"[covariate] {m.id}: derived input declares quantity '{d.quantity}' but "
                    f"equation '{d.equation}' computes '{rec.get('quantity')}'"
                )
            # (3) the equation's inputs are available to the model
            for var in rec.get("inputs", []):
                cov = _INPUT_TO_COVARIATE.get(var, var)
                if cov not in _STANDARD_COVARIATES and cov not in declared_covs:
                    problems.append(
                        f"[covariate] {m.id}: equation '{d.equation}' needs covariate '{var}' "
                        f"but the model does not declare it (covariates.required/optional)"
                    )
        # (2) used_for symbols resolve to real parameters
        for sym in d.used_for:
            if sym not in param_symbols:
                problems.append(
                    f"[covariate] {m.id}: derived input '{d.equation}' used_for '{sym}' "
                    "is not a parameter of this model"
                )
        # (5) covariate-layer citation resolves
        if d.primary_citation and d.primary_citation not in known_citations:
            problems.append(
                f"[cite] {m.id}: covariate_model '{d.equation}' cites unknown "
                f"'{d.primary_citation}'"
            )
    return problems


def _check_covariate_equations(ds, known_citations) -> List[str]:
    """Covariate-equation library checks (v0.7 spec §9): each record's citation and
    failure-mode citations resolve, its validity envelope is well-ordered, and it is
    actually implemented in the pure equation registry (data/code stay in sync — an
    equation a model could bind but the code cannot evaluate is a silent trap)."""
    from .covariates import EQUATIONS  # local import to avoid a cycle

    problems: List[str] = []
    for eid, rec in getattr(ds, "covariate_equations", {}).items():
        cid = rec.get("primary_citation")
        if cid and cid not in known_citations:
            problems.append(f"[cite] covariate_equation {eid}: primary_citation '{cid}' not in citations/")
        for fm in rec.get("known_failure_modes", []):
            fcid = fm.get("citation")
            if fcid and fcid not in known_citations:
                problems.append(f"[cite] covariate_equation {eid}: failure-mode cites unknown '{fcid}'")
        env = rec.get("validity_envelope") or {}
        for axis, rng in env.items():
            if isinstance(rng, dict):
                lo, hi = rng.get("min"), rng.get("max")
                if lo is not None and hi is not None and lo > hi:
                    problems.append(f"[covariate] equation {eid}: {axis} min {lo} > max {hi}")
        if eid not in EQUATIONS:
            problems.append(
                f"[covariate] equation {eid}: curated in the library but not implemented in "
                "hypnos.covariates.EQUATIONS (data/code out of sync)"
            )
    return problems


def _check_review_state(m) -> List[str]:
    """Review-state integrity — the ``pending_human_review`` governance line (the new state).

    An automated source cross-check may *populate evidence* (set ``pending_human_review`` +
    a ``source_review`` provenance block) but NEVER asserts human verification:
    ``human_verified`` is always false outside ``verified``. A ``pending_human_review``
    record must carry the provenance (what was compared, against which fetched source(s));
    a ``contested`` record may carry it to document a machine-found disagreement."""
    problems: List[str] = []
    sr = m.source_review
    if m.review_status == "pending_human_review":
        if sr is None:
            problems.append(
                f"[review] {m.id}: review_status 'pending_human_review' requires a "
                "source_review provenance block (what was cross-checked, against which source)")
        elif sr.get("human_verified") is not False:
            problems.append(
                f"[review] {m.id}: a pending_human_review source_review must have "
                "human_verified=false (an automated check never reaches human verification)")
        elif not sr.get("sources"):
            problems.append(
                f"[review] {m.id}: pending_human_review source_review must cite the source(s) "
                "actually compared against (sources[] is empty)")
    elif m.review_status == "contested":
        if sr is not None and sr.get("human_verified") is not False:
            problems.append(
                f"[review] {m.id}: a contested source_review must have human_verified=false")
    elif sr is not None and m.review_status != "verified":
        problems.append(
            f"[review] {m.id}: a source_review block is present but review_status is "
            f"'{m.review_status}' (expected 'pending_human_review', 'contested', or 'verified')")
    return problems


def _check_developmental(m, known_citations) -> List[str]:
    """Developmental-layer checks (v0.8 spec §8):

    * the block and its size/maturation components cite resolvably;
    * ``evidence_tier`` is D (an extrapolation defaults to D; a fitted-in-children model
      is a different record, labeled ``fitted_pediatric``);
    * ``applied_by_default`` is never true (an extrapolation is opt-in only);
    * a ``maturation`` block's ``driver`` is a PMA-class clock (Trap 1);
    * ``extrapolation_basis`` matches the curated contents (allometry_plus_maturation
      requires a maturation block; allometry_only forbids one).
    """
    problems: List[str] = []
    dev = m.developmental_model
    if dev is None:
        return problems
    if dev.primary_citation and dev.primary_citation not in known_citations:
        problems.append(f"[cite] {m.id}: developmental_model cites unknown '{dev.primary_citation}'")
    if dev.size and dev.size.primary_citation and dev.size.primary_citation not in known_citations:
        problems.append(f"[cite] {m.id}: developmental size cites unknown '{dev.size.primary_citation}'")
    if dev.maturation and dev.maturation.primary_citation \
            and dev.maturation.primary_citation not in known_citations:
        problems.append(
            f"[cite] {m.id}: developmental maturation cites unknown '{dev.maturation.primary_citation}'")
    if dev.extrapolation_basis != "fitted_pediatric" and dev.evidence_tier != "D":
        problems.append(
            f"[developmental] {m.id}: evidence_tier '{dev.evidence_tier}' for an extrapolation "
            "must be D (a fitted-in-children model is a different record; v0.8 §5)")
    if dev.applied_by_default:
        problems.append(
            f"[developmental] {m.id}: applied_by_default must be false — an allometric/maturation "
            "extrapolation is opt-in only, never a silent reparameterization (v0.8 §9)")
    if dev.maturation is not None and dev.maturation.driver != "pma_weeks":
        problems.append(
            f"[developmental] {m.id}: maturation driver '{dev.maturation.driver}' must be a PMA-class "
            "clock (pma_weeks), never chronological age (v0.8 §4 Trap 1)")
    if dev.extrapolation_basis == "allometry_plus_maturation" and dev.maturation is None:
        problems.append(
            f"[developmental] {m.id}: extrapolation_basis 'allometry_plus_maturation' but no "
            "maturation block is curated")
    if dev.extrapolation_basis == "allometry_only" and dev.maturation is not None:
        problems.append(
            f"[developmental] {m.id}: extrapolation_basis 'allometry_only' but a maturation block "
            "is present (should be 'allometry_plus_maturation')")
    return problems


def _check_pharmacogenomics(m, known_citations) -> List[str]:
    """Pharmacogenomic-layer checks (v0.9 spec §8), enforcing the flag/modifier separation:

    * a kinetic modifier is Tier-D, opt-in (``applied_by_default`` false), carries a
      non-empty ``substrate_scope`` that includes this model's drug (the substrate
      guardrail — a genetic effect never leaks to a non-substrate drug, Trap 4);
    * a safety flag carries ``safety_critical`` and resolvable citation;
    * every modifier/flag ``primary_citation`` resolves.

    The category-error guard (Trap 2 — a susceptibility is not a kinetic scale-factor) is
    enforced *structurally* by the schema: the two ``$defs`` are disjoint (a modifier
    cannot carry ``trigger_agents``, a flag cannot carry an ``adjustment``), so a flag has
    nothing to round-trip kinetically.
    """
    problems: List[str] = []
    for mod in m.pharmacogenomic_modifiers:
        if mod.evidence_tier != "D":
            problems.append(
                f"[pgx] {m.id}: kinetic modifier ({mod.gene}) evidence_tier '{mod.evidence_tier}' "
                "must be D unless backed by a linked fitted-by-genotype model (v0.9 §5)")
        if mod.applied_by_default:
            problems.append(
                f"[pgx] {m.id}: kinetic modifier ({mod.gene}) applied_by_default must be false "
                "(opt-in only; v0.9 §9)")
        if not mod.substrate_scope:
            problems.append(
                f"[pgx] {m.id}: kinetic modifier ({mod.gene}) needs a non-empty substrate_scope "
                "(Trap 4 — a mis-scoped genetic effect is silently wrong)")
        elif m.drug_name not in mod.substrate_scope:
            problems.append(
                f"[pgx] {m.id}: kinetic modifier ({mod.gene}) attaches to drug '{m.drug_name}' "
                f"which is not in its substrate_scope {mod.substrate_scope} (substrate guardrail, Trap 4)")
        if mod.primary_citation and mod.primary_citation not in known_citations:
            problems.append(f"[cite] {m.id}: pgx modifier cites unknown '{mod.primary_citation}'")
    for flag in m.pharmacogenomic_safety_flags:
        if not flag.safety_critical:
            problems.append(
                f"[pgx] {m.id}: safety flag ({flag.gene}) must carry safety_critical=true (v0.9 §4)")
        if flag.primary_citation and flag.primary_citation not in known_citations:
            problems.append(f"[cite] {m.id}: pgx safety flag cites unknown '{flag.primary_citation}'")
    return problems


# mode -> the targets it may carry (the Varvel quantity must match the modality)
_MODE_TARGETS = {
    "pk_concentration": {"cp", "ce"},
    "pd_bis": {"bis"},
    "pd_tof": {"tof"},
}


# In-envelope MDAPE band ceilings per declared tier (v0.1 §5, the inaccuracy thresholds).
# A computed in-envelope MDAPE above the ceiling FALSIFIES the tier claim (advisory only).
_TIER_MDAPE_CEILING = {"A": 30.0, "B": 40.0, "C": 50.0}


def _check_tier_falsification(m) -> List[str]:
    """Make the tier a falsifiable claim (v0.4 §5): if a model's Hypnos-computed **in-envelope**
    MDAPE exceeds its declared tier's band, raise an ADVISORY tier-mismatch flag — not an
    auto-demotion. A computed metric *informs* but never *sets* a tier (humans still gate
    promotion). Dormant until a model carries an in-envelope external_validation MDAPE."""
    problems: List[str] = []
    ceiling = _TIER_MDAPE_CEILING.get(m.tier)
    if ceiling is None:
        return problems
    for e in m.external_validation:
        if (e.get("cohort") or {}).get("in_envelope") is not True:
            continue  # only an IN-envelope metric falsifies the tier (out-of-envelope is expected to be worse)
        for met in e.get("metrics", []):
            if met.get("name") == "MDAPE" and met.get("value") is not None and met["value"] > ceiling:
                problems.append(
                    f"[tier] {m.id}: declared Tier {m.tier} but Hypnos-computed IN-ENVELOPE MDAPE "
                    f"{met['value']:g}% exceeds the Tier-{m.tier} band (~{ceiling:g}%) "
                    f"[dataset {e.get('dataset')}] — ADVISORY tier-mismatch flag for human review "
                    "(a computed metric informs but never sets a tier; v0.4 §5)")
    return problems


def _check_external_validation(m) -> List[str]:
    """External-validation consistency checks (v0.4 spec §4): validation_status
    matches the curated entries, and each entry's mode/target/CI are coherent.

    These are dormant until a model carries computed metrics; they enforce the
    block's invariants the moment one does, so a Hypnos-computed artifact can never
    be mislabeled (e.g. a BIS validation filed as a concentration validation)."""
    problems: List[str] = []
    entries = m.external_validation
    status = m.validation_status

    if entries and status == "none":
        problems.append(
            f"[validation] {m.id}: external_validation present but validation_status 'none'"
        )

    has_pk = any(e.get("mode") == "pk_concentration" for e in entries)
    has_pd = any(e.get("mode") in ("pd_bis", "pd_tof") for e in entries)
    if status in ("external_pk", "external_both") and not has_pk:
        problems.append(
            f"[validation] {m.id}: validation_status '{status}' requires a "
            "pk_concentration external_validation entry"
        )
    if status in ("external_pd", "external_both") and not has_pd:
        problems.append(
            f"[validation] {m.id}: validation_status '{status}' requires a "
            "pd_bis/pd_tof external_validation entry"
        )

    for i, e in enumerate(entries):
        mode = e.get("mode")
        target = e.get("target")
        allowed = _MODE_TARGETS.get(mode, set())
        if target not in allowed:
            problems.append(
                f"[validation] {m.id}: external_validation[{i}] target '{target}' "
                f"inconsistent with mode '{mode}' (expected one of {sorted(allowed)})"
            )
        for met in e.get("metrics", []):
            ci = met.get("ci95")
            if ci:
                lo, hi = ci.get("low"), ci.get("high")
                if lo is not None and hi is not None and lo > hi:
                    problems.append(
                        f"[validation] {m.id}: external_validation[{i}] metric "
                        f"{met.get('name')} ci95 low {lo} > high {hi}"
                    )
    return problems


def assert_valid(ds: Optional[Dataset] = None) -> None:
    problems = validate_dataset(ds)
    if problems:
        raise ValidationError(
            f"{len(problems)} dataset problem(s):\n  " + "\n  ".join(problems)
        )
