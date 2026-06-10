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


def _load_schema(root: Path) -> dict:
    schema_path = root / "schema" / "model.schema.json"
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
    except ImportError:  # pragma: no cover - jsonschema is a hard dependency
        problems.append("[schema] jsonschema not installed; skipped schema validation")

    # --- integrity layer --------------------------------------------------
    from .export.registry import (  # local import to avoid cycle
        INTERACTION_KERNELS,
        KERNELS,
        PD_KERNELS,
        VOLATILE_KERNELS,
    )

    known_kernels = set(KERNELS) | set(PD_KERNELS) | set(INTERACTION_KERNELS) | set(VOLATILE_KERNELS)

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
                     "crcl_ml_min", "albumin_g_dl", "ejection_fraction_pct"):
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

        # id prefix matches declared subsystem
        if m.id.split(".")[0] != m.subsystem:
            problems.append(f"[id] {m.id}: id prefix does not match subsystem '{m.subsystem}'")

        # --- variability layer (v0.2 spec §9) -----------------------------
        problems.extend(_check_variability(m, known_citations))

        # --- external-validation layer (v0.4 spec §4) ---------------------
        problems.extend(_check_external_validation(m))

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


# mode -> the targets it may carry (the Varvel quantity must match the modality)
_MODE_TARGETS = {
    "pk_concentration": {"cp", "ce"},
    "pd_bis": {"bis"},
    "pd_tof": {"tof"},
}


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
