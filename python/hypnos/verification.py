"""Verification workflow support — tooling for the single highest-leverage
contribution (spec §9): promoting ``unverified`` models to ``verified`` by
reading the source PDF, field by field.

This module **guides** human verification; it never performs it. Promotion to
``verified`` requires a human to confirm the parameters *and the covariate
equations* against the primary source and then edit the record's
``extraction.review_status`` (filling ``verified_by`` / ``verified_date``).
Nothing here writes that field — by design (``Humans verify; LLMs do not
promote``).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .load import Dataset, load
from .models import TIER_RANK, Model

REVIEW_STATES = ("unverified", "verified", "contested")


@dataclass
class ChecklistItem:
    group: str           # "structural" | "covariate" | "envelope" | "population" | "estimation" | "citation"
    label: str
    value: str           # the current value a verifier must confirm against the PDF
    confirmed: bool = False
    # curated provenance pointer (e.g. "Schnider 1998, Table 2") so a verifier goes
    # straight to the right table; None where no source_locator is curated yet — which
    # the renderers flag as a gap to fill, the same honesty the dataset applies to itself.
    locator: Optional[str] = None


@dataclass
class ModelVerification:
    model_id: str
    label: str
    tier: str
    review_status: str
    kernel_implemented: bool
    primary_citation: str
    doi: Optional[str]
    pmid: Optional[str]
    source_locator: Optional[str]
    checklist: List[ChecklistItem] = field(default_factory=list)
    blocking: List[str] = field(default_factory=list)   # what stands between this and `verified`

    @property
    def n_items(self) -> int:
        return len(self.checklist)


def _checklist_for(model: Model, ds: Dataset) -> List[ChecklistItem]:
    items: List[ChecklistItem] = []
    # 1. structural parameters
    for p in model.parameters:
        units = p.units or ""
        central = "covariate" if p.central is None else f"{p.central:g}"
        loc = (p.extraction or {}).get("source_locator")
        items.append(ChecklistItem("structural", f"parameter {p.symbol}",
                                   f"{central} {units}".strip(), locator=loc))
    # 2. covariate equations (where transcription errors hide)
    for p in model.parameters:
        if p.covariate_model:
            items.append(ChecklistItem("covariate", f"{p.symbol} equation", p.covariate_model))
    cov = model.covariates
    if cov.get("lbm_equation"):
        items.append(ChecklistItem("covariate", "LBM/FFM equation", cov["lbm_equation"]))
    # 3. derivation population and n
    env = model.applicability_envelope
    if env.populations:
        items.append(ChecklistItem("population", "derivation populations", ", ".join(env.populations)))
    if env.derivation_n is not None:
        items.append(ChecklistItem("population", "derivation n", str(env.derivation_n)))
    # 4. applicability range
    for name in ("age_years", "weight_kg", "height_cm", "bmi_kg_m2"):
        rng = getattr(env, name)
        if rng.min is not None or rng.max is not None:
            items.append(ChecklistItem("envelope", name, f"[{rng.min}, {rng.max}]"))
    # 5. estimation uncertainty (v0.3 §9) — the RSE-vs-CV disambiguation is the
    #    cardinal human line item: confirm each curated number is an estimation
    #    RSE/SE, not a between-subject CV, against the source column header.
    for p in model.parameters:
        e = p.estimation
        if e is None:
            continue
        loc = (e.extraction or {}).get("source_locator")
        bits = []
        if e.se is not None:
            bits.append(f"se={e.se:g} (scale={e.scale})")
        if e.rse_percent is not None:
            bits.append(f"RSE%={e.rse_percent:g}")
        bits.append(f"method={e.method}")
        items.append(ChecklistItem(
            "estimation", f"{p.symbol} estimation uncertainty (RSE/SE — NOT a BSV CV)",
            "; ".join(bits), locator=loc))
    if model.estimate_covariance is not None:
        ec = model.estimate_covariance
        items.append(ChecklistItem(
            "estimation", "estimate covariance ($COV)",
            f"{len(ec.correlations)} pair(s), complete={ec.complete}, "
            f"cov_step_succeeded={ec.covariance_step_succeeded}"))
    # 5b. local-anesthetic toxicity thresholds (v0.6 LA1 §8) — the safety-critical
    #     group; given the stakes, threshold RANGES specifically require human source
    #     confirmation before `verified`. The five checks the spec enumerates:
    drug = ds.drug(model.drug_name) if ds is not None else None
    saturable = bool(((drug or {}).get("protein_binding") or {}).get("saturable"))
    for th in model.toxicity_thresholds:
        loc = (th.extraction or {}).get("source_locator")
        # (1) basis total vs free, and conversion correct (Trap 1)
        items.append(ChecklistItem(
            "local_anesthetic", f"{th.endpoint} threshold BASIS (total vs free — Trap 1)",
            f"[{th.low:g}, {th.high:g}] {th.units} on {th.basis}", locator=loc))
        # (4) saturation caveat present for a binding-sensitive (saturable) agent
        if saturable and th.basis == "total_plasma":
            items.append(ChecklistItem(
                "local_anesthetic", f"{th.endpoint} saturation caveat (free fraction rises non-linearly)",
                (th.saturation_caveat or "⚠️ MISSING — required for a saturable drug on total basis")))
        # (3) speed-of-rise / method context recorded (Trap 3)
        items.append(ChecklistItem(
            "local_anesthetic", f"{th.endpoint} method/speed-of-rise context (Trap 3)",
            th.method_caveat or "⚠️ no method_caveat curated"))
    if model.has_toxicity_thresholds:
        # (2) salt vs base and units (Trap 2); (5) curated as a RANGE, not a line
        items.append(ChecklistItem(
            "local_anesthetic", "salt-vs-base & concentration units (Trap 2)",
            "confirm assay basis (base) vs dose salt form against the source"))
        items.append(ChecklistItem(
            "local_anesthetic", "every threshold curated as a RANGE, never a line (§3.3)",
            f"{len(model.toxicity_thresholds)} threshold range(s); none collapsed to a single value"))
    # 6. citation resolves to the right paper
    cit = ds.citation(model.primary_citation) if ds is not None else None
    if cit:
        ref = f"{cit.get('container', '')} {cit.get('year', '')}; doi:{cit.get('doi', '?')}"
        items.append(ChecklistItem("citation", "primary citation resolves", ref.strip()))
    return items


def model_verification(ds: Dataset, model_id: str) -> ModelVerification:
    m = ds[model_id]
    cit = ds.citation(m.primary_citation) or {}
    blocking: List[str] = []
    if m.review_status == "verified":
        pass
    else:
        blocking.append("a human must confirm every field below against the source PDF")
        if any(p.central is None for p in m.parameters):
            blocking.append("some parameters have no central value (transcription incomplete)")
        if not m.kernel_implemented and m.purpose in ("pk", "pd", "interaction", "physicochemical"):
            blocking.append("reference kernel pending — verify parameters before implementing it")
    return ModelVerification(
        model_id=m.id, label=m.label, tier=m.tier, review_status=m.review_status,
        kernel_implemented=m.kernel_implemented, primary_citation=m.primary_citation,
        doi=cit.get("doi"), pmid=cit.get("pmid"),
        source_locator=m.extraction.get("source_locator"),
        checklist=_checklist_for(m, ds), blocking=blocking,
    )


def _priority(m: Model) -> tuple:
    # Highest leverage first: implemented kernels (verifying unlocks trustworthy
    # simulation), then better tier, then by id for stability.
    return (0 if m.kernel_implemented else 1, TIER_RANK.get(m.tier, 9), m.id)


def next_to_verify(ds: Dataset, limit: Optional[int] = None) -> List[Model]:
    """Unverified models in highest-leverage order (implemented kernel, best tier first)."""
    pending = [m for m in ds if m.review_status != "verified"]
    pending.sort(key=_priority)
    return pending[:limit] if limit else pending


def verification_summary(ds: Optional[Dataset] = None) -> Dict[str, Any]:
    if ds is None:
        ds = load()
    by_status: Dict[str, int] = {s: 0 for s in REVIEW_STATES}
    for m in ds:
        by_status[m.review_status] = by_status.get(m.review_status, 0) + 1
    n = len(ds)
    verified = by_status.get("verified", 0)
    return {
        "n_models": n,
        "by_review_status": by_status,
        "verified_fraction": (verified / n) if n else 0.0,
        "next_to_verify": [
            {"model_id": m.id, "tier": m.tier, "kernel": m.kernel_implemented,
             "citation": m.primary_citation}
            for m in next_to_verify(ds, limit=5)
        ],
    }


def checklist_markdown(mv: ModelVerification) -> str:
    """A copy-pasteable verification checklist (e.g. for a PR description)."""
    lines = [
        f"# Verification checklist — `{mv.model_id}`",
        "",
        f"- **Model:** {mv.label}",
        f"- **Tier:** {mv.tier}  ·  **Status:** {mv.review_status}  ·  "
        f"**Kernel:** {'implemented' if mv.kernel_implemented else 'pending'}",
        f"- **Source:** doi:{mv.doi or '?'}  ·  PMID:{mv.pmid or '?'}"
        + (f"  ·  {mv.source_locator}" if mv.source_locator else ""),
        "",
        "Confirm every item below **against the source PDF**, then set "
        "`extraction.review_status = \"verified\"` and fill `verified_by` / "
        "`verified_date`. The covariate equations are the part most worth "
        "double-checking; that is where published-vs-implemented divergence hides.",
        "",
    ]
    groups = ["structural", "covariate", "population", "envelope", "estimation",
              "local_anesthetic", "citation"]
    titles = {
        "structural": "Structural parameters",
        "covariate": "Covariate equations (incl. exact LBM/FFM form)",
        "population": "Derivation population & n",
        "envelope": "Stated applicability range",
        "estimation": "Estimation uncertainty (RSE/SE vs BSV CV — read the column header)",
        "local_anesthetic": "LA toxicity thresholds (basis, units, method, saturation, range-not-line — safety-critical)",
        "citation": "Primary citation",
    }
    for g in groups:
        items = [it for it in mv.checklist if it.group == g]
        if not items:
            continue
        lines.append(f"## {titles[g]}")
        for it in items:
            where = (f"  — _{it.locator}_" if it.locator
                     else "  — ⚠️ **no source locator curated; add one**" if it.group == "structural"
                     else "")
            lines.append(f"- [ ] **{it.label}** = `{it.value}`{where}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
