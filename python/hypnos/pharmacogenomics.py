"""The pharmacogenomic-envelope overlay — the v0.9 genetic special-population axis.

Two record kinds, never blurred (v0.9 §2):

* a **kinetic modifier** scales a forward prediction (opt-in, Tier-D, substrate-scoped);
* a **safety flag** forbids a trigger (MH) or warns of a prolonged effect (atypical
  BCHE) — and carries **no** kinetic effect.

:func:`pgx_overlay` answers, for a *declared* genotype, "which drugs carry a cited
kinetic modifier, which carry a hard safety/avoidance flag, and which — the honest
majority — carry neither because no actionable evidence exists" (v0.9 §6). It leads
with avoidance, because that is the true register of anesthetic pharmacogenetics.

Two non-negotiable rules (v0.9 §5):

* **never-invent** — no modifier/flag without a curated, cited record;
* **never-infer** — a modifier or flag fires ONLY on a genotype/phenotype the caller
  explicitly declares; nothing is assumed from ancestry, frequency, or likelihood.
  Absent a declaration there is no genetic effect, shown as "genotype not declared",
  never "wild-type assumed". This is an ethics line as much as a scientific one.

There is no inverse control: a genotype shifts a *forward* prediction or raises a
contraindication; it never yields a dose (v0.9 §9).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .load import Dataset, load

# The declared genotype dimensions Hypnos understands. A patient/genotype dict carries a
# subset of these; anything else is ignored (never inferred).
GENOTYPE_DIMENSIONS = ("bche_phenotype", "cyp2d6_phenotype", "cyp3a_phenotype",
                       "mh_susceptibility")

# Atypical BCHE phenotypes that trigger the prolonged-block flag/modifier (the normal
# phenotype triggers nothing). Heterozygous is milder than homozygous, but both are
# clinically relevant (v0.9 §3).
_ATYPICAL_BCHE = ("heterozygous_atypical", "homozygous_atypical")

# Substrate-specificity guardrails (v0.9 §3 — the anti-trap facts). These are
# DEFINITIONAL metabolic-route facts (which enzyme hydrolyzes/metabolizes a drug), not
# fitted PK, so they live in code, named — the same stance as the organ-staging
# thresholds in models.py. They drive the "explicitly NOT affected" message that keeps a
# genetic effect from leaking across a drug class.
_SUBSTRATE_GUARDRAILS = {
    "bche_phenotype": [
        ("remifentanil", "hydrolyzed by nonspecific blood/tissue esterases, NOT by BCHE — "
                         "atypical BCHE does not prolong remifentanil"),
        ("lidocaine", "amide local anesthetic — hepatically metabolized, not a BCHE substrate"),
        ("bupivacaine", "amide local anesthetic — hepatically metabolized, not a BCHE substrate"),
        ("levobupivacaine", "amide local anesthetic — hepatically metabolized, not a BCHE substrate"),
        ("ropivacaine", "amide local anesthetic — hepatically metabolized, not a BCHE substrate"),
        ("rocuronium", "aminosteroid neuromuscular blocker — not a cholinesterase substrate "
                       "(sugammadex is its reversal path)"),
    ],
}


def _declared(genotype: Dict[str, Any], dimension: str) -> Optional[Any]:
    """The declared value for a genotype dimension, or None when not declared."""
    if genotype is None:
        return None
    v = genotype.get(dimension)
    return v if (v is not None and v != "" and v != "normal") else None


def genotype_triggers(genotype: Dict[str, Any], dimension: str, value: Any) -> bool:
    """Does the declared genotype trigger a record curated for (dimension, value)?

    Never-infer: only an explicitly declared dimension can trigger. For BCHE, any
    declared atypical phenotype triggers an atypical-curated record (heterozygous and
    homozygous both matter); for MH, a truthy declaration triggers. Other dimensions
    trigger on exact phenotype match."""
    declared = _declared(genotype, dimension)
    if declared is None:
        return False
    if dimension == "mh_susceptibility":
        return bool(declared)
    if dimension == "bche_phenotype":
        return str(declared) in _ATYPICAL_BCHE and str(value) in _ATYPICAL_BCHE
    return str(declared) == str(value)


@dataclass
class PgxFinding:
    """One pharmacogenomic record matched to a declared genotype, with its drug context."""

    drug: str
    gene: str
    model_id: Optional[str]
    record: Any                     # PharmacogenomicModifier | PharmacogenomicSafetyFlag


@dataclass
class PgxOverlay:
    """The declared-genotype overlay (v0.9 §6). Avoidance first, then kinetics, then the gaps."""

    genotype: Dict[str, Any]
    safety_flags: List[PgxFinding] = field(default_factory=list)
    modifiers: List[PgxFinding] = field(default_factory=list)
    substrate_guardrails: List[Dict[str, str]] = field(default_factory=list)
    no_evidence: List[str] = field(default_factory=list)

    @property
    def declared_dimensions(self) -> List[str]:
        return [d for d in GENOTYPE_DIMENSIONS if _declared(self.genotype, d) is not None]

    @property
    def avoidance_flags(self) -> List[PgxFinding]:
        return [f for f in self.safety_flags if getattr(f.record, "kind", "avoidance") == "avoidance"]

    @property
    def awareness_flags(self) -> List[PgxFinding]:
        return [f for f in self.safety_flags if getattr(f.record, "kind", "avoidance") == "awareness"]


def pgx_overlay(ds: Optional[Dataset] = None, *, genotype: Dict[str, Any]) -> PgxOverlay:
    """Build the declared-genotype overlay (v0.9 §6).

    Scans every curated model for pharmacogenomic safety flags and kinetic modifiers
    whose phenotype the ``genotype`` declares, leads with the avoidance flags, lists the
    opt-in Tier-D kinetic modifiers, surfaces the substrate-specificity guardrails (the
    explicitly NOT-affected drugs), and names the honest majority with no actionable
    evidence. Nothing is inferred — only declared dimensions can match (v0.9 §5)."""
    if ds is None:
        ds = load()

    overlay = PgxOverlay(genotype=dict(genotype or {}))
    matched_drugs: set = set()

    for m in ds:
        drug = m.drug_name
        for flag in m.pharmacogenomic_safety_flags:
            if genotype_triggers(genotype, flag.phenotype_dimension, flag.phenotype_value):
                overlay.safety_flags.append(PgxFinding(drug, flag.gene, m.id, flag))
                matched_drugs.add(drug)
        for mod in m.pharmacogenomic_modifiers:
            if genotype_triggers(genotype, mod.phenotype_dimension, mod.phenotype_value):
                overlay.modifiers.append(PgxFinding(drug, mod.gene, m.id, mod))
                matched_drugs.add(drug)

    # Substrate guardrails — the explicitly NOT-affected drugs, for each declared gene that
    # carries a guardrail table (today: BCHE). Only surface a guardrail for a drug actually
    # in the dataset, so the message reflects this formulary.
    dataset_drugs = {d.get("name") for d in ds.drugs.values()}
    guardrail_drugs: set = set()
    for dim, table in _SUBSTRATE_GUARDRAILS.items():
        if _declared(genotype, dim) is None:
            continue
        gene = "BCHE" if dim == "bche_phenotype" else dim
        for drug, reason in table:
            if drug in dataset_drugs:
                overlay.substrate_guardrails.append({"gene": gene, "drug": drug, "reason": reason})
                guardrail_drugs.add(drug)

    # The honest majority: dataset drugs with neither a matched record nor a guardrail note.
    overlay.no_evidence = sorted(
        d for d in dataset_drugs
        if d and d not in matched_drugs and d not in guardrail_drugs
    )
    return overlay
