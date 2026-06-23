"""Provenance & safety annotations attached to every export.

Two durability tiers, mirroring the spec (§7):

* **MIRIAM-style RDF** (``bqmodel:isDerivedFrom`` -> DOI/PMID) survives even if a
  downstream tool strips the custom ``hypnos:`` predicates;
* **custom ``hypnos:`` predicates** carry the propagated confidence tier, the
  dataset version, and the universal, machine-readable
  ``hypnos:clinicalUse = "PROHIBITED"`` flag.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import CLINICAL_USE, __version__
from ..models import Model

HYPNOS_NS = "https://w3id.org/hypnos/terms#"


def citation_uris(model: Model, ds) -> List[str]:
    """Return resolvable DOI/PMID URIs for a model's primary citation."""
    uris: List[str] = []
    cit = ds.citation(model.primary_citation) if ds is not None else None
    if cit:
        if cit.get("doi"):
            uris.append(f"https://doi.org/{cit['doi']}")
        if cit.get("pmid"):
            uris.append(f"https://identifiers.org/pubmed:{cit['pmid']}")
    return uris


def provenance(model: Model, ds=None, tier: Optional[str] = None) -> Dict[str, Any]:
    """Structured provenance block reused across all exporters."""
    prov = {
        "hypnos:datasetVersion": __version__,
        "hypnos:modelId": model.id,
        "hypnos:confidenceTier": tier or model.tier,
        "hypnos:reviewStatus": model.review_status,
        "hypnos:clinicalUse": CLINICAL_USE,
        "bqmodel:isDerivedFrom": citation_uris(model, ds),
    }
    # Local-anesthetic records are doubly safety-critical (v0.6 §7/§9): the
    # concentration-vs-threshold view is exactly the artifact that tempts the
    # forbidden "max safe dose?" question, so a downstream consumer is twice-warned.
    if model.is_safety_critical:
        prov["hypnos:safetyCritical"] = "true"
    return prov


def banner(model: Model, tier: Optional[str] = None) -> str:
    """Human-readable banner for text-format exports (NONMEM, R, etc.)."""
    t = tier or model.tier
    safety = (
        "  SAFETY-CRITICAL (local anesthetic): systemic concentration only — NOT a\n"
        "  toxicity margin, max-dose, or 'is this safe?' tool. Thresholds are RANGES.\n"
        if model.is_safety_critical else ""
    )
    return (
        "============================================================\n"
        "  HYPNOS EXPORT — NOT FOR CLINICAL USE\n"
        f"  clinicalUse = {CLINICAL_USE}\n"
        f"  model: {model.id}  (tier {t}, review: {model.review_status})\n"
        f"  dataset version: {__version__}\n"
        + safety
        + "  Research / education / simulation only. Not a dosing tool.\n"
        "============================================================"
    )


def variability_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates describing the curated random-effects layer (v0.2 §8).

    SBML core cannot express population random effects, so the Ω/Σ ride as annotation
    metadata: a deterministic SBML consumer (COPASI/Tellurium) still sees only the
    typical patient (the parameters + rate rules), but the curated NLME object travels
    with the model and survives a round-trip. Returns the inner predicate lines (to be
    embedded inside the model's single ``rdf:RDF`` block) — or ``""`` when the model
    publishes no variability (the never-synthesize rule: absence is a true gap).
    """
    from ._variability import omega_correlations, omega_diagonal, residual_spec

    if not model.has_published_variability:
        return ""
    pad = f"{indent}    "
    lines = [
        f"{pad}<hypnos:variabilityStatus>{model.variability_status}</hypnos:variabilityStatus>",
        f"{pad}<hypnos:bandTier>{model.band_tier}</hypnos:bandTier>",
        f"{pad}<hypnos:betweenSubjectVariability>",
        f"{pad}  <rdf:Bag>",
    ]
    for sym, om2, cv in omega_diagonal(model):
        lines.append(f'{pad}    <rdf:li hypnos:parameter="{sym}" '
                     f'hypnos:omega2="{om2:.10g}" hypnos:cvPercent="{cv:.4g}"/>')
    lines += [f"{pad}  </rdf:Bag>", f"{pad}</hypnos:betweenSubjectVariability>"]
    for a, b, r in omega_correlations(model):
        lines.append(f'{pad}<hypnos:omegaCorrelation hypnos:between="{a} {b}" '
                     f'hypnos:correlation="{r:.6g}"/>')
    spec = residual_spec(model)
    if spec is not None:
        attrs = f'hypnos:model="{spec.model}"'
        if spec.log_sd is not None:
            attrs += f' hypnos:logSd="{spec.log_sd:.10g}"'
        if spec.prop_var is not None:
            attrs += f' hypnos:proportionalVariance="{spec.prop_var:.10g}"'
        if spec.add_sd is not None:
            attrs += f' hypnos:additiveSd="{spec.add_sd:.10g}"'
        lines.append(f"{pad}<hypnos:residualError {attrs}/>")
    return "\n".join(lines)


def la_rdf(model: Model, drug: Optional[Dict[str, Any]] = None, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the local-anesthetic blocks (v0.6 §9).

    Carries the site-absorption rank order, the toxicity-threshold *ranges*, and (v0.6
    LA2) the drug's ``cardiotoxicity_class`` as metadata so they travel with an
    SBML/PharmML model whose core is deterministic. Threshold ranges export **as
    ranges** — there is no format projection that collapses them to a single value
    (v0.6 §9). ``drug`` is the resolved drug record (for the cardiotoxicity class);
    pass ``None`` to omit it. Returns ``""`` for a non-LA model (and emits nothing for
    blocks the model does not carry)."""
    if not model.is_safety_critical:
        return ""
    pad = f"{indent}    "
    lines = [f"{pad}<hypnos:safetyCritical>true</hypnos:safetyCritical>"]
    absn = model.absorption or {}
    rates = absn.get("site_rates", [])
    if rates:
        lines.append(f"{pad}<hypnos:siteAbsorption>")
        lines.append(f"{pad}  <rdf:Seq>")
        for s in sorted(rates, key=lambda r: r.get("rank", 99)):
            ka = s.get("ka")
            ka_attr = f' hypnos:ka="{ka:.10g}"' if ka is not None else ""
            lines.append(f'{pad}    <rdf:li hypnos:site="{s.get("site")}" '
                         f'hypnos:rank="{s.get("rank")}"{ka_attr}/>')
        lines += [f"{pad}  </rdf:Seq>", f"{pad}</hypnos:siteAbsorption>"]
    for th in model.toxicity_thresholds:
        lines.append(
            f'{pad}<hypnos:toxicityThresholdRange hypnos:endpoint="{th.endpoint}" '
            f'hypnos:basis="{th.basis}" hypnos:low="{th.low:.10g}" hypnos:high="{th.high:.10g}" '
            f'hypnos:units="{th.units}" hypnos:tier="{th.tier}"/>')
    cc = (drug or {}).get("cardiotoxicity_class") if drug else None
    if cc:
        lines.append(
            f'{pad}<hypnos:cardiotoxicityClass hypnos:rank="{cc.get("rank")}" '
            f'hypnos:stereochemistry="{cc.get("stereochemistry")}" '
            f'hypnos:cnsToCvsMargin="{cc.get("cns_to_cvs_margin")}"/>')
    ffm = ((drug or {}).get("protein_binding") or {}).get("free_fraction_model") if drug else None
    if ffm:
        lines.append(
            f'{pad}<hypnos:freeFractionModel hypnos:type="{ffm.get("type")}" '
            f'hypnos:bindingCapacityUgMl="{ffm.get("binding_capacity_ug_ml")}" '
            f'hypnos:tier="{ffm.get("tier")}"/>')
    return "\n".join(lines)


def developmental_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the developmental extrapolation block (v0.8 §7).

    The allometric size scaling and the PMA maturation function ride as metadata so they
    travel with an SBML/PharmML model whose core is deterministic, carrying the Tier-D flag,
    the opt-in flag, and (for ``allometry_only``) the over-dose caveat. The maturation driver
    is emitted as PMA with a comment forbidding chronological-age substitution (Trap 1).
    Returns ``""`` when the model carries no developmental block."""
    dev = model.developmental_model
    if dev is None:
        return ""
    pad = f"{indent}    "
    lines = [
        f'{pad}<hypnos:developmentalModel hypnos:extrapolationBasis="{dev.extrapolation_basis}" '
        f'hypnos:evidenceTier="{dev.evidence_tier}" hypnos:appliedByDefault="{str(dev.applied_by_default).lower()}">',
    ]
    if dev.size is not None:
        s = dev.size
        desc = f' hypnos:sizeDescriptor="{s.size_descriptor}"' if s.size_descriptor else ""
        lines.append(
            f'{pad}  <hypnos:allometricSize hypnos:clExponent="{s.cl_exponent:g}" '
            f'hypnos:vExponent="{s.v_exponent:g}" hypnos:referenceWeightKg="{s.reference_weight_kg:g}" '
            f'hypnos:exponentBasis="{s.exponent_basis}"{desc} hypnos:tier="{s.tier}"/>')
    if dev.maturation is not None:
        m = dev.maturation
        lines.append(
            f'{pad}  <!-- maturation driver is PMA (post-menstrual age); a downstream consumer '
            f'MUST NOT substitute chronological age (v0.8 Trap 1) -->')
        lines.append(
            f'{pad}  <hypnos:maturationFunction hypnos:function="{m.function}" '
            f'hypnos:tm50Weeks="{m.tm50_weeks:g}" hypnos:hill="{m.hill:g}" hypnos:driver="{m.driver}" '
            f'hypnos:affectedParameter="{m.affected_parameter}" hypnos:tier="{m.tier}"/>')
    else:
        lines.append(f'{pad}  <hypnos:maturationCaveat>allometry_only — maturation un-modeled; '
                     f'neonatal clearance OVER-stated, exposure UNDER-predicted</hypnos:maturationCaveat>')
    lines.append(f'{pad}  <hypnos:developmentalExtrapolation>true</hypnos:developmentalExtrapolation>')
    lines.append(f'{pad}</hypnos:developmentalModel>')
    return "\n".join(lines)


def pharmacogenomics_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the pharmacogenomic blocks (v0.9 §7).

    Kinetic modifiers ride as ``hypnos:pharmacogenomicModifier`` metadata (genotype is a
    covariate the core formats describe only weakly); safety flags ride as
    ``hypnos:safetyCritical`` + ``hypnos:pharmacogenomicSafetyFlag``, **explicitly not** as
    any parameter change — the RDF carrier keeps the avoidance fact lossless and portable
    without ever encoding it as a number a deterministic consumer might apply. Returns ``""``
    when the model carries no pharmacogenomic blocks."""
    if not model.has_pharmacogenomics:
        return ""
    pad = f"{indent}    "
    lines: List[str] = []
    for mod in model.pharmacogenomic_modifiers:
        sf = mod.scale_factor
        sf_attr = f' hypnos:scaleFactor="{sf:.10g}"' if sf is not None else ""
        lines.append(
            f'{pad}<hypnos:pharmacogenomicModifier hypnos:gene="{mod.gene}" '
            f'hypnos:phenotype="{mod.phenotype_value}" hypnos:affectedParameter="{mod.affected_parameter}"'
            f'{sf_attr} hypnos:substrateScope="{" ".join(mod.substrate_scope)}" '
            f'hypnos:evidenceTier="{mod.evidence_tier}" hypnos:appliedByDefault="false"/>')
    for flag in model.pharmacogenomic_safety_flags:
        lines.append(
            f'{pad}<hypnos:pharmacogenomicSafetyFlag hypnos:gene="{flag.gene}" '
            f'hypnos:kind="{flag.kind}" hypnos:triggers="{" ".join(flag.trigger_agents)}" '
            f'hypnos:affectsKinetics="false" hypnos:evidenceTier="{flag.evidence_tier}"/>')
    return "\n".join(lines)


def rdf_annotation_xml(model: Model, ds=None, tier: Optional[str] = None, indent: str = "  ",
                       extra_predicates: str = "") -> str:
    """A small RDF/MIRIAM block embeddable in SBML/PharmML <annotation> elements.

    ``extra_predicates`` (e.g. :func:`variability_rdf`) is injected inside the single
    ``rdf:RDF`` element — SBML allows only one ``<annotation>`` per element, so the
    random-effects layer rides alongside the provenance rather than in a second block.
    """
    t = tier or model.tier
    uris = citation_uris(model, ds)
    bag = "\n".join(
        f'{indent}      <rdf:li rdf:resource="{u}"/>' for u in uris
    ) or f'{indent}      <rdf:li rdf:resource="urn:hypnos:no-citation"/>'
    extra = (extra_predicates + "\n") if extra_predicates else ""
    return (
        f'{indent}<annotation>\n'
        f'{indent}  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n'
        f'{indent}           xmlns:bqmodel="http://biomodels.net/model-qualifiers/"\n'
        f'{indent}           xmlns:hypnos="{HYPNOS_NS}">\n'
        f'{indent}    <hypnos:clinicalUse>{CLINICAL_USE}</hypnos:clinicalUse>\n'
        f'{indent}    <hypnos:confidenceTier>{t}</hypnos:confidenceTier>\n'
        f'{indent}    <hypnos:datasetVersion>{__version__}</hypnos:datasetVersion>\n'
        f'{indent}    <bqmodel:isDerivedFrom>\n'
        f'{indent}      <rdf:Bag>\n{bag}\n{indent}      </rdf:Bag>\n'
        f'{indent}    </bqmodel:isDerivedFrom>\n'
        f'{extra}'
        f'{indent}  </rdf:RDF>\n'
        f'{indent}</annotation>'
    )
