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
    return {
        "hypnos:datasetVersion": __version__,
        "hypnos:modelId": model.id,
        "hypnos:confidenceTier": tier or model.tier,
        "hypnos:reviewStatus": model.review_status,
        "hypnos:clinicalUse": CLINICAL_USE,
        "bqmodel:isDerivedFrom": citation_uris(model, ds),
    }


def banner(model: Model, tier: Optional[str] = None) -> str:
    """Human-readable banner for text-format exports (NONMEM, R, etc.)."""
    t = tier or model.tier
    return (
        "============================================================\n"
        "  HYPNOS EXPORT — NOT FOR CLINICAL USE\n"
        f"  clinicalUse = {CLINICAL_USE}\n"
        f"  model: {model.id}  (tier {t}, review: {model.review_status})\n"
        f"  dataset version: {__version__}\n"
        "  Research / education / simulation only. Not a dosing tool.\n"
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
