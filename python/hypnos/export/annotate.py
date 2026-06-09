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


def rdf_annotation_xml(model: Model, ds=None, tier: Optional[str] = None, indent: str = "  ") -> str:
    """A small RDF/MIRIAM block embeddable in SBML/PharmML <annotation> elements."""
    t = tier or model.tier
    uris = citation_uris(model, ds)
    bag = "\n".join(
        f'{indent}      <rdf:li rdf:resource="{u}"/>' for u in uris
    ) or f'{indent}      <rdf:li rdf:resource="urn:hypnos:no-citation"/>'
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
        f'{indent}  </rdf:RDF>\n'
        f'{indent}</annotation>'
    )
