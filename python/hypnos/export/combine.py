"""COMBINE ``.omex`` archive exporter.

Bundles a model's representations — SBML (master) + PharmML + TCI-JSON — together
with a provenance ``metadata.rdf`` and a ``citations.bib`` into a single
COMBINE-compliant ``.omex`` archive (a ZIP with an ``omex-manifest`` describing
each entry's format). This is the durable, self-describing distribution unit:
one file carrying the model, its interop projections, its confidence tier, its
DOI/PMID provenance, and the universal ``clinicalUse = PROHIBITED`` flag.

Archives are written **deterministically** (fixed entry timestamps) so the same
dataset version always produces byte-identical bytes — a reproducibility
guarantee, not a nicety.
"""
from __future__ import annotations

import zipfile
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

from .. import CLINICAL_USE, __version__
from . import annotate, bibtex, pharmml, sbml, tci_json
from ._common import resolve_patient, safe_name

# COMBINE / MIRIAM format URIs
_FMT_MANIFEST = "http://identifiers.org/combine.specifications/omex-manifest"
_FMT_OMEX = "http://identifiers.org/combine.specifications/omex"
_FMT_SBML = "http://identifiers.org/combine.specifications/sbml"
_FMT_META = "http://identifiers.org/combine.specifications/omex-metadata"
_FMT_PHARMML = "http://purl.org/NET/mediatypes/application/pharmml+xml"
_FMT_JSON = "http://purl.org/NET/mediatypes/application/json"
_FMT_BIBTEX = "http://purl.org/NET/mediatypes/application/x-bibtex"
_FMT_TEXT = "http://purl.org/NET/mediatypes/text/plain"

# Fixed timestamp for byte-reproducible archives (ZIP epoch lower bound).
_FIXED_DATE = (1980, 1, 1, 0, 0, 0)


def _manifest(entries: List[Tuple[str, str, bool]]) -> str:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">',
        f'  <content location="." format="{_FMT_OMEX}"/>',
        f'  <content location="./manifest.xml" format="{_FMT_MANIFEST}"/>',
    ]
    for loc, fmt, master in entries:
        m = ' master="true"' if master else ""
        lines.append(f'  <content location="./{loc}" format="{fmt}"{m}/>')
    lines.append("</omexManifest>")
    return "\n".join(lines) + "\n"


def _metadata_rdf(models, ds) -> str:
    """A small RDF/XML provenance document for the archive."""
    blocks = []
    for m in models:
        prov = annotate.provenance(m, ds)
        cites = "".join(
            f'      <bqmodel:isDerivedFrom rdf:resource="{u}"/>\n'
            for u in prov["bqmodel:isDerivedFrom"]
        )
        blocks.append(
            f'  <rdf:Description rdf:about="#{safe_name(m)}">\n'
            f'    <hypnos:modelId>{m.id}</hypnos:modelId>\n'
            f'    <hypnos:confidenceTier>{m.tier}</hypnos:confidenceTier>\n'
            f'    <hypnos:reviewStatus>{m.review_status}</hypnos:reviewStatus>\n'
            f'    <hypnos:clinicalUse>{CLINICAL_USE}</hypnos:clinicalUse>\n'
            f'    <hypnos:datasetVersion>{__version__}</hypnos:datasetVersion>\n'
            f'{cites}'
            f'  </rdf:Description>'
        )
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n'
        '         xmlns:bqmodel="http://biomodels.net/model-qualifiers/"\n'
        '         xmlns:hypnos="https://w3id.org/hypnos/terms#">\n'
        + "\n".join(blocks) + "\n</rdf:RDF>\n"
    )


def _zip(files: Dict[str, str]) -> bytes:
    """Write a deterministic ZIP (fixed timestamps) from {name: text}."""
    buf = BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in files:  # insertion order is deterministic
            info = zipfile.ZipInfo(filename=name, date_time=_FIXED_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, files[name])
    return buf.getvalue()


def build_model_archive(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> bytes:
    """Build a single-model ``.omex`` archive (bytes)."""
    pat = resolve_patient(model, patient)
    name = safe_name(model)
    files: Dict[str, str] = {}
    entries: List[Tuple[str, str, bool]] = []

    sbml_name = f"{name}.sbml.xml"
    files[sbml_name] = sbml.build(model, ds, pat)
    entries.append((sbml_name, _FMT_SBML, True))  # master

    pharmml_name = f"{name}.pharmml.xml"
    files[pharmml_name] = pharmml.build(model, ds, pat)
    entries.append((pharmml_name, _FMT_PHARMML, False))

    tci_name = f"{name}.tci.json"
    files[tci_name] = tci_json.build(model, ds, pat)
    entries.append((tci_name, _FMT_JSON, False))

    files["metadata.rdf"] = _metadata_rdf([model], ds)
    entries.append(("metadata.rdf", _FMT_META, False))

    if ds is not None:
        files["citations.bib"] = bibtex.build_for_model(model, ds)
        entries.append(("citations.bib", _FMT_BIBTEX, False))

    # manifest first in iteration order, but ZIP order is fine either way
    out = {"manifest.xml": _manifest(entries)}
    out.update(files)
    return _zip(out)


def build_dataset_archive(ds, models=None, patient: Optional[Dict[str, Any]] = None) -> bytes:
    """Build a whole-dataset ``.omex`` bundling every PK model's SBML + shared provenance."""
    if models is None:
        models = [m for m in ds.models if m.purpose == "pk" and m.kernel_implemented]
    files: Dict[str, str] = {}
    entries: List[Tuple[str, str, bool]] = []

    for i, m in enumerate(models):
        name = safe_name(m)
        sbml_name = f"{name}.sbml.xml"
        files[sbml_name] = sbml.build(m, ds, resolve_patient(m, patient))
        entries.append((sbml_name, _FMT_SBML, i == 0))

    files["citations.bib"] = bibtex.build(ds, models)
    entries.append(("citations.bib", _FMT_BIBTEX, False))
    files["metadata.rdf"] = _metadata_rdf(models, ds)
    entries.append(("metadata.rdf", _FMT_META, False))
    files["README.txt"] = (
        f"Hypnos dataset COMBINE archive — v{__version__}\n"
        f"clinicalUse = {CLINICAL_USE}\n"
        f"Contains {len(models)} PK model(s) as SBML, with provenance and citations.\n"
        "NOT FOR CLINICAL USE. Research / education / simulation only.\n"
    )
    entries.append(("README.txt", _FMT_TEXT, False))

    out = {"manifest.xml": _manifest(entries)}
    out.update(files)
    return _zip(out)


def model_filename(model) -> str:
    return f"{safe_name(model)}.omex"
