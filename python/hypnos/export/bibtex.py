"""BibTeX citation exporter.

Flat citation export (spec §7). Emits one ``@article`` entry per citation record,
so a consumer who uses a Hypnos model can cite both Hypnos and the original
source. Used standalone and bundled inside the COMBINE ``.omex`` archive.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


def _entry(cit: Dict[str, Any]) -> str:
    fields: List[str] = []
    authors = cit.get("authors") or []
    if authors:
        fields.append(f"  author = {{{' and '.join(authors)}}}")
    if cit.get("title"):
        fields.append(f"  title = {{{cit['title']}}}")
    if cit.get("container"):
        fields.append(f"  journal = {{{cit['container']}}}")
    if cit.get("year") is not None:
        fields.append(f"  year = {{{cit['year']}}}")
    if cit.get("volume"):
        fields.append(f"  volume = {{{cit['volume']}}}")
    if cit.get("issue"):
        fields.append(f"  number = {{{cit['issue']}}}")
    if cit.get("pages"):
        fields.append(f"  pages = {{{cit['pages']}}}")
    if cit.get("doi"):
        fields.append(f"  doi = {{{cit['doi']}}}")
    if cit.get("pmid"):
        fields.append(f"  note = {{PMID: {cit['pmid']}}}")
    body = ",\n".join(fields)
    return f"@article{{{cit['id']},\n{body}\n}}\n"


def build_for_model(model, ds=None, patient=None) -> str:
    """BibTeX for a single model: its primary citation plus any per-parameter citations.

    ``patient`` is accepted (and ignored) for a uniform exporter signature.
    """
    if ds is None:
        return ""
    ids: List[str] = []
    for cid in [model.primary_citation] + [p.primary_citation for p in model.parameters]:
        if cid and cid not in ids:
            ids.append(cid)
    chunks = [f"% Citations for Hypnos model {model.id}\n"]
    for cid in ids:
        cit = ds.citation(cid)
        if cit:
            chunks.append(_entry(cit))
    return "\n".join(chunks)


def build(ds, models: Optional[list] = None) -> str:
    """BibTeX export. With ``models=None`` emits the entire citation library; with
    a model subset, emits only the citations those models reference (primary,
    per-parameter, failure-mode, predictive-performance, organ-tolerance, and the
    drug's protein-binding source)."""
    if models is None:
        ids: List[str] = list(ds.citations.keys())
    else:
        ids = []
        for m in models:
            refs = [m.primary_citation]
            refs += [p.primary_citation for p in m.parameters]
            refs += [fm.citation for fm in m.known_failure_modes]
            refs += [pp.get("citation") for pp in m.predictive_performance]
            # v0.5: organ-failure standing (model-level) + protein binding (drug-level)
            refs += [ot.get("citation") for ot in m.applicability_envelope.organ_tolerance]
            refs.append(((ds.drug(m.drug_name) or {}).get("protein_binding") or {}).get("citation"))
            for cid in refs:
                if cid and cid not in ids:
                    ids.append(cid)
    chunks = ["% Hypnos citation export — cite Hypnos AND the original source of each model used.\n"]
    for cid in ids:
        cit = ds.citation(cid)
        if cit:
            chunks.append(_entry(cit))
    return "\n".join(chunks)


def filename(model=None) -> str:
    return "citations.bib"
