"""CSV flat-parameter exporter (spec §7).

A spreadsheet-friendly projection of the dataset: one row per parameter across
all models, with the record-level tier/review status and the resolved
DOI/PMID joined in. Complements the structured exporters (NONMEM/PharmML/SBML)
for quick filtering, pivoting, and provenance lookup. Citation export is the
sibling (`hypnos.export.bibtex`).
"""
from __future__ import annotations

import csv
import io
from typing import Any, List, Optional

_HEADER = [
    "model_id", "subsystem", "drug", "purpose", "record_tier", "review_status",
    "symbol", "label", "central", "low", "high", "units", "param_tier",
    "covariate_model", "primary_citation", "doi", "pmid",
]


def _rows_for(model, ds) -> List[List[Any]]:
    rows: List[List[Any]] = []
    for p in model.parameters:
        cite_id = p.primary_citation or model.primary_citation
        cit = ds.citation(cite_id) if ds is not None else None
        rows.append([
            model.id, model.subsystem, model.drug_name, model.purpose,
            model.tier, model.review_status,
            p.symbol, p.label or "",
            "" if p.central is None else p.central,
            p.value.get("low", "") if p.value.get("low") is not None else "",
            p.value.get("high", "") if p.value.get("high") is not None else "",
            p.units or "", p.tier,
            p.covariate_model or "", cite_id,
            (cit or {}).get("doi", ""), (cit or {}).get("pmid", ""),
        ])
    return rows


def _write(rows: List[List[Any]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(_HEADER)
    w.writerows(rows)
    return buf.getvalue()


def build_for_model(model, ds=None, patient=None) -> str:
    """Flat CSV of one model's parameters (header + one row per parameter).

    ``patient`` is accepted (and ignored) for a uniform exporter signature.
    """
    return _write(_rows_for(model, ds))


def build(ds, models: Optional[list] = None) -> str:
    """Flat CSV of every parameter across the dataset (or a subset)."""
    if models is None:
        models = ds.models
    rows: List[List[Any]] = []
    for m in models:
        rows.extend(_rows_for(m, ds))
    return _write(rows)


def filename(model=None) -> str:
    return "parameters.csv"
