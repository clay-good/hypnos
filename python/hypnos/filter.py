"""Filtering and summary helpers over a :class:`~hypnos.load.Dataset`."""
from __future__ import annotations

from collections import Counter
from typing import Iterable, List, Optional

from .load import Dataset
from .models import Model


def select(
    ds: Dataset,
    *,
    drug: Optional[str] = None,
    subsystem: Optional[str] = None,
    purpose: Optional[str] = None,
    tier: Optional[str] = None,
    review_status: Optional[str] = None,
    kernel_only: bool = False,
) -> List[Model]:
    """Return models matching all supplied criteria."""
    out: List[Model] = []
    for m in ds:
        if drug is not None and m.drug_name != drug:
            continue
        if subsystem is not None and m.subsystem != subsystem:
            continue
        if purpose is not None and m.purpose != purpose:
            continue
        if tier is not None and m.tier != tier:
            continue
        if review_status is not None and m.review_status != review_status:
            continue
        if kernel_only and not m.kernel_implemented:
            continue
        out.append(m)
    return out


def pk_drugs(ds: Dataset) -> List[str]:
    """Drugs with at least one executable PK kernel (i.e. simulatable / comparable)."""
    out: List[str] = []
    for m in ds:
        if m.purpose == "pk" and m.kernel_implemented and m.drug_name not in out:
            out.append(m.drug_name)
    return sorted(out)


def performance_table(ds: Dataset, *, drug: Optional[str] = None) -> List[dict]:
    """Published predictive-performance metrics (Varvel's MDPE/MDAPE, wobble,
    divergence) across the dataset — one row per metric per model.

    These are the numeric counterpart to the editorial confidence tier (spec §5:
    "tier assignment can be partly numeric"). MDPE is bias (signed), MDAPE is
    inaccuracy. Each row carries its derivation/validation population and resolves
    the citation's DOI so the number is traceable to a source, never asserted bare.
    """
    rows: List[dict] = []
    for m in sorted(ds, key=lambda x: x.id):
        if drug is not None and m.drug_name != drug:
            continue
        for pp in m.predictive_performance:
            cid = pp.get("citation")
            cit = ds.citation(cid) if cid else None
            rows.append({
                "model_id": m.id,
                "tier": m.tier,
                "metric": pp["metric"],
                "value": pp["value"],
                "units": pp.get("units", ""),
                "population": pp.get("population"),
                "citation": cid,
                "doi": (cit or {}).get("doi"),
            })
    return rows


def _counter(items: Iterable[str]) -> dict:
    return dict(sorted(Counter(items).items()))


def summary(ds: Dataset) -> dict:
    """Counts by subsystem / tier / review status / purpose — what ``hypnos info`` prints."""
    models = ds.models
    return {
        "version": ds.version,
        "n_models": len(models),
        "n_drugs": len(ds.drugs),
        "n_citations": len(ds.citations),
        "by_subsystem": _counter(m.subsystem for m in models),
        "by_purpose": _counter(m.purpose for m in models),
        "by_tier": _counter(m.tier for m in models),
        "by_review_status": _counter(m.review_status for m in models),
        "kernels_implemented": sum(1 for m in models if m.kernel_implemented),
        "models_with_predictive_performance": sum(1 for m in models if m.predictive_performance),
    }
