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
    }
