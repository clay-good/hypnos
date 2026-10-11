"""Shared estimation-uncertainty projection for the exporters (v0.3 §8).

The estimation layer is the SE/RSE on each *typical value* θ plus (where published)
the estimate correlation matrix — the ``$COV`` output. It is **reducible** (shrinks
with more data) and a different object from the between-subject Ω that
:mod:`._variability` projects; every carrier built on these helpers says so in its
own idiom, because conflating an RSE with a CV is the cardinal trap (v0.3 §4 Trap 1).

Two facts every format must keep visible:

* the SEs are curated on the **curated typical value** (``value.central`` — the
  source's reference individual), *not* on the parameters an export instantiates at
  its own reference patient. Each carrier therefore emits the curated point estimate
  beside its SE, and the scale-free RSE that survives a covariate rescaling;
* a ``marginal`` model publishes no estimate correlations, so its covariance is the
  diagonal ``se²`` under a *recorded* independence assumption (v0.3 §5) — never a
  correlation invented to fill the matrix.

Nothing here invents a number: a parameter without a curated SE is absent, and a model
with ``uncertainty_status: none`` yields nothing (the never-synthesize rule).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from ._variability import VC_ORDER


@dataclass(frozen=True)
class ThetaEstimate:
    """One curated θ with its estimation uncertainty, normalized for export."""

    symbol: str
    central: float                 # the curated typical value the SE is ON
    se: float                      # as curated, on ``scale``
    scale: str                     # "natural" | "log"
    se_natural: float              # natural-scale SE at ``central`` (log: |θ|·se)
    rse_percent: Optional[float]   # recomputed from ``se`` — never the transcribed copy
    method: str
    tier: str
    ci95: Optional[Tuple[float, float]] = None


def theta_estimates(model) -> List[ThetaEstimate]:
    """Curated per-θ estimation uncertainty in the canonical η/θ order (empty if none).

    ``rse_percent`` is recomputed from ``se`` on its declared scale so an export can
    never disagree with the canonical SE (the estimation analog of v0.2 Trap 1).
    """
    if not model.has_published_estimation:
        return []
    found = {}
    for p in model.parameters:
        e = p.estimation
        if e is None or e.se is None or p.central is None:
            continue
        central = float(p.central)
        scale = e.scale or "natural"
        se = float(e.se)
        ci = e.ci95 or {}
        found[p.symbol] = ThetaEstimate(
            symbol=p.symbol, central=central, se=se, scale=scale,
            se_natural=abs(central) * se if scale == "log" else se,
            rse_percent=e.rse_from_se(central), method=e.method, tier=e.tier,
            ci95=((float(ci["low"]), float(ci["high"]))
                  if ci.get("low") is not None and ci.get("high") is not None else None),
        )
    order = [s for s in VC_ORDER if s in found] + [s for s in found if s not in VC_ORDER]
    return [found[s] for s in order]


def estimate_correlations(model) -> List[Tuple[str, str, float]]:
    """Published estimate correlations as ``[(sym_a, sym_b, r)]`` (empty if none)."""
    ec = model.estimate_covariance
    if ec is None:
        return []
    out: List[Tuple[str, str, float]] = []
    for c in ec.correlations:
        pair = c.get("between")
        r = c.get("correlation")
        if pair and len(pair) == 2 and r is not None:
            out.append((pair[0], pair[1], float(r)))
    return out


def theta_covariance(model) -> Optional[Tuple[List[str], List[List[float]]]]:
    """The estimate covariance of θ on the natural scale: ``(symbols, cov)`` or ``None``.

    Diagonal ``se²``; off-diagonals ``r·seᵢ·seⱼ`` only for *published* pairs whose
    members both carry a curated SE. Unpublished pairs stay 0 — the independence
    assumption a ``marginal`` model records, not a claim that the estimates are
    uncorrelated.
    """
    ests = theta_estimates(model)
    if not ests:
        return None
    syms = [e.symbol for e in ests]
    se = {e.symbol: e.se_natural for e in ests}
    pos = {s: i for i, s in enumerate(syms)}
    cov = [[0.0] * len(syms) for _ in syms]
    for s, i in pos.items():
        cov[i][i] = se[s] ** 2
    for a, b, r in estimate_correlations(model):
        if a in pos and b in pos:
            cov[pos[a]][pos[b]] = cov[pos[b]][pos[a]] = r * se[a] * se[b]
    return syms, cov


def covariance_step_note(model) -> Optional[str]:
    """State the covariance-step outcome explicitly when the source reports one (Trap 4)."""
    ec = model.estimate_covariance
    if ec is None or ec.covariance_step_succeeded is None:
        return None
    return ("covariance step succeeded" if ec.covariance_step_succeeded
            else "covariance step FAILED — absent SEs are MISSING, not zero")


def independence_caveat(model) -> Optional[str]:
    """The recorded caveat for a covariance carried without published correlations."""
    if not theta_estimates(model) or estimate_correlations(model):
        return None
    return ("estimate correlations not published — off-diagonals are 0 under a RECORDED "
            "independence assumption (v0.3 §5), not a claim the estimates are uncorrelated")
