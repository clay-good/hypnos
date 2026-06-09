"""Shared random-effects (Ω/Σ) projection for the population exporters (v0.2 §8).

Each pharmacometric target expresses the NLME random-effects layer in its own
idiom, but they all draw on the same three curated pieces:

* the **Ω diagonal** — per-structural-parameter between-subject variance ``omega2``
  (η-scale), ordered canonically so every export wires the η's the same way;
* the **Σ residual** — the observation-error model, normalized so a format can
  render it without re-deriving variance-vs-SD (the most common transcription
  trap, spec §4 / §3.2);
* the **off-diagonal Ω correlations** — published correlated η's, projected to a
  covariance matrix only over a *contiguous* η span (so NONMEM's ``$OMEGA BLOCK``
  stays valid), and otherwise left as an honest diagonal-plus-caveat.

Nothing here invents a band: a parameter absent from ``bsv_omegas`` carries no
η, and a model with no ``residual_error`` gets ``None`` — the never-synthesize
rule (spec §5) lives in the dataset, and these helpers only ever read it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple

# Canonical structural order shared by every population export: clearances and
# volumes interleaved as a clinician reads them, then the effect-site link. The
# η wired onto each parameter follows this order so NONMEM ETA(k), the PharmML
# RandomEffect list, the rxode2 ``lotri`` rows, and the Pumas ``@random`` block
# all agree.
VC_ORDER = ["Cl1", "V1", "Cl2", "V2", "Cl3", "V3", "ke0"]


def omega_diagonal(model) -> List[Tuple[str, float, float]]:
    """Return ``[(symbol, omega2, cv_percent)]`` for structural params carrying BSV.

    Ordered by :data:`VC_ORDER` (any unrecognized symbol trails, stably). ``cv_percent``
    is recomputed from ``omega2`` via the exact log-normal relation so it never
    silently disagrees with the canonical variance (spec §4 Trap 1).
    """
    omegas = model.bsv_omegas()
    ordered = [s for s in VC_ORDER if s in omegas] + [s for s in omegas if s not in VC_ORDER]
    return [(s, omegas[s], 100.0 * math.sqrt(math.exp(omegas[s]) - 1.0)) for s in ordered]


@dataclass(frozen=True)
class ResidualSpec:
    """Σ normalized to canonical scales (variance for proportional, SD for additive/log)."""

    model: str  # log | proportional | additive | combined
    label: str
    log_sd: Optional[float] = None
    prop_var: Optional[float] = None
    add_sd: Optional[float] = None


def _as_sd(d) -> Optional[float]:
    if not d:
        return None
    if d.get("sd") is not None:
        return float(d["sd"])
    if d.get("variance") is not None:
        return math.sqrt(float(d["variance"]))
    return None


def _as_var(d) -> Optional[float]:
    if not d:
        return None
    if d.get("variance") is not None:
        return float(d["variance"])
    if d.get("cv_percent") is not None:
        return (float(d["cv_percent"]) / 100.0) ** 2
    if d.get("sd") is not None:
        return float(d["sd"]) ** 2
    return None


def residual_spec(model) -> Optional[ResidualSpec]:
    """Normalize a curated ``residual_error`` into a :class:`ResidualSpec` (or None)."""
    re_ = model.residual_error
    if re_ is None:
        return None
    m = re_.model
    if m == "log":
        sd = _as_sd(re_.log)
        if sd is None:
            return None
        return ResidualSpec("log", "log-additive (≈ proportional on natural scale)", log_sd=sd)
    if m == "proportional":
        var = _as_var(re_.proportional)
        if var is None:
            return None
        return ResidualSpec("proportional", "proportional", prop_var=var)
    if m == "additive":
        sd = _as_sd(re_.additive)
        if sd is None:
            return None
        return ResidualSpec("additive", "additive", add_sd=sd)
    if m == "combined":
        pv = _as_var(re_.proportional)
        asd = _as_sd(re_.additive)
        if pv is None or asd is None:
            return None
        return ResidualSpec("combined", "combined (proportional + additive)",
                            prop_var=pv, add_sd=asd)
    return None


def omega_correlations(model) -> List[Tuple[str, str, float]]:
    """Published off-diagonal η correlations as ``[(sym_a, sym_b, r)]`` (empty if none)."""
    ob = model.omega_block
    if ob is None:
        return []
    out: List[Tuple[str, str, float]] = []
    for c in ob.correlations:
        pair = c.get("between")
        r = c.get("correlation")
        if pair and len(pair) == 2 and r is not None:
            out.append((pair[0], pair[1], float(r)))
    return out


def contiguous_block(model) -> Optional[Tuple[List[str], List[List[float]]]]:
    """Project a curated ``omega_block`` to a dense lower-triangular covariance matrix.

    Returns ``(symbols, cov)`` where ``symbols`` is a **contiguous** prefix of the BSV
    diagonal (in :data:`VC_ORDER`) that exactly spans every correlated parameter, and
    ``cov`` is the symmetric covariance matrix over them — the form NONMEM ``$OMEGA
    BLOCK`` and PharmML's correlated ``VariabilityLevel`` need. Off-diagonals are
    ``r·√(ωᵢ²·ωⱼ²)``; pairs not published within the block are 0 covariance, which
    a ``complete: true`` block asserts (spec §3.3).

    Returns ``None`` — i.e. *fall back to an honest diagonal + caveat* — when the
    block is incomplete, or when the correlated parameters are not a contiguous η
    span (a BLOCK over non-adjacent η's would be invalid), or when any block member
    lacks a diagonal ω². No covariance is ever invented to force a block.
    """
    corrs = omega_correlations(model)
    if not corrs:
        return None
    ob = model.omega_block
    if ob is None or not ob.complete:
        return None

    omegas = model.bsv_omegas()
    members = set()
    for a, b, _ in corrs:
        members.update((a, b))
    if any(m not in omegas for m in members):
        return None
    if any(m not in VC_ORDER for m in members):
        return None

    # Work in the full η-index space (VC_ORDER), not the BSV-only diagonal, because
    # a $OMEGA BLOCK must span an unbroken run of η's with no fixed-zero gap.
    positions = sorted(VC_ORDER.index(m) for m in members)
    span = list(range(positions[0], positions[-1] + 1))
    # front-anchored (so the BLOCK precedes the remaining diagonal η's) and gap-free
    # (every parameter in the span carries its own published BSV — no 0-FIX inside).
    if positions[0] != 0 or any(VC_ORDER[i] not in omegas for i in span):
        return None

    syms = [VC_ORDER[i] for i in span]
    n = len(syms)
    cov = [[0.0] * n for _ in range(n)]
    for i, s in enumerate(syms):
        cov[i][i] = omegas[s]
    pos = {s: i for i, s in enumerate(syms)}
    for a, b, r in corrs:
        i, j = pos[a], pos[b]
        c = r * math.sqrt(omegas[a] * omegas[b])
        cov[i][j] = cov[j][i] = c
    return syms, cov
