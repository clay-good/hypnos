"""The covariate-equation library — the v0.7 covariate sublayer.

Pure, validity-bounded body-size/composition equations (lean body mass, fat-free
mass) that a covariate-scaled model is *derived with*. :func:`evaluate` dispatches
to the named equation, checks the equation's **own** validity envelope, and attaches
the equation's known-failure-mode warning when violated — the James LBM inversion is
the canonical entry, surfaced (never silently "fixed") at its source.

This generalizes :func:`hypnos.reference.lbm_james` into a small registered library,
each member a verbatim transcription of a curated ``dataset/covariate_equations/``
record's ``form``. The choice of equation (James vs. Janmahasatian) silently
reparameterizes a model, so it is a first-class, validated, cited object here.

There is no inverse control (spec §10): an equation maps a patient → a derived
covariate. Nothing searches over covariates to hit a target.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from .load import Dataset, load
from .models import worst_tier
from .reference import bmi as _bmi, ffm_al_sallami, ffm_janmahasatian, lbm_james


# --------------------------------------------------------------------------- #
# Covariate-value access (point value, or the mean of a caller-supplied dist)
# --------------------------------------------------------------------------- #
def _num(patient: Dict[str, Any], key: str) -> float:
    """Read a covariate that may be a scalar or a ``{mean, sd, ...}`` distribution.

    A caller-supplied distribution (the v0.7 covariate band, future C2) collapses to
    its mean for a point evaluation; a scalar is used directly. Never invents a value.
    """
    v = patient.get(key)
    if isinstance(v, dict):
        v = v.get("mean", v.get("central"))
    if v is None:
        raise ValueError(f"covariate equation requires patient covariate '{key}'")
    return float(v)


def _sex(patient: Dict[str, Any]) -> str:
    return str(patient.get("sex", "M"))


# Pure equation implementations, keyed by id — verbatim transcriptions of each
# library record's `form`. Each maps a patient dict to the derived value (kg).
def _james_1976(p: Dict[str, Any]) -> float:
    return lbm_james(_num(p, "weight"), _num(p, "height"), _sex(p))


def _janmahasatian_2005(p: Dict[str, Any]) -> float:
    return ffm_janmahasatian(_num(p, "weight"), _num(p, "height"), _sex(p))


def _al_sallami_2015(p: Dict[str, Any]) -> float:
    return ffm_al_sallami(_num(p, "weight"), _num(p, "height"), _num(p, "age"), _sex(p))


EQUATIONS: Dict[str, Callable[[Dict[str, Any]], float]] = {
    "james_1976": _james_1976,
    "janmahasatian_2005": _janmahasatian_2005,
    "al_sallami_2015": _al_sallami_2015,
}


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #
@dataclass
class EquationResult:
    """The value of one covariate equation for one patient, plus its honesty flags."""

    equation_id: str
    quantity: str
    value: float
    tier: str
    out_of_envelope: bool = False
    inverted: bool = False              # the James-style peak-then-decline non-physical inversion
    warnings: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Envelope evaluation
# --------------------------------------------------------------------------- #
def _patient_axis(patient: Dict[str, Any], axis: str) -> Optional[float]:
    if axis == "bmi_kg_m2":
        try:
            return _bmi(_num(patient, "weight"), _num(patient, "height"))
        except ValueError:
            return None
    key = {"age_years": "age", "weight_kg": "weight", "height_cm": "height"}.get(axis)
    if key is None:
        return None
    try:
        return _num(patient, key)
    except ValueError:
        return None


def _check_envelope(record: Dict[str, Any], patient: Dict[str, Any]) -> List[str]:
    """Return human-readable violation messages for the equation's OWN envelope."""
    out: List[str] = []
    env = record.get("validity_envelope") or {}
    for axis, rng in env.items():
        if not isinstance(rng, dict):
            continue
        val = _patient_axis(patient, axis)
        if val is None:
            continue
        lo, hi = rng.get("min"), rng.get("max")
        if lo is not None and val < lo:
            out.append(f"{axis}={val:g} below the equation's validity range (>= {lo:g})")
        if hi is not None and val > hi:
            out.append(f"{axis}={val:g} above the equation's validity range (<= {hi:g})")
    return out


def _is_inverted(equation_id: str, patient: Dict[str, Any]) -> bool:
    """Tested-property check: is the derived value *decreasing* in weight at this patient?

    The James LBM inversion is exactly this — above its envelope, computed LBM peaks and
    then declines with rising weight (non-physical). A one-sided finite difference in
    weight (holding height/age) detects it generically, for any equation in the library.
    """
    fn = EQUATIONS[equation_id]
    try:
        w = _num(patient, "weight")
    except ValueError:
        return False
    p = dict(patient)
    base = fn({**p, "weight": w})
    bumped = fn({**p, "weight": w * 1.01})
    return bumped < base


# --------------------------------------------------------------------------- #
# evaluate — the public entry point
# --------------------------------------------------------------------------- #
def evaluate(
    equation_id: str,
    patient: Dict[str, Any],
    *,
    ds: Optional[Dataset] = None,
) -> EquationResult:
    """Evaluate a named covariate equation for a patient, honestly.

    Dispatches to the pure equation, checks the equation's own ``validity_envelope``,
    and — when violated — tiers the result to D and attaches the equation's
    ``known_failure_modes`` warning (the James inversion is the canonical one). The
    inversion is *surfaced*, never silently substituted away (spec §10): a model that
    is derived with James in obesity rests on an inverted LBM, and that is the honest
    result a downstream consumer must see.
    """
    if ds is None:
        ds = load()
    if equation_id not in EQUATIONS:
        raise KeyError(
            f"no implemented covariate equation '{equation_id}' "
            f"(known: {sorted(EQUATIONS)})"
        )
    record = ds.covariate_equations.get(equation_id, {})
    value = EQUATIONS[equation_id](patient)
    tier = record.get("tier", "D")

    warnings: List[str] = []
    viols = _check_envelope(record, patient)
    out_of_envelope = bool(viols)
    inverted = _is_inverted(equation_id, patient)

    if out_of_envelope:
        tier = worst_tier([tier, "D"])
        detail = "; ".join(viols)
        for fm in record.get("known_failure_modes", []):
            if fm.get("action") in ("tier_down_to_D", "exclude"):
                warnings.append(
                    f"{equation_id} ({record.get('quantity', '?')}) outside its validity "
                    f"envelope ({detail}): {fm.get('behavior')} -> covariate layer tier D"
                )
                break
        else:
            warnings.append(
                f"{equation_id} ({record.get('quantity', '?')}) outside its validity "
                f"envelope ({detail}) -> covariate layer tier D"
            )
    if inverted:
        warnings.append(
            f"{equation_id} has INVERTED at this patient: the derived {record.get('quantity', 'value')} "
            f"is DECREASING with weight (non-physical) — value {value:.1f} kg is not trustworthy"
        )

    return EquationResult(
        equation_id=equation_id,
        quantity=record.get("quantity", "?"),
        value=value,
        tier=tier,
        out_of_envelope=out_of_envelope,
        inverted=inverted,
        warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# Model-level helpers
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Covariate-VALUE uncertainty — caller-supplied input distributions (v0.7 C2)
# --------------------------------------------------------------------------- #
# A covariate may be supplied as a scalar (exact) or as a distribution dict
# {"mean": .., "sd": ..} (or {"central", "cv"}); a distribution is what the v0.7
# covariate band propagates. Hypnos NEVER invents a distribution — absent one, the
# covariate is treated as exact (spec §2.3/§5).
_NONNEGATIVE = ("weight", "height", "age", "crcl_ml_min", "albumin_g_dl")


def _is_distribution(v: Any) -> bool:
    return isinstance(v, dict) and ("sd" in v or "cv" in v)


def has_covariate_distribution(patient: Dict[str, Any]) -> bool:
    """True when the patient supplies at least one covariate-value distribution."""
    return any(_is_distribution(v) for v in patient.values())


def distribution_keys(patient: Dict[str, Any]) -> List[str]:
    """The covariate names carrying a value distribution (drives the §7.3 readout)."""
    return [k for k, v in patient.items() if _is_distribution(v)]


def point_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Collapse any covariate distribution to its point (mean) value.

    Scalars pass through unchanged, so a scalar-covariate patient is byte-identical —
    every deterministic path (kernel, envelope, PD link) runs on this point vector."""
    out: Dict[str, Any] = {}
    for k, v in patient.items():
        out[k] = v.get("mean", v.get("central")) if isinstance(v, dict) else v
    return out


def sample_covariate_vector(patient: Dict[str, Any], rng) -> Dict[str, Any]:
    """Draw one perturbed covariate vector from caller-supplied ``{mean, sd, dist}`` marginals.

    Covariates given as scalars are held fixed (never invented). Marginals are perturbed
    **independently** — Hypnos assumes no covariate covariance absent a caller-supplied one,
    the same honest default as v0.2's ``omega_block.complete = false`` (spec §14). Draws come
    from the seeded ``rng`` so a covariate band is byte-reproducible (spec §6).
    """
    out: Dict[str, Any] = {}
    for k, v in patient.items():
        if _is_distribution(v):
            mean = v.get("mean", v.get("central"))
            sd = v.get("sd")
            if sd is None and v.get("cv") is not None:
                sd = abs(float(mean)) * float(v["cv"])
            dist = v.get("dist", "normal")
            if dist == "lognormal" and mean:
                draw = float(mean) * math.exp(rng.normal(0.0, float(sd) / float(mean)))
            else:
                draw = rng.normal(float(mean), float(sd))
            out[k] = max(float(draw), 1e-6) if k in _NONNEGATIVE else float(draw)
        elif isinstance(v, dict):
            out[k] = v.get("mean", v.get("central"))
        else:
            out[k] = v
    return out


def covariate_layer_tier(model, patient: Dict[str, Any], ds: Dataset) -> Optional[str]:
    """Worst tier the model's bound covariate equations contribute at this patient.

    Returns None when the model declares no covariate_model (an explicit gap, never an
    assumption). An equation evaluated outside its own envelope forces D (spec §5), so
    this composes into the worst-tier chain like any other component.
    """
    cm = model.covariate_model
    if cm is None:
        return None
    tiers: List[str] = []
    for d in cm.derived_inputs:
        if d.equation in EQUATIONS:
            tiers.append(evaluate(d.equation, patient, ds=ds).tier)
        else:
            tiers.append(d.tier)
    return worst_tier(tiers) if tiers else None
