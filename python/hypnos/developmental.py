"""Developmental-pharmacology kernels — the v0.8 size + maturation layer.

Two pure, forward-only primitives compose an adult disposition model into a
neonate/infant by the two mechanistic facts developmental PK turns on (v0.8 §3):

* :func:`allometric_scale` — theory-based size scaling (clearance ∝ weight^¾,
  volume ∝ weight¹), relative to a reference individual;
* :func:`maturation_factor` — the Anderson–Holford PMA-driven clearance-maturation
  sigmoid.

Both are applied to a model's *volumes/clearances* view (before the structural→micro
conversion in :mod:`hypnos.reference`), so the matrix-exponential solver and every
exporter keep their single canonical representation — the same composition point
v0.2's ``sample_individual`` and v0.7's covariate sampling use.

There is **no inverse control** (spec §10): nothing here computes a pediatric dose,
a per-kg scaling to a target, or a maturation-corrected dose. The honest output is
the extrapolation *spread* and its Tier-D labels, never a number to administer.
"""
from __future__ import annotations

from typing import Dict, Optional

from .models import DevelopmentalModel

# Volumes/clearances symbols and which allometric exponent each takes: clearances
# (including inter-compartmental) scale on the ¾-type exponent, volumes on the 1.0-type.
_CLEARANCE_SYMBOLS = ("Cl1", "Cl2", "Cl3")
_VOLUME_SYMBOLS = ("V1", "V2", "V3")


def allometric_scale(
    theta_ref: Dict[str, float],
    weight: float,
    *,
    cl_exponent: float = 0.75,
    v_exponent: float = 1.0,
    reference_weight_kg: float = 70.0,
) -> Dict[str, float]:
    """Scale reference (adult) volumes/clearances to a patient's size by allometry.

    ``theta_ref`` is a volumes/clearances dict (``V1``/``V2``/``V3``/``Cl1``/``Cl2``/
    ``Cl3``/``ke0``). Clearances are multiplied by ``(weight/ref)^cl_exponent``, volumes
    by ``(weight/ref)^v_exponent``; ``ke0`` (a first-order rate, 1/min) is left untouched
    (it is not a size-scaled disposition parameter). Returns a new dict — never mutates
    the input.

    The theory: physiological rates are non-linear in weight, so a 3 kg neonate does not
    clear at 3/70 of the adult rate. A naive linear-per-kg rescaling (see
    :func:`linear_per_kg_scale`) gets this wrong — that divergence is the safety message
    (v0.8 §6)."""
    if reference_weight_kg <= 0:
        raise ValueError("reference_weight_kg must be positive")
    ratio = float(weight) / float(reference_weight_kg)
    out = dict(theta_ref)
    for sym in _CLEARANCE_SYMBOLS:
        if sym in out and out[sym]:
            out[sym] = out[sym] * ratio ** cl_exponent
    for sym in _VOLUME_SYMBOLS:
        if sym in out and out[sym]:
            out[sym] = out[sym] * ratio ** v_exponent
    return out


def linear_per_kg_scale(
    theta_ref: Dict[str, float],
    weight: float,
    *,
    reference_weight_kg: float = 70.0,
) -> Dict[str, float]:
    """The naive, *visibly wrong* linear-per-kg adult rescaling (v0.8 §6).

    Scales every disposition parameter linearly with weight — the dangerous-but-intuitive
    "just scale the adult mg/kg down" mental model. It is rendered ONLY as a labeled
    reference line diverging from the mechanistic extrapolation, never offered as a number
    to administer: ignoring allometry mis-states volume and ignoring maturation over-doses
    the neonate. Exposed so the divergence can be *shown*, which is the teaching/safety point."""
    ratio = float(weight) / float(reference_weight_kg)
    out = dict(theta_ref)
    for sym in (*_CLEARANCE_SYMBOLS, *_VOLUME_SYMBOLS):
        if sym in out and out[sym]:
            out[sym] = out[sym] * ratio
    return out


def maturation_factor(pma_weeks: float, tm50_weeks: float, hill: float) -> float:
    """The Anderson–Holford clearance-maturation sigmoid (v0.8 §3).

        MF(PMA) = PMA^Hill / (TM50^Hill + PMA^Hill)

    A verbatim transcription of the curated ``function`` field. Driven by **PMA**
    (post-menstrual age, weeks) — never chronological age (the cardinal pediatric trap,
    v0.8 §4 Trap 1). Returns a value in (0, 1]: a term neonate (PMA ≈ 40 wk) has only a
    fraction of the size-adjusted adult clearance. Never silently defaults to 1.0 when
    maturation is *uncurated* — that is a gap (handled by the caller), not a mature patient."""
    p = float(pma_weeks)
    if p <= 0:
        return 0.0
    num = p ** hill
    return num / (float(tm50_weeks) ** hill + num)


def apply_developmental(
    theta_ref: Dict[str, float],
    dev: DevelopmentalModel,
    patient: Dict[str, float],
) -> Dict[str, float]:
    """Compose an adult model's volumes/clearances into a child by size + maturation.

    Applies the curated :class:`~hypnos.models.SizeModel` allometry, then (where curated)
    multiplies the maturation-``affected_parameter`` clearances by ``MF(PMA)``. Returns a
    new volumes/clearances dict. ``allometry_only`` (no maturation block) is honest about
    its over-dose direction at the call site (v0.8 §5) — this function applies exactly what
    is curated and nothing more (never-invent)."""
    out = dict(theta_ref)
    weight = float(patient.get("weight", 70.0))
    if dev.size is not None:
        out = allometric_scale(
            out, weight,
            cl_exponent=dev.size.cl_exponent, v_exponent=dev.size.v_exponent,
            reference_weight_kg=dev.size.reference_weight_kg,
        )
    if dev.maturation is not None:
        pma = patient.get("pma_weeks")
        if pma is None:
            raise ValueError(
                "developmental maturation needs patient['pma_weeks'] (PMA, not chronological "
                "age) — the maturation clock is PMA by construction (v0.8 §4 Trap 1)")
        mf = maturation_factor(pma, dev.maturation.tm50_weeks, dev.maturation.hill)
        affected = dev.maturation.affected_parameter
        # 'Cl' is shorthand for all clearance pathways; a specific symbol scales only it.
        targets = _CLEARANCE_SYMBOLS if affected in ("Cl", "CL", "clearance") else (affected,)
        for sym in targets:
            if sym in out and out[sym]:
                out[sym] = out[sym] * mf
    return out


def maturation_value(dev: DevelopmentalModel, patient: Dict[str, float]) -> Optional[float]:
    """The maturation factor MF(PMA) for this patient, or None when maturation is uncurated."""
    if dev.maturation is None:
        return None
    pma = patient.get("pma_weeks")
    if pma is None:
        return None
    return maturation_factor(pma, dev.maturation.tm50_weeks, dev.maturation.hill)
