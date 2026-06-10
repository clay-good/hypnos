"""Inhalational-agent MAC API — the volatile (physicochemical) convention.

Volatile anaesthetics are not described by compartmental PK; they are
characterized by MAC (minimum alveolar concentration), its age correction, and
partition coefficients. This module evaluates age-corrected MAC, the MAC
fraction (a depth surrogate) given an end-tidal concentration, and the *additive*
combined MAC fraction when nitrous oxide is co-administered.

Like the rest of Hypnos, it enforces the applicability envelope (the MAC
age-correction is valid for age > 1 y) and reports the model tier. It is a
research/education tool, NOT a dosing or depth-of-anaesthesia monitor.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .load import Dataset
from .models import worst_tier
from .reference import alveolar_washin, alveolar_washout, mac_age_corrected
from .reference import mac_fraction as _mac_fraction


@dataclass
class MacResult:
    agent_id: str
    age: float
    mac40: float
    mac_age: float                 # age-corrected MAC (vol%)
    mac_awake_age: float           # age-corrected MAC-awake (vol%)
    blood_gas: float
    oil_gas: float
    tier: str
    warnings: List[str] = field(default_factory=list)
    end_tidal_pct: Optional[float] = None
    mac_fraction: Optional[float] = None         # end-tidal / MAC(age)
    combined_mac_fraction: Optional[float] = None  # incl. nitrous oxide additivity


@dataclass
class WashinResult:
    agent_id: str
    blood_gas: float               # λ, blood:gas partition coefficient
    plateau: float                 # early FA/FI "knee" = V̇_A / (V̇_A + λ·Q̇)
    tau_min: float                 # time constant τ = FRC / (V̇_A + λ·Q̇)
    t_min: float                   # the time point reported below
    fa_fi: float                   # FA/FI at t_min
    tier: str
    # the (stated, overridable) ventilation assumptions this result rests on
    alveolar_ventilation: float
    frc: float
    cardiac_output: float


@dataclass
class WashoutResult:
    agent_id: str
    blood_gas: float               # λ, blood:gas partition coefficient
    floor: float                   # early elimination floor = λ·Q̇ / (V̇_A + λ·Q̇) = 1 − plateau
    tau_min: float                 # time constant τ = FRC / (V̇_A + λ·Q̇)
    t_min: float                   # the time point reported below
    fa_fa0: float                  # FA/FA₀ at t_min
    tier: str
    # the (stated, overridable) ventilation assumptions this result rests on
    alveolar_ventilation: float
    frc: float
    cardiac_output: float


def _params(model) -> Dict[str, float]:
    return {p.symbol: p.central for p in model.parameters}


def washin(
    ds: Dataset,
    agent_id: str,
    *,
    t_min: float = 3.0,
    alveolar_ventilation: float = 4.0,
    frc: float = 2.5,
    cardiac_output: float = 5.0,
) -> WashinResult:
    """Single-compartment alveolar wash-in (FA/FI) for one volatile agent.

    A comparative, solubility-driven characterization of inhalational *uptake*
    (spec §6) from the curated blood:gas partition coefficient. Reports the early
    FA/FI plateau (the wash-in "knee"), the time constant τ, and FA/FI at ``t_min``
    minutes. Lower blood:gas solubility → higher plateau → faster wash-in. Standard
    70-kg-adult ventilation constants, all overridable; NOT a per-patient predictor.
    """
    model = ds[agent_id]
    if model.purpose != "physicochemical":
        raise ValueError(f"{agent_id} has purpose '{model.purpose}', expected 'physicochemical'")
    lam = _params(model).get("blood_gas")
    if lam is None:
        raise ValueError(f"{agent_id} has no blood_gas partition coefficient")
    fa_fi, plateau, tau = alveolar_washin(
        lam, np.array([t_min]), alveolar_ventilation=alveolar_ventilation,
        frc=frc, cardiac_output=cardiac_output,
    )
    return WashinResult(
        agent_id=agent_id, blood_gas=lam, plateau=plateau, tau_min=tau,
        t_min=t_min, fa_fi=float(fa_fi[0]), tier=model.tier,
        alveolar_ventilation=alveolar_ventilation, frc=frc, cardiac_output=cardiac_output,
    )


def washin_comparison(ds: Dataset, agents: Optional[List[str]] = None, **kwargs) -> List[WashinResult]:
    """Wash-in for every volatile agent, sorted fastest-first (highest plateau).

    The honest counterpart to the model-divergence view, for inhalational uptake:
    it makes the textbook solubility ordering (desflurane / nitrous oxide fast,
    isoflurane slow) computable from the curated partition coefficients.
    """
    if agents is None:
        agents = [m.id for m in ds if m.purpose == "physicochemical"]
    results = [washin(ds, a, **kwargs) for a in agents]
    return sorted(results, key=lambda r: r.plateau, reverse=True)


def washout(
    ds: Dataset,
    agent_id: str,
    *,
    t_min: float = 3.0,
    alveolar_ventilation: float = 4.0,
    frc: float = 2.5,
    cardiac_output: float = 5.0,
) -> WashoutResult:
    """Single-compartment alveolar wash-out (FA/FA₀) for one volatile agent.

    The offset mirror of :func:`washin`: a comparative, solubility-driven
    characterization of inhalational *emergence* from the curated blood:gas
    partition coefficient. Reports the early elimination floor (= 1 − the wash-in
    plateau), the time constant τ, and FA/FA₀ at ``t_min`` minutes. Lower blood:gas
    solubility → lower floor → faster, more complete wash-out. Standard
    70-kg-adult ventilation constants, all overridable; NOT a per-patient predictor.
    """
    model = ds[agent_id]
    if model.purpose != "physicochemical":
        raise ValueError(f"{agent_id} has purpose '{model.purpose}', expected 'physicochemical'")
    lam = _params(model).get("blood_gas")
    if lam is None:
        raise ValueError(f"{agent_id} has no blood_gas partition coefficient")
    fa_fa0, floor, tau = alveolar_washout(
        lam, np.array([t_min]), alveolar_ventilation=alveolar_ventilation,
        frc=frc, cardiac_output=cardiac_output,
    )
    return WashoutResult(
        agent_id=agent_id, blood_gas=lam, floor=floor, tau_min=tau,
        t_min=t_min, fa_fa0=float(fa_fa0[0]), tier=model.tier,
        alveolar_ventilation=alveolar_ventilation, frc=frc, cardiac_output=cardiac_output,
    )


def washout_comparison(ds: Dataset, agents: Optional[List[str]] = None, **kwargs) -> List[WashoutResult]:
    """Wash-out for every volatile agent, sorted fastest-first (lowest floor).

    The offset counterpart to :func:`washin_comparison`: it makes the textbook
    emergence ordering (desflurane / nitrous oxide fast, isoflurane slow) computable
    from the curated partition coefficients — the clinical reason desflurane is
    preferred for long cases.
    """
    if agents is None:
        agents = [m.id for m in ds if m.purpose == "physicochemical"]
    results = [washout(ds, a, **kwargs) for a in agents]
    return sorted(results, key=lambda r: r.floor)


def mac(
    ds: Dataset,
    agent_id: str,
    *,
    age: float,
    end_tidal_pct: Optional[float] = None,
    n2o_end_tidal_pct: Optional[float] = None,
    n2o_agent_id: str = "volatiles.nitrous_oxide.mac",
) -> MacResult:
    """Evaluate age-corrected MAC for a volatile agent.

    If ``end_tidal_pct`` is given, also returns the MAC fraction (a depth
    surrogate). If ``n2o_end_tidal_pct`` is given, adds nitrous oxide's MAC
    fraction (MAC fractions are additive).
    """
    from .simulate import evaluate_safety  # local import to avoid cycle

    model = ds[agent_id]
    if model.purpose != "physicochemical":
        raise ValueError(f"{agent_id} has purpose '{model.purpose}', expected 'physicochemical'")
    p = _params(model)
    mac40 = p["MAC40"]
    mac_age = mac_age_corrected(mac40, age)
    mac_awake_age = mac_age_corrected(p.get("MAC_awake", 0.34 * mac40), age)

    tier_floor, warnings, _ = evaluate_safety(model, {"age": age})
    tier = worst_tier([model.tier, tier_floor])

    res = MacResult(
        agent_id=agent_id, age=age, mac40=mac40, mac_age=mac_age,
        mac_awake_age=mac_awake_age, blood_gas=p.get("blood_gas", float("nan")),
        oil_gas=p.get("oil_gas", float("nan")), tier=tier, warnings=list(warnings),
    )

    if end_tidal_pct is not None:
        res.end_tidal_pct = end_tidal_pct
        res.mac_fraction = _mac_fraction(end_tidal_pct, mac40, age)
        total = res.mac_fraction
        if n2o_end_tidal_pct is not None:
            n2o = ds[n2o_agent_id]
            n2o_mac40 = _params(n2o)["MAC40"]
            total += _mac_fraction(n2o_end_tidal_pct, n2o_mac40, age)
            # composed estimate inherits the worst tier of the two agents
            res.tier = worst_tier([res.tier, n2o.tier])
        res.combined_mac_fraction = total

    return res
