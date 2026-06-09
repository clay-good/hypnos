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
from typing import Any, Dict, List, Optional

from .load import Dataset
from .models import worst_tier
from .reference import mac_age_corrected
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


def _params(model) -> Dict[str, float]:
    return {p.symbol: p.central for p in model.parameters}


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
