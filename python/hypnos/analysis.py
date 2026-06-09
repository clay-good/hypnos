"""Derived PK/PD characterizations (forward-only, safe by construction).

Currently: **time to peak effect** (``tpeak``), the spec's ``effect_link``
"time-to-peak-effect parameterization" (§3). After a bolus, the effect-site
concentration peaks when its rate of change is zero, i.e. exactly when the
effect-site concentration equals the plasma concentration (Ce = Cp). ``tpeak``
characterizes a drug/model's *onset* and is independent of dose magnitude (the
system is linear). It is computed by pure forward simulation.

**Why no context-sensitive half-time (CSHT).** The classic CSHT is the time for
plasma to fall 50% after a target-controlled infusion that held plasma
*constant* for some duration. Computing it requires solving for the infusion that
maintains a target concentration, which is inverse control: exactly the step
Hypnos refuses to take (spec §10). Onset (``tpeak``) and the decline after a
*fixed* dose history are forward problems and in scope; constant-concentration
CSHT is not.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List

import numpy as np

from .load import Dataset
from .models import worst_tier
from .reference import Dosing, simulate as _simulate_ref


@dataclass
class PeakEffect:
    model_id: str
    tpeak_min: float                 # time of the effect-site peak after a bolus
    ke0: float                       # the effect-site equilibration rate constant
    ce_cp_ratio_at_peak: float       # ~1.0 by definition (Ce = Cp at the peak); a sanity check
    tier: str
    warnings: List[str] = field(default_factory=list)


def time_to_peak_effect(
    ds: Dataset, model_id: str, *, patient: Dict[str, Any], tmax: float = 20.0, n: int = 8001
) -> PeakEffect:
    """Time to peak effect-site concentration after a bolus (onset characterization).

    Forward-only and dose-independent. Requires a PK model with an effect
    compartment (ke0 > 0); raises for PK-only models (e.g. Kim remifentanil,
    Paedfusor) where effect-site onset is undefined.
    """
    from .export.registry import KERNELS
    from .simulate import evaluate_safety

    model = ds[model_id]
    if model.purpose != "pk":
        raise ValueError(f"time_to_peak_effect expects a PK model; {model_id} is '{model.purpose}'")
    if not model.kernel_implemented or model.kernel_function not in KERNELS:
        raise NotImplementedError(f"{model_id}: reference kernel not implemented")
    params = KERNELS[model.kernel_function](patient)
    if params.ke0 <= 0:
        raise ValueError(
            f"{model_id} has no effect compartment (ke0 = 0); time-to-peak-effect is undefined"
        )

    t = np.linspace(0.0, tmax, n)
    traj = _simulate_ref(params, Dosing(boluses=((0.0, 1.0),)), t)  # unit bolus; tpeak is dose-independent
    i = int(np.argmax(traj.ce))

    tier_floor, warnings, _ = evaluate_safety(model, patient)
    warnings = list(warnings)
    if i >= n - 1:
        warnings.append(f"WARNING: effect-site peak not reached within tmax={tmax:g} min; increase tmax")
    ratio = float(traj.ce[i] / traj.cp[i]) if traj.cp[i] > 0 else float("nan")

    return PeakEffect(
        model_id=model_id, tpeak_min=float(t[i]), ke0=params.ke0,
        ce_cp_ratio_at_peak=ratio, tier=worst_tier([model.tier, tier_floor]),
        warnings=warnings,
    )
