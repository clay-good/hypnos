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
class DecrementTime:
    model_id: str
    infusion: str
    duration_min: float
    fraction: float
    conc_at_stop: float          # plasma conc at the moment infusion stops (internal ug/mL)
    decrement_min: float         # time after stop for plasma to fall by `fraction` (inf if not reached)
    tier: str
    warnings: List[str] = field(default_factory=list)


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


def decrement_time(
    ds: Dataset, model_id: str, *, patient: Dict[str, Any], infusion: str,
    duration: float, fraction: float = 0.5,
) -> DecrementTime:
    """Plasma decrement time after a **constant-rate** infusion of given duration.

    Runs a fixed-rate infusion for ``duration`` minutes, stops it, and returns the
    time for the plasma concentration to fall by ``fraction`` (default 50%) from
    its value at the stop. Forward-only and a function of the *fixed dose history*.

    This is **not** the classic context-sensitive half-time, which is defined for a
    target-controlled infusion holding plasma *constant* — that requires inverse
    control, which Hypnos does not do (spec §10). This metric captures the same
    context-sensitivity qualitatively (it lengthens with infusion duration for
    accumulating drugs, and is near-flat for remifentanil) from a forward,
    constant-rate regimen.
    """
    from .export.registry import instantiate
    from .simulate import build_dosing, evaluate_safety

    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must be in (0, 1)")
    model = ds[model_id]
    if model.purpose != "pk":
        raise ValueError(f"decrement_time expects a PK model; {model_id} is '{model.purpose}'")
    params = instantiate(model, patient)  # raises NotImplementedError if kernel pending

    weight = float(patient.get("weight", 70.0))
    schedule = [("infusion", 0.0, infusion), ("infusion", float(duration), "0 mg/kg/h")]
    dosing = build_dosing(schedule, weight)

    horizon = duration + max(120.0, 6.0 * duration)
    grid = np.linspace(0.0, horizon, int(horizon * 4) + 1)
    t = np.union1d(grid, [float(duration)])  # ensure the stop instant is sampled
    traj = _simulate_ref(params, dosing, t)

    i_stop = int(np.searchsorted(t, duration))
    c_stop = float(traj.cp[i_stop])
    tier_floor, warnings, _ = evaluate_safety(model, patient)
    warnings = list(warnings)

    decrement = float("inf")
    if c_stop <= 0:
        warnings.append("WARNING: plasma concentration at stop is ~0; decrement undefined")
    else:
        threshold = (1.0 - fraction) * c_stop
        after = np.where((t >= duration) & (traj.cp <= threshold))[0]
        if len(after):
            decrement = float(t[after[0]] - duration)
        else:
            warnings.append(f"WARNING: {100*fraction:g}% decrement not reached within {horizon:g} min")

    return DecrementTime(
        model_id=model_id, infusion=infusion, duration_min=float(duration),
        fraction=fraction, conc_at_stop=c_stop, decrement_min=decrement,
        tier=worst_tier([model.tier, tier_floor]), warnings=warnings,
    )
