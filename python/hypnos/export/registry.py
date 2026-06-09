"""Model registry: the executable binding between a dataset record and its
covariate equations.

Each PK model id maps (via ``kernel.function`` in the record) to a function that
turns a patient covariate dict into a :class:`~hypnos.reference.MicroParams`.
The verbatim ``covariate_model`` strings in the dataset are the *human-readable*
provenance of these equations; the functions here are their tested
implementation. Validation (``hypnos validate``) asserts every implemented
kernel binding resolves to a function registered in :data:`KERNELS`.
"""
from __future__ import annotations

import re
from typing import Callable, Dict

from ..reference import MicroParams, lbm_james

# Covariate-function signature: dict(patient) -> MicroParams
KernelFn = Callable[[dict], MicroParams]


def _req(patient: dict, key: str) -> float:
    if key not in patient or patient[key] is None:
        raise ValueError(f"model requires patient covariate '{key}'")
    return float(patient[key])


# --------------------------------------------------------------------------- #
# Propofol — Marsh (1991): weight-proportional V1, fixed micro-rate constants
# --------------------------------------------------------------------------- #
def propofol_marsh_1991(patient: dict) -> MicroParams:
    wgt = _req(patient, "weight")
    V1 = 0.228 * wgt
    return MicroParams(
        V1=V1,
        k10=0.119,
        k12=0.112,
        k21=0.055,
        k13=0.0419,
        k31=0.0033,
        ke0=0.26,
        n_compartments=3,
        derived={
            "V2": V1 * 0.112 / 0.055,
            "V3": V1 * 0.0419 / 0.0033,
        },
    )


# --------------------------------------------------------------------------- #
# Propofol — Schnider (1998): volumes/clearances with age/weight/height/LBM
# --------------------------------------------------------------------------- #
def propofol_schnider_1998(patient: dict) -> MicroParams:
    age = _req(patient, "age")
    wgt = _req(patient, "weight")
    hgt = _req(patient, "height")
    sex = patient.get("sex", "M")
    lbm = lbm_james(wgt, hgt, sex)

    V1 = 4.27
    V2 = 18.9 - 0.391 * (age - 53.0)
    V3 = 238.0
    Cl1 = 1.89 + 0.0456 * (wgt - 77.0) - 0.0681 * (lbm - 59.0) + 0.0264 * (hgt - 177.0)
    Cl2 = 1.29 - 0.024 * (age - 53.0)
    Cl3 = 0.836
    ke0 = 0.456
    return MicroParams.from_volumes_clearances(
        V1=V1, Cl1=Cl1, V2=V2, Cl2=Cl2, V3=V3, Cl3=Cl3, ke0=ke0
    )


# --------------------------------------------------------------------------- #
# Remifentanil — Minto (1997): volumes/clearances with age + LBM (James 1976)
# --------------------------------------------------------------------------- #
def remifentanil_minto_1997(patient: dict) -> MicroParams:
    age = _req(patient, "age")
    wgt = _req(patient, "weight")
    hgt = _req(patient, "height")
    sex = patient.get("sex", "M")
    lbm = lbm_james(wgt, hgt, sex)

    V1 = 5.1 - 0.0201 * (age - 40.0) + 0.072 * (lbm - 55.0)
    V2 = 9.82 - 0.0811 * (age - 40.0) + 0.108 * (lbm - 55.0)
    V3 = 5.42
    Cl1 = 2.6 - 0.0162 * (age - 40.0) + 0.0191 * (lbm - 55.0)
    Cl2 = 2.05 - 0.0301 * (age - 40.0)
    Cl3 = 0.076 - 0.00113 * (age - 40.0)
    ke0 = 0.595 - 0.007 * (age - 40.0)
    return MicroParams.from_volumes_clearances(
        V1=V1, Cl1=Cl1, V2=V2, Cl2=Cl2, V3=V3, Cl3=Cl3, ke0=ke0
    )


KERNELS: Dict[str, KernelFn] = {
    "propofol_marsh_1991": propofol_marsh_1991,
    "propofol_schnider_1998": propofol_schnider_1998,
    "remifentanil_minto_1997": remifentanil_minto_1997,
}


def instantiate(model, patient: dict) -> MicroParams:
    """Compute numeric :class:`MicroParams` for a model + virtual patient."""
    if not model.kernel_implemented or model.kernel_function not in KERNELS:
        raise NotImplementedError(f"{model.id}: no implemented PK kernel to instantiate")
    return KERNELS[model.kernel_function](patient)


# --------------------------------------------------------------------------- #
# PD kernels: (effect-site concentration, parameter dict) -> effect
# --------------------------------------------------------------------------- #
def _sigmoid_params(model) -> dict:
    return {p.symbol: p.central for p in model.parameters}


PD_KERNELS = {"propofol_bis_sigmoid"}

# Two-drug response-surface kernels (purpose == "interaction").
INTERACTION_KERNELS = {"greco_response_surface"}


# --------------------------------------------------------------------------- #
# Dose-schedule parsing (mg / mcg, per-kg, per-h / per-min)
# --------------------------------------------------------------------------- #
_AMOUNT_RE = re.compile(
    r"^\s*([0-9]*\.?[0-9]+)\s*(mg|mcg|ug|µg)\s*(?:/\s*kg)?\s*$", re.IGNORECASE
)
_RATE_RE = re.compile(
    r"^\s*([0-9]*\.?[0-9]+)\s*(mg|mcg|ug|µg)\s*(/\s*kg)?\s*/\s*(h|hr|min)\s*$",
    re.IGNORECASE,
)


def _mass_to_mg(value: float, unit: str) -> float:
    unit = unit.lower()
    if unit == "mg":
        return value
    return value / 1000.0  # mcg/ug/µg -> mg


def parse_amount(spec: str, weight: float) -> float:
    """Parse a bolus amount like '2 mg/kg' or '100 mg' into milligrams."""
    m = _AMOUNT_RE.match(spec)
    if not m:
        raise ValueError(f"cannot parse bolus amount {spec!r}")
    val = float(m.group(1))
    mg = _mass_to_mg(val, m.group(2))
    if "/kg" in spec.replace(" ", "").lower():
        mg *= weight
    return mg


def parse_rate(spec: str, weight: float) -> float:
    """Parse an infusion rate like '6 mg/kg/h' or '10 mg/min' into mg/min."""
    m = _RATE_RE.match(spec)
    if not m:
        raise ValueError(f"cannot parse infusion rate {spec!r}")
    val = float(m.group(1))
    mg = _mass_to_mg(val, m.group(2))
    if m.group(3):  # per kg
        mg *= weight
    per = m.group(4).lower()
    if per in ("h", "hr"):
        return mg / 60.0
    return mg  # per min
