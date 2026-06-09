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

import math
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


# --------------------------------------------------------------------------- #
# Propofol — Paedfusor (2005): pediatric, weight-proportional V1, scaled k10
# --------------------------------------------------------------------------- #
def propofol_paedfusor_2005(patient: dict) -> MicroParams:
    wgt = _req(patient, "weight")
    V1 = 0.4584 * wgt
    k10 = 0.1527 * wgt ** (-0.3)
    return MicroParams(
        V1=V1,
        k10=k10,
        k12=0.114,
        k21=0.055,
        k13=0.0419,
        k31=0.0033,
        ke0=0.0,
        n_compartments=3,
        derived={"V2": V1 * 0.114 / 0.055, "V3": V1 * 0.0419 / 0.0033},
    )


# --------------------------------------------------------------------------- #
# Dexmedetomidine — Hannivoort (2015): allometric three-compartment
# --------------------------------------------------------------------------- #
def dexmedetomidine_hannivoort_2015(patient: dict) -> MicroParams:
    wgt = _req(patient, "weight")
    wr = wgt / 70.0
    V1 = 1.78 * wr
    V2 = 30.3 * wr
    V3 = 52.0 * wr
    Cl1 = 0.686 * wr ** 0.75
    Cl2 = 2.98 * wr ** 0.75
    Cl3 = 0.602 * wr ** 0.75
    return MicroParams.from_volumes_clearances(
        V1=V1, Cl1=Cl1, V2=V2, Cl2=Cl2, V3=V3, Cl3=Cl3, ke0=0.0
    )


# --------------------------------------------------------------------------- #
# Propofol — Eleveld (2018): general-purpose model, broad covariate structure
# --------------------------------------------------------------------------- #
def propofol_eleveld_2018(patient: dict) -> MicroParams:
    """Eleveld 2018 general-purpose propofol PK model (PK arm, arterial sampling).

    Transcribed from the published equations (cross-checked against the ``tci`` R
    package, ``pkmod_eleveld_ppf``) and validated to reproduce the reference
    individual (35 y, 70 kg, 170 cm, male, no concomitant anaesthetics):
    V1=6.28, V2=25.5, V3=273, CL=1.79, Q3=1.11, ke0=0.146. The model record's
    ``review_status`` remains ``unverified`` pending human PDF confirmation.
    """
    age = _req(patient, "age")
    wgt = _req(patient, "weight")
    hgt = _req(patient, "height")
    male = 1.0 if str(patient.get("sex", "M")).upper().startswith("M") else 0.0
    opiate = bool(patient.get("opiate_coadministration", False))
    arterial = bool(patient.get("arterial", True))
    pma = patient.get("postmenstrual_age")
    if pma is None:
        pma = age + 40.0 / 52.0

    # fixed effects (theta), 1-indexed to match the source for readability
    th = [None, 6.28, 25.5, 273.0, 1.79, 1.75, 1.11, 0.191, 42.3, 9.06, -0.0156,
          -0.00286, 33.6, -0.0138, 68.3, 2.10, 1.30, 1.42, 0.68]
    ke0_arterial, ke0_venous = 0.146, 1.24

    AGEref, WGTref, HGTref = 35.0, 70.0, 170.0
    PMAref = AGEref + 40.0 / 52.0
    bmi = 10000.0 * wgt / (hgt * hgt)
    bmiref = 10000.0 * WGTref / (HGTref * HGTref)

    def faging(x: float) -> float:
        return math.exp(x * (age - AGEref))

    def fsig(x: float, e50: float, lam: float) -> float:
        return x ** lam / (x ** lam + e50 ** lam)

    def fcentral(x: float) -> float:
        return fsig(x, th[12], 1.0)

    def fopiates(x: float) -> float:
        return math.exp(x * age) if opiate else 1.0

    fcl_mat = fsig(pma * 52.0, th[8], th[9])
    fcl_mat_ref = fsig(PMAref * 52.0, th[8], th[9])
    fq3_mat = fsig(age * 52.0 + 40.0, th[14], 1.0)
    fq3_mat_ref = fsig(AGEref * 52.0 + 40.0, th[14], 1.0)

    if male:
        ffm = (0.88 + (1 - 0.88) / (1 + (age / 13.4) ** (-12.7))) * (9270 * wgt) / (6680 + 216 * bmi)
    else:
        ffm = (1.11 + (1 - 1.11) / (1 + (age / 7.1) ** (-1.1))) * (9270 * wgt) / (8780 + 244 * bmi)
    ffm_ref = (0.88 + (1 - 0.88) / (1 + (AGEref / 13.4) ** (-12.7))) * (9270 * WGTref) / (6680 + 216 * bmiref)

    v1_art = th[1] * fcentral(wgt) / fcentral(WGTref)
    V1 = v1_art if arterial else v1_art * (1 + th[17] * (1 - fcentral(wgt)))
    V2 = th[2] * wgt / WGTref * faging(th[10])
    V3 = th[3] * (ffm / ffm_ref) * fopiates(th[13])
    CL = (male * th[4] + (1 - male) * th[15]) * (wgt / WGTref) ** 0.75 * (fcl_mat / fcl_mat_ref) * fopiates(th[11])
    q2_art = th[5] * (V2 / th[2]) ** 0.75 * (1 + th[16] * (1 - fq3_mat))
    Q2 = q2_art if arterial else q2_art * th[18]
    Q3 = th[6] * (V3 / th[3]) ** 0.75 * (fq3_mat / fq3_mat_ref)
    ke0 = (ke0_arterial if arterial else ke0_venous) * (wgt / 70.0) ** (-0.25)

    return MicroParams.from_volumes_clearances(
        V1=V1, Cl1=CL, V2=V2, Cl2=Q2, V3=V3, Cl3=Q3, ke0=ke0
    )


# --------------------------------------------------------------------------- #
# Remifentanil — Eleveld (2017): allometric general-purpose model
# --------------------------------------------------------------------------- #
def remifentanil_eleveld_2017(patient: dict) -> MicroParams:
    """Eleveld 2017 allometric remifentanil PK model.

    Transcribed from the published equations (cross-checked against the ``tci`` R
    package, ``pkmod_eleveld_remi``) and validated to reproduce the reference
    individual (35 y, 70 kg, 170 cm, male): V1=5.81, V2=8.82, V3=5.03, CL=2.58,
    Q2=1.72, Q3=0.124, ke0=1.09.

    NOTE: the ``tci`` source computes V3 from ``V2ref`` (8.82), which does not
    reproduce the published V3 reference (5.03); that is an apparent transcription
    typo in ``tci``. Hypnos uses ``V3ref`` (5.03), which matches the source paper's
    reference individual. ``review_status`` remains ``unverified`` pending human
    PDF confirmation.
    """
    age = _req(patient, "age")
    wgt = _req(patient, "weight")
    hgt = _req(patient, "height")
    male = 1.0 if str(patient.get("sex", "M")).upper().startswith("M") else 0.0

    AGEref, TBWref, HGTref = 35.0, 70.0, 170.0
    bmi = 10000.0 * wgt / (hgt * hgt)
    bmiref = 10000.0 * TBWref / (HGTref * HGTref)

    if male:
        ffm = (0.88 + (1 - 0.88) / (1 + (age / 13.4) ** (-12.7))) * (9270 * wgt) / (6680 + 216 * bmi)
    else:
        ffm = (1.11 + (1 - 1.11) / (1 + (age / 7.1) ** (-1.1))) * (9270 * wgt) / (8780 + 244 * bmi)
    ffm_ref = (0.88 + (1 - 0.88) / (1 + (AGEref / 13.4) ** (-12.7))) * (9270 * TBWref) / (6680 + 216 * bmiref)

    def faging(x: float) -> float:
        return math.exp(x * (age - AGEref))

    def fsig(x: float, e50: float, lam: float) -> float:
        return x ** lam / (x ** lam + e50 ** lam)

    th = [None, 2.88, -0.00554, -0.00327, -0.0315, 0.470, -0.0260]
    V1ref, V2ref, V3ref, CLref, Q2ref, Q3ref = 5.81, 8.82, 5.03, 2.58, 1.72, 0.124

    size = ffm / ffm_ref
    kmat = fsig(wgt, th[1], 2.0)
    kmat_ref = fsig(TBWref, th[1], 2.0)
    ksex = 1.0 if male else (1 + th[5] * fsig(age, 12.0, 6.0) * (1 - fsig(age, 45.0, 6.0)))

    V1 = V1ref * size * faging(th[2])
    V2 = V2ref * size * faging(th[3]) * ksex
    V3 = V3ref * size * faging(th[4]) * math.exp(th[6] * (wgt - 70.0))
    CL = CLref * size ** 0.75 * (kmat / kmat_ref) * ksex * faging(th[3])
    Q2 = Q2ref * (V2 / V2ref) ** 0.75 * faging(th[2]) * ksex
    Q3 = Q3ref * (V3 / V3ref) ** 0.75 * faging(th[2])
    ke0 = 1.09 * faging(-0.0289)

    return MicroParams.from_volumes_clearances(
        V1=V1, Cl1=CL, V2=V2, Cl2=Q2, V3=V3, Cl3=Q3, ke0=ke0
    )


KERNELS: Dict[str, KernelFn] = {
    "propofol_marsh_1991": propofol_marsh_1991,
    "propofol_schnider_1998": propofol_schnider_1998,
    "propofol_paedfusor_2005": propofol_paedfusor_2005,
    "propofol_eleveld_2018": propofol_eleveld_2018,
    "remifentanil_minto_1997": remifentanil_minto_1997,
    "remifentanil_eleveld_2017": remifentanil_eleveld_2017,
    "dexmedetomidine_hannivoort_2015": dexmedetomidine_hannivoort_2015,
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


# Single-drug effect-site sigmoid kernels (purpose == "pd"); all use sigmoid_emax
# with parameters {E0, Emax, Ce50, gamma} read from the record.
PD_KERNELS = {"propofol_bis_sigmoid", "nmb_tof_sigmoid"}

# Two-drug response-surface kernels (purpose == "interaction").
INTERACTION_KERNELS = {"greco_response_surface"}

# Inhalational physicochemical kernels (purpose == "physicochemical"); MAC age
# correction + partition coefficients, read directly from the record.
VOLATILE_KERNELS = {"volatile_mac"}


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
