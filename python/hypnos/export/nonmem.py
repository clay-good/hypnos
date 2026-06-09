"""NONMEM control-stream exporter.

Emits an ADVAN11/TRANS4 three-compartment control stream instantiated at a
reference individual. The verbatim population covariate equations are preserved
as comments, and the Hypnos banner + provenance ride along as ``;`` comment lines
so they survive any text round-trip.

When the model carries a curated random-effects layer (v0.2), the exporter emits
the real ``$OMEGA`` diagonal (from each parameter's ``omega2``) and ``$SIGMA``
(from ``residual_error``), wiring the matching ``EXP(ETA(.))`` into ``$PK`` — the
single most natural upgrade from v0.1's ``$OMEGA 0 FIX`` placeholder (spec §8).
Otherwise it keeps ``0 FIX`` and names the missing component in a comment.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Optional

from . import annotate
from ._common import resolve_patient, safe_name
from .registry import KERNELS

# $PK / $THETA structural order (TRANS4), paired with the dataset symbol whose
# omega2 supplies the matching ETA's between-subject variance.
_ETA_ORDER = [
    ("CL", "Cl1"), ("V1", "V1"), ("Q2", "Cl2"), ("V2", "V2"),
    ("Q3", "Cl3"), ("V3", "V3"), ("KE0", "ke0"),
]


def _comment_block(text: str) -> str:
    return "\n".join("; " + line for line in text.splitlines())


def _residual_sigma(model) -> Optional[Dict[str, Any]]:
    """Translate a curated residual_error into $ERROR lines + $SIGMA variances.

    Returns ``{"error": [lines], "sigma": [variances], "label": str}`` or None.
    """
    re_ = model.residual_error
    if re_ is None:
        return None
    m = re_.model
    if m == "log":
        sd = (re_.log or {}).get("sd")
        var = (re_.log or {}).get("variance")
        v = var if var is not None else (sd ** 2 if sd is not None else None)
        if v is None:
            return None
        return {"error": ["  Y = IPRED*EXP(EPS(1))  ; log-additive residual"],
                "sigma": [v], "label": "log-additive (proportional on natural scale)"}
    if m == "proportional":
        v = (re_.proportional or {}).get("variance")
        if v is None and (re_.proportional or {}).get("cv_percent") is not None:
            v = (re_.proportional["cv_percent"] / 100.0) ** 2
        if v is None:
            return None
        return {"error": ["  Y = IPRED*(1 + EPS(1))  ; proportional residual"],
                "sigma": [v], "label": "proportional"}
    if m == "additive":
        sd = (re_.additive or {}).get("sd")
        var = (re_.additive or {}).get("variance")
        v = var if var is not None else (sd ** 2 if sd is not None else None)
        if v is None:
            return None
        return {"error": ["  Y = IPRED + EPS(1)  ; additive residual"],
                "sigma": [v], "label": "additive"}
    if m == "combined":
        pv = (re_.proportional or {}).get("variance")
        if pv is None and (re_.proportional or {}).get("cv_percent") is not None:
            pv = (re_.proportional["cv_percent"] / 100.0) ** 2
        asd = (re_.additive or {}).get("sd")
        av = (re_.additive or {}).get("variance")
        av = av if av is not None else (asd ** 2 if asd is not None else None)
        if pv is None or av is None:
            return None
        return {"error": ["  Y = IPRED*(1 + EPS(1)) + EPS(2)  ; combined residual"],
                "sigma": [pv, av], "label": "combined (proportional + additive)"}
    return None


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    pat = resolve_patient(model, patient)
    tier = model.tier
    lines = []
    lines.append(_comment_block(annotate.banner(model, tier)))
    lines.append(";")
    lines.append(f"; Instantiated for reference patient: {pat}")
    lines.append("; Population covariate equations (verbatim from source):")
    for p in model.parameters:
        if p.covariate_model:
            lines.append(f";   {p.symbol}: {p.covariate_model}")
    lines.append(";")

    if not (model.kernel_implemented and model.kernel_function in KERNELS):
        lines.append("; KERNEL PENDING — numeric THETA block omitted (no verified kernel).")
        return "\n".join(lines) + "\n"

    vc = KERNELS[model.kernel_function](pat).as_volumes_clearances()
    omegas = model.bsv_omegas() if model.has_published_variability else {}
    sigma = _residual_sigma(model) if model.has_published_variability else None
    theta_n = {"CL": 1, "V1": 2, "Q2": 3, "V2": 4, "Q3": 5, "V3": 6, "KE0": 7}

    lines.append(f"$PROBLEM {model.id} (Hypnos export — NOT FOR CLINICAL USE)")
    lines.append("$INPUT ID TIME AMT RATE DV CMT EVID")
    lines.append("$DATA data.csv IGNORE=@")
    lines.append("$SUBROUTINES ADVAN11 TRANS4")
    lines.append("$PK")
    # Wire EXP(ETA(.)) onto each structural parameter that carries a curated omega2;
    # parameters without BSV keep their fixed ETA at 0 (honest "partial" structure).
    for k, (nm, sym) in enumerate(_ETA_ORDER, start=1):
        eta = f"*EXP(ETA({k}))" if omegas else ""
        lines.append(f"  {nm} = THETA({theta_n[nm]}){eta}")
    lines.append("  S1 = V1")
    lines.append("$ERROR")
    lines.append("  IPRED = A(1)/V1")
    if sigma:
        lines.extend(sigma["error"])
    else:
        lines.append("  Y = IPRED")
    # THETA in TRANS4 order: CL V1 Q2 V2 Q3 V3, then KE0
    lines.append("$THETA")
    lines.append(f"  {vc['Cl1']:.6g}   ; 1 CL  (L/min)")
    lines.append(f"  {vc['V1']:.6g}    ; 2 V1  (L)")
    lines.append(f"  {vc['Cl2']:.6g}   ; 3 Q2  (L/min)")
    lines.append(f"  {vc['V2']:.6g}    ; 4 V2  (L)")
    lines.append(f"  {vc['Cl3']:.6g}   ; 5 Q3  (L/min)")
    lines.append(f"  {vc['V3']:.6g}    ; 6 V3  (L)")
    lines.append(f"  {vc['ke0']:.6g}   ; 7 KE0 (1/min)")
    if omegas:
        lines.append(f"; $OMEGA — between-subject variance (eta-scale); band-tier {model.band_tier}")
        lines.append("$OMEGA")
        for k, (nm, sym) in enumerate(_ETA_ORDER, start=1):
            om2 = omegas.get(sym)
            if om2 is not None:
                cv = 100.0 * math.sqrt(math.exp(om2) - 1.0)
                lines.append(f"  {om2:.6g}       ; ETA({k}) {nm}  (BSV, CV~{cv:.0f}%)")
            else:
                lines.append(f"  0 FIX        ; ETA({k}) {nm}  (no published BSV — fixed)")
        if model.omega_block is not None:
            lines.append("; NOTE: off-diagonal Omega correlations are curated but emitted as a "
                         "diagonal here; see omega_block in the source record.")
    else:
        lines.append("$OMEGA 0 FIX  ; no published between-subject variability (variability_status: "
                     f"{model.variability_status})")
    if sigma:
        lines.append(f"$SIGMA  ; residual error — {sigma['label']}")
        for v in sigma["sigma"]:
            lines.append(f"  {v:.6g}")
    else:
        lines.append("$SIGMA 0 FIX  ; residual error not curated")
    lines.append("$ESTIMATION MAXEVAL=0  ; simulation/instantiation only")
    lines.append(";")
    prov = annotate.provenance(model, ds, tier)
    lines.append(f"; hypnos:clinicalUse = {prov['hypnos:clinicalUse']}")
    lines.append(f"; hypnos:confidenceTier = {prov['hypnos:confidenceTier']}")
    for u in prov["bqmodel:isDerivedFrom"]:
        lines.append(f"; bqmodel:isDerivedFrom = {u}")
    return "\n".join(lines) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.ctl"
