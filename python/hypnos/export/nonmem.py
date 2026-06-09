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
from ._variability import contiguous_block, residual_spec
from .registry import KERNELS

# $PK / $THETA structural order (TRANS4), paired with the dataset symbol whose
# omega2 supplies the matching ETA's between-subject variance. The sym sequence
# matches _variability.VC_ORDER, so ETA(k) corresponds to VC_ORDER[k-1].
_ETA_ORDER = [
    ("CL", "Cl1"), ("V1", "V1"), ("Q2", "Cl2"), ("V2", "V2"),
    ("Q3", "Cl3"), ("V3", "V3"), ("KE0", "ke0"),
]


def _comment_block(text: str) -> str:
    return "\n".join("; " + line for line in text.splitlines())


def _residual_sigma(model) -> Optional[Dict[str, Any]]:
    """Translate a curated residual_error into $ERROR lines + $SIGMA variances.

    Returns ``{"error": [lines], "sigma": [variances], "label": str}`` or None.
    NONMEM's $SIGMA carries variances, so SD-scale terms are squared here.
    """
    spec = residual_spec(model)
    if spec is None:
        return None
    if spec.model == "log":
        return {"error": ["  Y = IPRED*EXP(EPS(1))  ; log-additive residual"],
                "sigma": [spec.log_sd ** 2], "label": spec.label}
    if spec.model == "proportional":
        return {"error": ["  Y = IPRED*(1 + EPS(1))  ; proportional residual"],
                "sigma": [spec.prop_var], "label": spec.label}
    if spec.model == "additive":
        return {"error": ["  Y = IPRED + EPS(1)  ; additive residual"],
                "sigma": [spec.add_sd ** 2], "label": spec.label}
    if spec.model == "combined":
        return {"error": ["  Y = IPRED*(1 + EPS(1)) + EPS(2)  ; combined residual"],
                "sigma": [spec.prop_var, spec.add_sd ** 2], "label": spec.label}
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
        block = contiguous_block(model)
        block_n = len(block[0]) if block else 0
        if block:
            syms, cov = block
            lines.append(f"$OMEGA BLOCK({block_n})  ; correlated BSV (off-diagonal Omega, "
                         "covariance scale = r*sqrt(om_i*om_j))")
            for i, sym in enumerate(syms):
                nm = _ETA_ORDER[i][0]
                row = "  ".join(f"{cov[i][j]:.6g}" for j in range(i + 1))
                lines.append(f"  {row}       ; ETA({i + 1}) {nm}")
        if block_n < len(_ETA_ORDER):
            lines.append("$OMEGA  ; diagonal BSV" if block else "$OMEGA")
        for k, (nm, sym) in enumerate(_ETA_ORDER, start=1):
            if k <= block_n:
                continue  # already carried by the $OMEGA BLOCK above
            om2 = omegas.get(sym)
            if om2 is not None:
                cv = 100.0 * math.sqrt(math.exp(om2) - 1.0)
                lines.append(f"  {om2:.6g}       ; ETA({k}) {nm}  (BSV, CV~{cv:.0f}%)")
            else:
                lines.append(f"  0 FIX        ; ETA({k}) {nm}  (no published BSV — fixed)")
        if model.omega_block is not None and block is None:
            lines.append("; NOTE: off-diagonal Omega correlations are curated but not emitted as a "
                         "$OMEGA BLOCK here (incomplete or non-contiguous); see omega_block in the "
                         "source record. The diagonal above assumes independence (with that caveat).")
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
