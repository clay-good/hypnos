"""NONMEM control-stream exporter.

Emits an ADVAN11/TRANS4 three-compartment control stream instantiated at a
reference individual (point parameters; inter-individual variability is out of
scope for a deterministic reference export). The verbatim population covariate
equations are preserved as comments, and the Hypnos banner + provenance ride
along as ``;`` comment lines so they survive any text round-trip.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from . import annotate
from ._common import resolve_patient, safe_name
from .registry import KERNELS


def _comment_block(text: str) -> str:
    return "\n".join("; " + line for line in text.splitlines())


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
    lines.append(f"$PROBLEM {model.id} (Hypnos export — NOT FOR CLINICAL USE)")
    lines.append("$INPUT ID TIME AMT RATE DV CMT EVID")
    lines.append("$DATA data.csv IGNORE=@")
    lines.append("$SUBROUTINES ADVAN11 TRANS4")
    lines.append("$PK")
    lines.append("  CL = THETA(1)")
    lines.append("  V1 = THETA(2)")
    lines.append("  Q2 = THETA(3)")
    lines.append("  V2 = THETA(4)")
    lines.append("  Q3 = THETA(5)")
    lines.append("  V3 = THETA(6)")
    lines.append("  KE0 = THETA(7)")
    lines.append("  S1 = V1")
    lines.append("$ERROR")
    lines.append("  IPRED = A(1)/V1")
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
    lines.append("$OMEGA 0 FIX")
    lines.append("$SIGMA 0 FIX")
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
