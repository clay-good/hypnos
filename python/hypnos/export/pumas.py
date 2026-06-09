"""Pumas (Julia) exporter.

Emits a Pumas ``@model`` block instantiated at a reference individual: the
micro-rate constants in ``@pre``, the three-compartment + effect-site ODEs in
``@dynamics``, and ``cp`` in ``@derived``. The constants are echoed in a
machine-readable ``# hypnos.params:`` comment so the export round-trips against
the reference kernel.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from . import annotate
from ._common import resolve_patient, safe_name
from .registry import KERNELS


def _comment_block(text: str) -> str:
    return "\n".join("# " + line for line in text.splitlines())


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    pat = resolve_patient(model, patient)
    head = _comment_block(annotate.banner(model, model.tier))
    if not (model.kernel_implemented and model.kernel_function in KERNELS):
        return head + "\n# KERNEL PENDING — no instantiated Pumas model emitted.\n"

    p = KERNELS[model.kernel_function](pat)
    prov = annotate.provenance(model, ds, model.tier)
    params = (
        f"k10={p.k10:.10g} k12={p.k12:.10g} k21={p.k21:.10g} "
        f"k13={p.k13:.10g} k31={p.k31:.10g} ke0={p.ke0:.10g} V1={p.V1:.10g}"
    )
    lines = [
        head,
        f"# instantiated for: {pat}",
        f"# hypnos.params: {params}",
        f"# hypnos:clinicalUse = {prov['hypnos:clinicalUse']}",
    ]
    for u in prov["bqmodel:isDerivedFrom"]:
        lines.append(f"# bqmodel:isDerivedFrom = {u}")
    lines += [
        "using Pumas",
        "",
        f"# {model.label}",
        f"{safe_name(model)} = @model begin",
        "    @pre begin",
        f"        V1  = {p.V1:.10g}",
        f"        k10 = {p.k10:.10g}",
        f"        k12 = {p.k12:.10g}",
        f"        k21 = {p.k21:.10g}",
        f"        k13 = {p.k13:.10g}",
        f"        k31 = {p.k31:.10g}",
        f"        ke0 = {p.ke0:.10g}",
        "    end",
        "    @dynamics begin",
        "        A1' = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3",
        "        A2' =  k12 * A1 - k21 * A2",
        "        A3' =  k13 * A1 - k31 * A3",
        "        Ce' =  ke0 * (A1 / V1 - Ce)",
        "    end",
        "    @derived begin",
        "        cp = A1 / V1",
        "    end",
        "end",
    ]
    return "\n".join(lines) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.pumas.jl"
