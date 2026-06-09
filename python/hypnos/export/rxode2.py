"""nlmixr2 / rxode2 (R) exporter.

Emits an ``rxode2()`` model block instantiated at a reference individual, with
the three-compartment mammillary ODEs + the effect-site link written out
explicitly. The micro-rate constants are emitted as named assignments (and
echoed in a machine-readable ``# hypnos.params:`` line) so the export round-trips
against the reference kernel.
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
        return head + "\n# KERNEL PENDING — no instantiated rxode2 model emitted.\n"

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
        "library(rxode2)",
        "",
        f"# {model.label}",
        f"{safe_name(model)} <- rxode2({{",
        f"  V1  <- {p.V1:.10g}",
        f"  k10 <- {p.k10:.10g}",
        f"  k12 <- {p.k12:.10g}",
        f"  k21 <- {p.k21:.10g}",
        f"  k13 <- {p.k13:.10g}",
        f"  k31 <- {p.k31:.10g}",
        f"  ke0 <- {p.ke0:.10g}",
        "  d/dt(A1) = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3",
        "  d/dt(A2) =  k12 * A1 - k21 * A2",
        "  d/dt(A3) =  k13 * A1 - k31 * A3",
        "  d/dt(Ce) =  ke0 * (A1 / V1 - Ce)",
        "  cp = A1 / V1",
        "})",
    ]
    return "\n".join(lines) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.rxode2.R"
