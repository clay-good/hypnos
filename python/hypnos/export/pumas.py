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
from ._variability import omega_diagonal, residual_spec
from .registry import KERNELS


def _comment_block(text: str) -> str:
    return "\n".join("# " + line for line in text.splitlines())


def _residual_derived(spec) -> str:
    """Pumas ``@derived`` dv distribution for ``cp`` from a normalized ResidualSpec."""
    if spec.model == "log":
        return (f"        dv ~ @. LogNormal(log(cp), {spec.log_sd:.6g})"
                f"   # residual: {spec.label}")
    if spec.model == "proportional":
        return (f"        dv ~ @. Normal(cp, cp*{spec.prop_var ** 0.5:.6g})"
                f"   # residual: {spec.label}")
    if spec.model == "additive":
        return f"        dv ~ @. Normal(cp, {spec.add_sd:.6g})   # residual: {spec.label}"
    # combined
    return (f"        dv ~ @. Normal(cp, sqrt((cp*{spec.prop_var ** 0.5:.6g})^2 "
            f"+ {spec.add_sd:.6g}^2))   # residual: {spec.label}")


def _population_block(model, name: str, p) -> list:
    """Native Pumas population model: @param Ω, @random η, V/Cl with log-normal BSV, Σ.

    Emitted only when the model publishes BSV. The structural ``@model`` above stays
    the typical-value reference; this is the runnable NLME companion.
    """
    vc = p.as_volumes_clearances()
    diag = omega_diagonal(model)
    omega_syms = {s for s, _, _ in diag}
    lines = [
        "",
        "# ── population variability (v0.2): @random Omega + residual Sigma ──",
        f"#    band-tier {model.band_tier} (worst of structural/variability/residual; spec §5)",
        f"{name}_pop = @model begin",
        "    @param begin",
    ]
    for sym, om2, cv in diag:
        lines.append(f"        ω²_{sym} ∈ ConstDomain({om2:.6g})        # BSV variance, CV~{cv:.0f}%")
    lines.append("    end")
    lines.append("    @random begin")
    for sym, om2, _ in diag:
        lines.append(f"        η_{sym} ~ Normal(0.0, sqrt(ω²_{sym}))")
    lines.append("    end")
    lines.append("    @pre begin")
    for sym in ("Cl1", "V1", "Cl2", "V2", "Cl3", "V3", "ke0"):
        val = vc[sym]
        if sym in omega_syms:
            lines.append(f"        {sym} = {val:.10g} * exp(η_{sym})")
        else:
            lines.append(f"        {sym} = {val:.10g}")
    lines += [
        "        k10 = Cl1 / V1",
        "        k12 = Cl2 / V1",
        "        k21 = Cl2 / V2",
        "        k13 = Cl3 / V1",
        "        k31 = Cl3 / V3",
        "    end",
        "    @dynamics begin",
        "        A1' = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3",
        "        A2' =  k12 * A1 - k21 * A2",
        "        A3' =  k13 * A1 - k31 * A3",
        "        Ce' =  ke0 * (A1 / V1 - Ce)",
        "    end",
        "    @derived begin",
        "        cp = @. A1 / V1",
    ]
    spec = residual_spec(model)
    if spec is not None:
        lines.append(_residual_derived(spec))
    lines += ["    end", "end"]
    return lines


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
    if model.has_published_variability:
        lines += _population_block(model, safe_name(model), p)
    return "\n".join(lines) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.pumas.jl"
