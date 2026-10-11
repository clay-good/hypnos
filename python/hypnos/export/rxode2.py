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
from ._common import instantiated_standing, resolve_patient, safe_name
from ._covariate import function_covariate_block
from ._developmental import function_block as developmental_functions
from ._estimation import (
    covariance_step_note, independence_caveat, theta_covariance, theta_estimates,
)
from ._variability import omega_diagonal, residual_spec
from .registry import KERNELS


def _comment_block(text: str) -> str:
    return "\n".join("# " + line for line in text.splitlines())


def _residual_endpoint(spec) -> str:
    """rxode2 residual endpoint for ``cp`` from a normalized :class:`ResidualSpec`."""
    if spec.model == "log":
        return f"  cp ~ lnorm({spec.log_sd:.6g})           # residual: {spec.label}"
    if spec.model == "proportional":
        return f"  cp ~ prop({spec.prop_var ** 0.5:.6g})           # residual: {spec.label}"
    if spec.model == "additive":
        return f"  cp ~ add({spec.add_sd:.6g})           # residual: {spec.label}"
    # combined
    return (f"  cp ~ prop({spec.prop_var ** 0.5:.6g}) + add({spec.add_sd:.6g})"
            f"   # residual: {spec.label}")


def _population_block(model, name: str, p) -> list:
    """Native nlmixr2/rxode2 population model: V/Cl with log-normal η + Σ + the Ω matrix.

    Emitted only for models that publish between-subject variability (spec §5's
    never-synthesize rule lives upstream in ``has_published_variability``). The
    structural block above stays the typical-value reference; this is the runnable
    population companion — ``rxSolve(<name>_pop, ev, omega=<name>_omega, nSub=...)``.
    """
    vc = p.as_volumes_clearances()
    omegas = {s: om2 for s, om2, _ in omega_diagonal(model)}
    lines = [
        "",
        "# ── population variability (v0.2): between-subject Omega + residual Sigma ──",
        f"#    band-tier {model.band_tier} (worst of structural/variability/residual; spec §5)",
        "#    log-normal BSV on volumes/clearances:  P_i = P_typ * exp(eta),  eta ~ N(0, omega2)",
        f"{name}_pop <- rxode2({{",
    ]
    # structural parameters in V/Cl form, each carrying its η where BSV is published
    for sym in ("Cl1", "V1", "Cl2", "V2", "Cl3", "V3", "ke0"):
        val = vc[sym]
        if sym in omegas:
            lines.append(f"  {sym:<3} <- {val:.10g} * exp(eta.{sym})")
        else:
            lines.append(f"  {sym:<3} <- {val:.10g}")
    lines += [
        "  k10 <- Cl1 / V1",
        "  k12 <- Cl2 / V1",
        "  k21 <- Cl2 / V2",
        "  k13 <- Cl3 / V1",
        "  k31 <- Cl3 / V3",
        "  d/dt(A1) = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3",
        "  d/dt(A2) =  k12 * A1 - k21 * A2",
        "  d/dt(A3) =  k13 * A1 - k31 * A3",
        "  d/dt(Ce) =  ke0 * (A1 / V1 - Ce)",
        "  cp = A1 / V1",
    ]
    spec = residual_spec(model)
    if spec is not None:
        lines.append(_residual_endpoint(spec))
    lines.append("})")
    # the Omega matrix (diagonal here; off-diagonals would extend the lotri block)
    lines.append(f"{name}_omega <- lotri({{")
    for sym, om2, cv in omega_diagonal(model):
        lines.append(f"  eta.{sym} ~ {om2:.6g}        # BSV, CV~{cv:.0f}%")
    lines.append("})")
    return lines


def _lotri_cov(name: str, prefix: str, syms, cov) -> list:
    """A ``lotri`` covariance block: one diagonal line per parameter, or — when any
    off-diagonal is published — a single lower-triangular block over all of them."""
    lines = [f"{name} <- lotri({{"]
    if any(cov[i][j] for i in range(len(syms)) for j in range(i)):
        tri = ", ".join(f"{cov[i][j]:.10g}" for i in range(len(syms)) for j in range(i + 1))
        lines.append("  " + " + ".join(f"{prefix}{s}" for s in syms) + f" ~ c({tri})")
    else:
        lines += [f"  {prefix}{s} ~ {cov[i][i]:.10g}" for i, s in enumerate(syms)]
    lines.append("})")
    return lines


def _estimation_block(model, name: str) -> list:
    """The curated ESTIMATE covariance of θ, usable to seed SIR / parametric resampling of
    the typical values in a simulation script (v0.3 §8). Emitted only when the model
    publishes estimation uncertainty; deliberately NOT named ``omega``."""
    ests = theta_estimates(model)
    if not ests:
        return []
    syms, cov = theta_covariance(model)
    lines = [
        "",
        "# ── estimation uncertainty (v0.3): the estimate covariance of theta ──",
        "#    ESTIMATION covariance (how well the data pin the typical value; reducible) —",
        "#    NOT Omega (between-subject spread). Do not pass this as `omega=`.",
        f"#    hypnos:uncertaintyStatus = {model.uncertainty_status}",
        f"#    hypnos:estimationBandTier = {model.estimation_band_tier}",
        f"#    method {ests[0].method}; natural scale, on the CURATED typical values below",
        "#    (not the patient-instantiated parameters above; RSE% is the scale-free form).",
        f"{name}_theta_ref <- c(" + ", ".join(f"{e.symbol} = {e.central:.10g}" for e in ests) + ")",
        f"{name}_theta_rse <- c(" + ", ".join(
            f"{e.symbol} = {e.rse_percent:.6g}" for e in ests if e.rse_percent is not None) + ")",
    ]
    lines += _lotri_cov(f"{name}_theta_cov", "theta.", syms, cov)
    for note in (covariance_step_note(model), independence_caveat(model)):
        if note:
            lines.append(f"# NOTE: {note}")
    return lines


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    pat = resolve_patient(model, patient)
    tier, env_warns = instantiated_standing(model, ds, pat)
    head = _comment_block(annotate.banner(model, tier))
    if env_warns:
        head += "\n" + _comment_block("\n".join(annotate.envelope_comments(env_warns)))
    # v0.5 §D: reversal coupling + opt-in disease modifiers, ahead of the kernel gate.
    breadth = annotate.breadth_comments(model, ds)
    if breadth:
        head += "\n" + _comment_block("\n".join(breadth))
    if not (model.kernel_implemented and model.kernel_function in KERNELS):
        return head + "\n# KERNEL PENDING — no instantiated rxode2 model emitted.\n"

    p = KERNELS[model.kernel_function](pat)
    prov = annotate.provenance(model, ds, tier)
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
    if model.has_published_variability:
        lines += _population_block(model, safe_name(model), p)
    lines += _estimation_block(model, safe_name(model))
    lines += function_covariate_block(model, ds, safe_name(model), "r")
    lines += developmental_functions(model, safe_name(model), "r")
    return "\n".join(lines) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.rxode2.R"
