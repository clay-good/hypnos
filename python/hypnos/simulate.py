"""Forward simulation API — the safe side of the dataset/simulator boundary.

``simulate`` maps (dose history -> predicted concentration/effect) for one model
and one virtual patient, enforcing the two anesthesia-specific load-bearing
ideas:

* **applicability envelope** — out-of-envelope requests are tiered down to D and
  warned (you cannot accidentally get an A-looking number from extrapolation);
* **tier propagation** — a composed PK + ke0 + PD simulation inherits the worst
  contributing tier.

There is deliberately **no inverse control**: nothing computes the infusion
required to reach a target concentration or BIS (spec §10).
"""
from __future__ import annotations

import ast
import operator
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .load import Dataset, load
from .models import Model, worst_tier
from .reference import (
    Dosing,
    MicroParams,
    Trajectory,
    bmi as _bmi,
    greco_response_surface,
    sigmoid_emax,
    simulate as _simulate_ref,
)
from .export.registry import INTERACTION_KERNELS, KERNELS, parse_amount, parse_rate

Schedule = Sequence[Tuple[str, float, str]]


# --------------------------------------------------------------------------- #
# Safe predicate evaluation for known_failure_modes
# --------------------------------------------------------------------------- #
_CMP_OPS = {
    ast.Lt: operator.lt, ast.LtE: operator.le, ast.Gt: operator.gt,
    ast.GtE: operator.ge, ast.Eq: operator.eq, ast.NotEq: operator.ne,
}
_BIN_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
            ast.Div: operator.truediv, ast.Pow: operator.pow}
_BOOL_OPS = {ast.And: all, ast.Or: any}


def _eval_predicate(expr: str, env: Dict[str, float]) -> bool:
    """Evaluate a simple comparison predicate (e.g. 'bmi > 42') over covariates."""
    tree = ast.parse(expr, mode="eval").body

    def ev(node):
        if isinstance(node, ast.BoolOp):
            vals = [ev(v) for v in node.values]
            return _BOOL_OPS[type(node.op)](vals)
        if isinstance(node, ast.Compare):
            left = ev(node.left)
            ok = True
            for op, comp in zip(node.ops, node.comparators):
                right = ev(comp)
                ok = ok and _CMP_OPS[type(op)](left, right)
                left = right
            return ok
        if isinstance(node, ast.BinOp):
            return _BIN_OPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -ev(node.operand)
        if isinstance(node, ast.Name):
            if node.id not in env:
                raise ValueError(f"predicate references unknown covariate '{node.id}'")
            return env[node.id]
        if isinstance(node, ast.Constant):
            return node.value
        raise ValueError(f"unsupported predicate construct: {ast.dump(node)}")

    return bool(ev(tree))


def _covariate_env(patient: Dict[str, Any]) -> Dict[str, float]:
    env = {k: v for k, v in patient.items() if isinstance(v, (int, float))}
    if "weight" in env and "height" in env and env.get("height"):
        env.setdefault("bmi", _bmi(env["weight"], env["height"]))
    return env


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #
@dataclass
class SimulationResult:
    model_id: str
    t: np.ndarray
    cp: np.ndarray
    ce: np.ndarray
    tier: str
    warnings: List[str] = field(default_factory=list)
    excluded: bool = False
    effect: Optional[np.ndarray] = None
    effect_label: Optional[str] = None
    pd_model_id: Optional[str] = None
    params: Optional[MicroParams] = None
    patient: Dict[str, Any] = field(default_factory=dict)

    # convenience scalars
    @property
    def cp_peak(self) -> float:
        return float(np.max(self.cp))

    @property
    def ce_peak(self) -> float:
        return float(np.max(self.ce))


# --------------------------------------------------------------------------- #
# Schedule -> Dosing
# --------------------------------------------------------------------------- #
def build_dosing(schedule: Schedule, weight: float) -> Dosing:
    boluses: List[Tuple[float, float]] = []
    infusions: List[Tuple[float, float]] = []
    for kind, t0, spec in schedule:
        if kind == "bolus":
            boluses.append((float(t0), parse_amount(spec, weight)))
        elif kind == "infusion":
            infusions.append((float(t0), parse_rate(spec, weight)))
        else:
            raise ValueError(f"unknown dosing kind {kind!r} (expected 'bolus' or 'infusion')")
    return Dosing(boluses=tuple(boluses), infusions=tuple(infusions))


# --------------------------------------------------------------------------- #
# Envelope + failure-mode evaluation -> (tier_floor, warnings)
# --------------------------------------------------------------------------- #
def evaluate_safety(model: Model, patient: Dict[str, Any]) -> Tuple[str, List[str], bool]:
    """Return (tier_floor, warnings, envelope_violated)."""
    warnings: List[str] = []
    tier_floor = model.tier
    envelope_violated = False

    viols = model.applicability_envelope.check(patient)
    if viols:
        envelope_violated = True
        tier_floor = "D"
        for v in viols:
            warnings.append(f"ENVELOPE: {v} -> tiered down to D")

    env = _covariate_env(patient)
    for fm in model.known_failure_modes:
        triggered = False
        if fm.predicate:
            try:
                triggered = _eval_predicate(fm.predicate, env)
            except ValueError:
                triggered = False
        if triggered:
            if fm.action == "tier_down_to_D":
                tier_floor = worst_tier([tier_floor, "D"])
                warnings.append(f"FAILURE MODE [{fm.condition}]: {fm.behavior} -> tiered down to D")
            elif fm.action == "exclude":
                envelope_violated = True
                tier_floor = "D"
                warnings.append(f"FAILURE MODE [{fm.condition}]: {fm.behavior} -> excluded")
            else:  # warn
                warnings.append(f"WARNING [{fm.condition}]: {fm.behavior}")
    return tier_floor, warnings, envelope_violated


# --------------------------------------------------------------------------- #
# simulate
# --------------------------------------------------------------------------- #
def simulate(
    ds: Dataset,
    model_id: str,
    *,
    patient: Dict[str, Any],
    schedule: Schedule,
    t: np.ndarray,
    pd_model: Optional[str] = None,
) -> SimulationResult:
    """Forward-simulate one PK (optionally + PD) model for one virtual patient."""
    model = ds[model_id]
    if model.purpose != "pk":
        raise ValueError(
            f"simulate() expects a PK model; {model_id} has purpose '{model.purpose}'. "
            "Pass a PK model and attach a PD model via pd_model=..."
        )
    if not model.kernel_implemented:
        raise NotImplementedError(
            f"{model_id}: reference kernel is not implemented (kernel pending verified "
            "transcription). Hypnos refuses to simulate rather than risk a mis-transcribed "
            "covariate equation. See the record's notes field."
        )
    kernel = KERNELS[model.kernel_function]
    params = kernel(patient)

    t = np.asarray(t, dtype=float)
    weight = float(patient.get("weight", 70.0))
    dosing = build_dosing(schedule, weight)
    traj: Trajectory = _simulate_ref(params, dosing, t)

    tier_floor, warnings, excluded = evaluate_safety(model, patient)
    tier = worst_tier([model.tier, tier_floor])

    result = SimulationResult(
        model_id=model_id, t=t, cp=traj.cp, ce=traj.ce, tier=tier,
        warnings=warnings, excluded=excluded, params=params, patient=dict(patient),
    )

    if pd_model is not None:
        pdm = ds[pd_model]
        if not pdm.kernel_implemented:
            raise NotImplementedError(f"{pd_model}: PD kernel not implemented")
        pp = {p.symbol: p.central for p in pdm.parameters}
        effect = sigmoid_emax(traj.ce, pp["E0"], pp["Emax"], pp["Ce50"], pp["gamma"])
        result.effect = effect
        result.effect_label = pdm.label
        result.pd_model_id = pd_model
        # composed simulation inherits the worst tier among PK, PD, and envelope floor
        result.tier = worst_tier([result.tier, pdm.tier])
        pd_floor, pd_warn, _ = evaluate_safety(pdm, patient)
        result.tier = worst_tier([result.tier, pd_floor])
        result.warnings.extend(w for w in pd_warn if w not in result.warnings)

    return result


# --------------------------------------------------------------------------- #
# compare — the model-divergence headline feature
# --------------------------------------------------------------------------- #
@dataclass
class Comparison:
    drug: str
    purpose: str
    t: np.ndarray
    included: List[SimulationResult] = field(default_factory=list)
    excluded: List[Dict[str, Any]] = field(default_factory=list)     # envelope-violating
    unavailable: List[Dict[str, Any]] = field(default_factory=list)  # kernel pending
    divergence: Dict[str, Any] = field(default_factory=dict)


def _divergence(results: List[SimulationResult], key: str) -> Dict[str, float]:
    """Pointwise spread across models for cp or ce. Reports peak absolute and relative spread."""
    if len(results) < 2:
        return {}
    stack = np.vstack([getattr(r, key) for r in results])  # (n_models, n_t)
    spread = stack.max(axis=0) - stack.min(axis=0)
    mean = stack.mean(axis=0)
    rel = np.divide(spread, mean, out=np.zeros_like(spread), where=mean > 1e-9)
    return {
        "max_abs": float(spread.max()),
        "max_rel": float(rel.max()),
        "mean_rel": float(rel.mean()),
    }


def compare(
    ds: Dataset,
    *,
    drug: str,
    patient: Dict[str, Any],
    schedule: Schedule,
    t: np.ndarray,
    purpose: str = "pk",
    pd_model: Optional[str] = None,
) -> Comparison:
    """Overlay every eligible model for a drug; grey out envelope-violators; quantify divergence."""
    from .filter import select

    t = np.asarray(t, dtype=float)
    cmp = Comparison(drug=drug, purpose=purpose, t=t)
    for m in select(ds, drug=drug, purpose=purpose):
        if not m.kernel_implemented:
            cmp.unavailable.append({"model_id": m.id, "tier": m.tier,
                                    "reason": "reference kernel pending verified transcription"})
            continue
        res = simulate(ds, m.id, patient=patient, schedule=schedule, t=t, pd_model=pd_model)
        if res.excluded:
            cmp.excluded.append({"model_id": m.id, "tier": res.tier, "reasons": res.warnings,
                                 "result": res})
        else:
            cmp.included.append(res)

    cmp.divergence = {
        "cp": _divergence(cmp.included, "cp"),
        "ce": _divergence(cmp.included, "ce"),
    }
    return cmp


# --------------------------------------------------------------------------- #
# simulate_interaction — two-drug response surface (hypnotic synergy)
# --------------------------------------------------------------------------- #
@dataclass
class InteractionResult:
    surface_id: str
    pk_a_id: str
    pk_b_id: str
    t: np.ndarray
    ce_a: np.ndarray
    ce_b: np.ndarray
    effect: np.ndarray
    tier: str
    effect_label: str
    warnings: List[str] = field(default_factory=list)

    @property
    def effect_min(self) -> float:
        return float(np.min(self.effect))


def simulate_interaction(
    ds: Dataset,
    surface_id: str,
    *,
    pk_a: str,
    pk_b: str,
    patient: Dict[str, Any],
    schedule_a: Schedule,
    schedule_b: Schedule,
    t: np.ndarray,
) -> InteractionResult:
    """Forward-simulate a two-drug response surface.

    Drug A is the hypnotic (e.g. propofol), drug B the opioid (e.g.
    remifentanil). Each PK model is simulated independently to its effect-site
    concentration, then the interaction surface maps (Ce_a, Ce_b) -> effect.
    The composed result inherits the **worst** tier among PK-A, PK-B, the
    surface, and any envelope floor (worst input wins).
    """
    surf = ds[surface_id]
    if surf.purpose != "interaction":
        raise ValueError(f"{surface_id} has purpose '{surf.purpose}', expected 'interaction'")
    if not surf.kernel_implemented or surf.kernel_function not in INTERACTION_KERNELS:
        raise NotImplementedError(f"{surface_id}: interaction kernel not implemented")

    t = np.asarray(t, dtype=float)
    res_a = simulate(ds, pk_a, patient=patient, schedule=schedule_a, t=t)
    res_b = simulate(ds, pk_b, patient=patient, schedule=schedule_b, t=t)

    sp = {p.symbol: p.central for p in surf.parameters}
    effect = greco_response_surface(
        res_a.ce, res_b.ce,
        E0=sp["E0"], Emax=sp["Emax"],
        Ce50_a=sp["Ce50_prop"], Ce50_b=sp["Ce50_remi"],
        alpha=sp["alpha"], gamma=sp["gamma"],
    )

    surf_floor, surf_warn, _ = evaluate_safety(surf, patient)
    tier = worst_tier([res_a.tier, res_b.tier, surf.tier, surf_floor])
    warnings: List[str] = []
    warnings += [f"[{pk_a}] {w}" for w in res_a.warnings]
    warnings += [f"[{pk_b}] {w}" for w in res_b.warnings]
    warnings += surf_warn

    return InteractionResult(
        surface_id=surface_id, pk_a_id=pk_a, pk_b_id=pk_b, t=t,
        ce_a=res_a.ce, ce_b=res_b.ce, effect=effect, tier=tier,
        effect_label=surf.label, warnings=warnings,
    )
