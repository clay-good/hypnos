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
import math
import operator
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .load import Dataset, load
from .models import Model, concentration_factor, worst_tier
from .reference import (
    Dosing,
    MicroParams,
    Trajectory,
    apply_residual,
    bmi as _bmi,
    greco_response_surface,
    residual_std,
    sample_individual,
    sigmoid_emax,
    sigmoid_emax_twoslope,
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
    concentration_unit: str = "ug/mL"   # conventional display unit for this drug

    # --- v0.2 prediction band (None unless bands requested AND the model carries
    # published between-subject variability — the never-synthesize rule) --------
    cp_quantiles: Optional[Dict[int, np.ndarray]] = None
    ce_quantiles: Optional[Dict[int, np.ndarray]] = None
    band_tier: Optional[str] = None
    band_percentile: Optional[Tuple[int, int]] = None
    band_includes_residual: bool = False
    # per-time variance components, for the compare() variance decomposition (§7.2)
    cp_bsv_var: Optional[np.ndarray] = None
    ce_bsv_var: Optional[np.ndarray] = None
    cp_resid_var: Optional[np.ndarray] = None
    ce_resid_var: Optional[np.ndarray] = None

    # convenience scalars — cp/ce arrays are always internal ug/mL (== mg/L)
    @property
    def cp_peak(self) -> float:
        return float(np.max(self.cp))

    @property
    def ce_peak(self) -> float:
        return float(np.max(self.ce))

    @property
    def conc_factor(self) -> float:
        return concentration_factor(self.concentration_unit)

    @property
    def cp_peak_display(self) -> float:
        """Peak plasma concentration in the drug's conventional unit (e.g. ng/mL for opioids)."""
        return self.cp_peak * self.conc_factor

    @property
    def ce_peak_display(self) -> float:
        return self.ce_peak * self.conc_factor


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
def _classify_age_extrapolation(model: Model, patient: Dict[str, Any]) -> Optional[str]:
    """Name an age extrapolation (pediatric/geriatric) when a patient falls outside
    the model's derivation age range. Returns a labeled Tier-D message or None."""
    age = patient.get("age")
    if age is None:
        return None
    rng = model.applicability_envelope.age_years
    pops = [p.lower() for p in model.applicability_envelope.populations]
    is_pediatric_model = any(p in ("child", "pediatric", "neonate", "infant") for p in pops)
    if rng.min is not None and age < rng.min and age < 18:
        return (f"PEDIATRIC EXTRAPOLATION: age {age:g} y is below the model's "
                f"derivation range (>= {rng.min:g} y); an adult model used in a child "
                "is not predictive -> Tier D")
    if rng.max is not None and age > rng.max and is_pediatric_model:
        return (f"EXTRAPOLATION: age {age:g} y exceeds this pediatric model's range "
                f"(<= {rng.max:g} y); a pediatric model used in an adult is not predictive -> Tier D")
    if rng.max is not None and age > rng.max and age >= 65:
        return (f"GERIATRIC EXTRAPOLATION: age {age:g} y exceeds the model's derivation "
                f"range (<= {rng.max:g} y) -> Tier D")
    return None


def _apply_pd(pdm: Model, ce: np.ndarray, patient: Dict[str, Any]) -> np.ndarray:
    """Map effect-site concentration to effect via the PD model's kernel.

    Dispatches on the PD record's ``kernel.function``: the Eleveld BIS model uses
    a two-slope sigmoid with an age-corrected Ce50; the others use a single-slope
    sigmoid E_max with parameters {E0, Emax, Ce50, gamma} read from the record.
    """
    pp = {p.symbol: p.central for p in pdm.parameters}
    if pdm.kernel_function == "eleveld_bis_twoslope":
        age = float(patient.get("age", 35.0))
        ce50 = pp["Ce50"] * float(np.exp(pp.get("Ce50_age_coeff", 0.0) * (age - 35.0)))
        return sigmoid_emax_twoslope(ce, pp["E0"], pp["Emax"], ce50, pp["gamma_low"], pp["gamma_high"])
    return sigmoid_emax(ce, pp["E0"], pp["Emax"], pp["Ce50"], pp["gamma"])


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
        # Explicit, categorized extrapolation labeling (spec §11, Phase C):
        # an adult model used in a child, or a pediatric model used in an adult,
        # is not merely "out of envelope" — name the extrapolation.
        extrap = _classify_age_extrapolation(model, patient)
        if extrap:
            warnings.append(extrap)

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
# Prediction bands — seeded Monte-Carlo over the curated random effects (§6)
# --------------------------------------------------------------------------- #
def _residual_primitives(re_) -> Tuple[str, Dict[str, float]]:
    """Translate a curated ResidualError into the reference helper's primitives."""
    kw: Dict[str, float] = {}
    if re_.proportional:
        v = re_.proportional.get("variance")
        if v is None and re_.proportional.get("cv_percent") is not None:
            v = (re_.proportional["cv_percent"] / 100.0) ** 2
        kw["prop_var"] = v
    if re_.additive:
        sd = re_.additive.get("sd")
        if sd is None and re_.additive.get("variance") is not None:
            sd = math.sqrt(re_.additive["variance"])
        kw["add_sd"] = sd
    if re_.log:
        sd = re_.log.get("sd")
        if sd is None and re_.log.get("variance") is not None:
            sd = math.sqrt(re_.log["variance"])
        kw["log_sd"] = sd
    return re_.model, kw


def _attach_bands(
    result: SimulationResult,
    model: Model,
    params: MicroParams,
    dosing: Dosing,
    t: np.ndarray,
    *,
    percentile: Tuple[int, int],
    samples: int,
    seed: int,
    residual: bool,
) -> None:
    """Populate the prediction-band fields on ``result`` in place.

    Honors the never-synthesize rule (spec §5): a model with no published BSV
    draws no band and instead records why. Bands are BSV-only by default (where is
    the individual's *true* curve?); ``residual`` adds Σ for observation-level bands
    (where will a *measured sample* land?).
    """
    lo, hi = int(percentile[0]), int(percentile[1])

    if not model.has_published_variability:
        result.warnings.append(
            "BAND: no published between-subject variability for this model — no band "
            "drawn (never-synthesize rule; the median line keeps its tier)"
        )
        return

    omegas = model.bsv_omegas()
    if not omegas:
        result.warnings.append(
            "BAND: variability_status declares random effects but no parameter carries "
            "a usable omega2 — no band drawn"
        )
        return

    # which structural parameters got NO BSV (band is a lower bound on true spread)
    fixed = [p.symbol for p in model.parameters
             if p.central is not None and p.symbol not in omegas
             and p.symbol in ("V1", "V2", "V3", "Cl1", "Cl2", "Cl3", "ke0")]

    rng = np.random.default_rng(seed)
    n = len(t)
    cp_draws = np.empty((samples, n))
    ce_draws = np.empty((samples, n))
    for i in range(samples):
        mp = sample_individual(params, omegas, rng)
        traj = _simulate_ref(mp, dosing, t)
        cp_draws[i] = traj.cp
        ce_draws[i] = traj.ce

    result.cp_bsv_var = cp_draws.var(axis=0)
    result.ce_bsv_var = ce_draws.var(axis=0)

    # residual variance from Σ evaluated at the deterministic (median) curve
    if model.residual_error is not None:
        rmodel, rkw = _residual_primitives(model.residual_error)
        result.cp_resid_var = residual_std(result.cp, rmodel, **rkw) ** 2
        result.ce_resid_var = residual_std(result.ce, rmodel, **rkw) ** 2
        if residual:
            cp_draws = apply_residual(cp_draws, rmodel, rng, **rkw)
            ce_draws = apply_residual(ce_draws, rmodel, rng, **rkw)
    else:
        result.cp_resid_var = np.zeros(n)
        result.ce_resid_var = np.zeros(n)

    qs = [lo, 50, hi]
    cp_q = np.percentile(cp_draws, qs, axis=0)
    ce_q = np.percentile(ce_draws, qs, axis=0)
    result.cp_quantiles = {q: cp_q[k] for k, q in enumerate(qs)}
    result.ce_quantiles = {q: ce_q[k] for k, q in enumerate(qs)}
    result.band_percentile = (lo, hi)
    result.band_includes_residual = bool(residual)

    # band tier = worst of the structural/variability tiers and the envelope floor
    bt = model.band_tier or model.tier
    result.band_tier = worst_tier([bt, result.tier])

    if model.variability_status == "partial" and fixed:
        result.warnings.append(
            "BAND: partial variability — parameters without published BSV held fixed "
            f"({', '.join(fixed)}); the band is a LOWER BOUND on true between-subject spread"
        )
    if model.variability_status == "diagonal":
        result.warnings.append(
            "BAND: diagonal Omega — off-diagonal correlations not published; η's drawn "
            "independently (recorded caveat, spec §5)"
        )
    if residual:
        result.warnings.append(
            "BAND: includes residual error Σ — this is an OBSERVATION-level band "
            "(where a measured sample lands), not the individual's true curve"
        )


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
    bands: bool = False,
    percentile: Tuple[int, int] = (5, 95),
    samples: int = 2000,
    seed: Optional[int] = None,
    residual: bool = False,
) -> SimulationResult:
    """Forward-simulate one PK (optionally + PD) model for one virtual patient.

    With ``bands=True`` (and a mandatory integer ``seed``), also draws a seeded
    Monte-Carlo prediction band from the model's curated between-subject variability
    (v0.2 §6). A model that publishes no BSV draws no band — the never-synthesize
    rule (§5): a missing band is honest; a borrowed one is a lie with error bars.
    """
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

    drug_meta = ds.drug(model.drug_name) or {}
    result = SimulationResult(
        model_id=model_id, t=t, cp=traj.cp, ce=traj.ce, tier=tier,
        warnings=warnings, excluded=excluded, params=params, patient=dict(patient),
        concentration_unit=drug_meta.get("concentration_unit", "ug/mL"),
    )

    if pd_model is not None:
        pdm = ds[pd_model]
        if not pdm.kernel_implemented:
            raise NotImplementedError(f"{pd_model}: PD kernel not implemented")
        effect = _apply_pd(pdm, traj.ce, patient)
        result.effect = effect
        result.effect_label = pdm.label
        result.pd_model_id = pd_model
        # composed simulation inherits the worst tier among PK, PD, and envelope floor
        result.tier = worst_tier([result.tier, pdm.tier])
        pd_floor, pd_warn, _ = evaluate_safety(pdm, patient)
        result.tier = worst_tier([result.tier, pd_floor])
        result.warnings.extend(w for w in pd_warn if w not in result.warnings)

    if bands:
        if seed is None:
            raise ValueError(
                "simulate(bands=True) requires an explicit integer seed — every "
                "band-producing call is seeded so quantiles are byte-reproducible (spec §6)."
            )
        _attach_bands(result, model, params, dosing, t, percentile=percentile,
                      samples=samples, seed=seed, residual=residual)

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
    concentration_unit: str = "ug/mL"   # conventional display unit for this drug
    bands: bool = False
    # band-eligible models contribute to the separation index & variance decomposition;
    # models with variability_status == "none" are named here, never silently dropped (§7.2)
    excluded_from_bands: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def conc_factor(self) -> float:
        return concentration_factor(self.concentration_unit)


def _divergence(results: List[SimulationResult], key: str) -> Dict[str, Any]:
    """Pointwise spread across models for cp or ce. Reports peak absolute and
    relative spread, plus the **driver pair** — the two models furthest apart at the
    instant of peak disagreement. The pooled spread says *how much* the models
    disagree; the driver says *which* model is the outlier (e.g. Schnider vs the
    rest), the actionable half of model-selection risk."""
    if len(results) < 2:
        return {}
    stack = np.vstack([getattr(r, key) for r in results])  # (n_models, n_t)
    spread = stack.max(axis=0) - stack.min(axis=0)
    mean = stack.mean(axis=0)
    rel = np.divide(spread, mean, out=np.zeros_like(spread), where=mean > 1e-9)
    t_star = int(np.argmax(spread))               # the instant of widest disagreement
    col = stack[:, t_star]
    hi, lo = int(np.argmax(col)), int(np.argmin(col))
    return {
        "max_abs": float(spread.max()),
        "max_rel": float(rel.max()),
        "mean_rel": float(rel.mean()),
        "driver": {
            "high": results[hi].model_id,
            "low": results[lo].model_id,
            "gap": float(col[hi] - col[lo]),       # == max_abs, in internal ug/mL
        },
    }


def _band_divergence(results: List[SimulationResult], key: str) -> Dict[str, Any]:
    """Uncertainty-aware divergence (v0.2 §7) for the band-eligible subset.

    Adds two readouts to the point-estimate divergence:

    * **separation index** — at the instant of peak median spread t*, are the two
      driver models' percentile bands disjoint? ``separation > 0`` => a genuine
      structural disagreement neither model's stated BSV explains away.
    * **variance decomposition** — what share of the total predictive variance at t*
      is between-model (structural), within-model (BSV), and residual (Σ)?
    """
    eligible = [r for r in results if getattr(r, f"{key}_quantiles") is not None]
    if not eligible:
        return {}
    lo, hi = eligible[0].band_percentile
    med = np.vstack([r.__dict__[f"{key}_quantiles"][50] for r in eligible])  # (n,t)
    qlo = np.vstack([r.__dict__[f"{key}_quantiles"][lo] for r in eligible])
    qhi = np.vstack([r.__dict__[f"{key}_quantiles"][hi] for r in eligible])

    out: Dict[str, Any] = {}

    # --- variance decomposition at t* (peak median spread) -----------------
    spread = med.max(axis=0) - med.min(axis=0)
    t_star = int(np.argmax(spread)) if med.shape[0] > 1 else int(np.argmax(med.mean(axis=0)))
    var_structural = float(med[:, t_star].var()) if med.shape[0] > 1 else 0.0
    bsv_var = np.mean([r.__dict__[f"{key}_bsv_var"][t_star] for r in eligible])
    resid_var = np.mean([r.__dict__[f"{key}_resid_var"][t_star] for r in eligible])
    total = var_structural + bsv_var + resid_var
    if total > 0:
        out["variance_share"] = {
            "structural": round(float(var_structural / total), 4),
            "bsv": round(float(bsv_var / total), 4),
            "residual": round(float(resid_var / total), 4),
            "t_star_min": float(eligible[0].t[t_star]),
        }

    # --- separation index (needs >= 2 band-eligible models) ----------------
    if len(eligible) >= 2:
        # disjointness across the whole grid, ranking high/low by median per instant
        hi_idx = med.argmax(axis=0)
        lo_idx = med.argmin(axis=0)
        cols = np.arange(med.shape[1])
        gap = qlo[hi_idx, cols] - qhi[lo_idx, cols]    # >0 => bands disjoint
        frac_disjoint = float(np.mean(gap > 0))

        hh, ll = int(med[:, t_star].argmax()), int(med[:, t_star].argmin())
        gap_star = float(qlo[hh, t_star] - qhi[ll, t_star])
        width_hi = qhi[hh, t_star] - qlo[hh, t_star]
        width_lo = qhi[ll, t_star] - qlo[ll, t_star]
        pooled = 0.5 * (width_hi + width_lo)
        out["separation"] = {
            "value": round(float(gap_star / pooled), 4) if pooled > 1e-12 else None,
            "bands_disjoint_at_tstar": bool(gap_star > 0),
            "fraction_trajectory_disjoint": round(frac_disjoint, 4),
            "percentile": [lo, hi],
            "driver_high": eligible[hh].model_id,
            "driver_low": eligible[ll].model_id,
            "band_tier": worst_tier([eligible[hh].band_tier, eligible[ll].band_tier]),
        }
    return out


def compare(
    ds: Dataset,
    *,
    drug: str,
    patient: Dict[str, Any],
    schedule: Schedule,
    t: np.ndarray,
    purpose: str = "pk",
    pd_model: Optional[str] = None,
    bands: bool = False,
    percentile: Tuple[int, int] = (5, 95),
    samples: int = 2000,
    seed: Optional[int] = None,
    residual: bool = False,
) -> Comparison:
    """Overlay every eligible model for a drug; grey out envelope-violators; quantify divergence.

    With ``bands=True`` (and a ``seed``), each band-eligible model also carries a
    seeded prediction band, and the divergence view answers the v0.2 headline
    question: *are the models distinguishable beyond their own stated variability?*
    (separation index) and *what dominates the uncertainty here?* (variance
    decomposition). Models that publish no BSV are named in ``excluded_from_bands``,
    never silently dropped.
    """
    from .filter import select

    if bands and seed is None:
        raise ValueError("compare(bands=True) requires an explicit integer seed (spec §6).")

    t = np.asarray(t, dtype=float)
    cmp = Comparison(drug=drug, purpose=purpose, t=t, bands=bands,
                     concentration_unit=(ds.drug(drug) or {}).get("concentration_unit", "ug/mL"))
    for m in select(ds, drug=drug, purpose=purpose):
        if not m.kernel_implemented:
            cmp.unavailable.append({"model_id": m.id, "tier": m.tier,
                                    "reason": "reference kernel pending verified transcription"})
            continue
        res = simulate(ds, m.id, patient=patient, schedule=schedule, t=t, pd_model=pd_model,
                       bands=bands, percentile=percentile, samples=samples, seed=seed,
                       residual=residual)
        if res.excluded:
            cmp.excluded.append({"model_id": m.id, "tier": res.tier, "reasons": res.warnings,
                                 "result": res})
        else:
            cmp.included.append(res)
            if bands and not m.has_published_variability:
                cmp.excluded_from_bands.append(
                    {"model_id": m.id, "tier": res.tier,
                     "reason": "variability_status: none — no published between-subject variability"})

    # Plasma divergence spans every included PK model; effect-site divergence is
    # only meaningful among models that actually carry an effect compartment
    # (a ke0 link). PK-only models (e.g. Kim remifentanil, Paedfusor) have ce==0
    # and must not pollute the ce spread.
    ce_results = [r for r in cmp.included if r.params is not None and r.params.ke0 > 0]
    cmp.divergence = {
        "cp": _divergence(cmp.included, "cp"),
        "ce": _divergence(ce_results, "ce"),
    }
    if bands:
        for k, rs in (("cp", cmp.included), ("ce", ce_results)):
            band_part = _band_divergence(rs, k)
            if band_part:
                cmp.divergence.setdefault(k, {}).update(band_part)
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
