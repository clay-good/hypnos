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

from .load import Dataset
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
    sample_parameter_vector,
    sigmoid_emax,
    sigmoid_emax_twoslope,
    simulate as _simulate_ref,
)
from .developmental import apply_developmental, linear_per_kg_scale, maturation_value
from .pharmacogenomics import genotype_triggers
from .export.registry import INTERACTION_KERNELS, KERNELS, parse_amount, parse_rate
from .covariates import (
    covariate_layer_tier as _cov_layer_tier,
    distribution_keys as _cov_dist_keys,
    evaluate as _cov_eval,
    has_covariate_distribution as _has_cov_dist,
    point_patient as _point_patient,
    sample_covariate_vector as _sample_cov,
)

# Body-size equations a covariate slot can be substituted by (the C1 override set);
# shared with covariate_divergence for the equation-choice variance component (§7.3).
_BAND_KINDS = ("prediction", "covariate")


def _normalize_bands(bands) -> set:
    """Normalize the ``bands`` argument to a set of kinds.

    Backward-compatible: ``True`` -> {"prediction"} (the v0.2 BSV band), ``False`` -> {},
    a str -> {str}, a sequence -> set(sequence). So ``bands=True`` keeps its v0.2 meaning
    while ``bands=["prediction","covariate"]`` opts into the v0.7 covariate band too.
    """
    if not bands:
        return set()
    if bands is True:
        return {"prediction"}
    if isinstance(bands, str):
        return {bands}
    return set(bands)


# Band kinds Hypnos can draw: v0.2 prediction (BSV), v0.3 confidence (estimation
# uncertainty on θ), v0.7 covariate (covariate-value / equation-choice). Distinct objects.
_KNOWN_BAND_KINDS = ("prediction", "confidence", "covariate")

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
    # PD effect band: the PK between-subject variability propagated through the
    # (deterministic) PD link. Present only when bands AND a pd_model are composed.
    # PD-parameter BSV (Ce50, gamma) is not curated, so this is a LOWER BOUND on the
    # true effect spread — labeled as such in the warnings (spec §14).
    effect_quantiles: Optional[Dict[int, np.ndarray]] = None
    band_tier: Optional[str] = None
    band_percentile: Optional[Tuple[int, int]] = None
    band_includes_residual: bool = False
    # --- v0.3 confidence band (None unless requested AND the model carries published
    # estimation uncertainty on θ — the never-synthesize rule). A NARROW band around the
    # typical curve (how well the data pins the mean down), distinct from the wide BSV band.
    cp_confidence_quantiles: Optional[Dict[int, np.ndarray]] = None
    ce_confidence_quantiles: Optional[Dict[int, np.ndarray]] = None
    confidence_band_tier: Optional[str] = None
    # per-time estimation variance, for the compare() four-way decomposition (v0.3 E2)
    cp_est_var: Optional[np.ndarray] = None
    ce_est_var: Optional[np.ndarray] = None
    # per-time variance components, for the compare() variance decomposition (§7.2)
    cp_bsv_var: Optional[np.ndarray] = None
    ce_bsv_var: Optional[np.ndarray] = None
    cp_resid_var: Optional[np.ndarray] = None
    ce_resid_var: Optional[np.ndarray] = None

    # --- v0.7 C2 covariate band (None unless requested AND a covariate-value
    # distribution is supplied — the never-invent rule) --------------------------
    cp_covariate_band: Optional[Dict[int, np.ndarray]] = None
    ce_covariate_band: Optional[Dict[int, np.ndarray]] = None
    effect_covariate_band: Optional[Dict[int, np.ndarray]] = None
    covariate_band_tier: Optional[str] = None
    covariate_warnings: List[str] = field(default_factory=list)
    # per-time covariate variance components for the §7.3 fifth-component decomposition:
    # value uncertainty (the supplied input distribution) and equation choice (which
    # body-size equation), kept distinct — the two halves of §1.
    cp_cov_value_var: Optional[np.ndarray] = None
    ce_cov_value_var: Optional[np.ndarray] = None
    cp_cov_equation_var: Optional[np.ndarray] = None
    ce_cov_equation_var: Optional[np.ndarray] = None

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


def evaluate_safety(
    model: Model, patient: Dict[str, Any], drug_meta: Optional[Dict[str, Any]] = None
) -> Tuple[str, List[str], bool]:
    """Return (tier_floor, warnings, envelope_violated).

    ``drug_meta`` (the drug record) is optional; when supplied and the patient is
    hypoalbuminemic, a *binding-sensitive* drug surfaces the cited protein-binding /
    free-fraction failure mode (v0.5 §B3) — the free-fraction shift is named, never
    silently modeled.
    """
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

    # Physiological (organ-function) envelope (v0.5 §B): hepatic/renal/cardiac/albumin.
    # An organ-failure patient is greyed + Tier-D for every model with no cited standing;
    # a model with standing (e.g. remifentanil's esterase clearance) carries an explaining
    # note instead. Independent of the demographic check above — organ failure can violate
    # the envelope even when age/weight/BMI are in range.
    albumin_flagged = False
    for f in model.applicability_envelope.organ_check(patient):
        if f.axis == "albumin":
            albumin_flagged = True
        if f.extrapolation:
            envelope_violated = True
            tier_floor = "D"
            warnings.append(f"ORGAN ENVELOPE: {f.message}")
        else:
            warnings.append(f"ORGAN NOTE: {f.message}")

    # Protein-binding / free-fraction failure mode (v0.5 §B3): for a binding-sensitive
    # drug, hypoalbuminemia raises the free (active) fraction, so the TOTAL-concentration
    # prediction under-estimates effect. Surfaced with its citation, never modeled.
    pb = (drug_meta or {}).get("protein_binding") or {}
    if albumin_flagged and pb.get("binding_sensitive"):
        fb = pb.get("fraction_bound")
        bound = f"~{fb * 100:g}% protein-bound" if fb is not None else "highly protein-bound"
        cite = pb.get("citation")
        warnings.append(
            f"BINDING-SENSITIVE: {model.drug_name} is {bound}; in hypoalbuminemia the free "
            "(active) fraction rises, so the total-concentration prediction UNDER-estimates "
            "effect — a documented failure mode, surfaced not modeled"
            + (f" [{cite}]" if cite else "")
        )

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

    # Pharmacogenomic safety flags (v0.9 §6): a DECLARED genotype that matches a curated
    # avoidance/awareness flag is surfaced as a contraindication-style warning with NO
    # numeric effect (a susceptibility is never a dose change — v0.9 §2 Trap 2). Never
    # inferred: only an explicitly declared genotype dimension can trigger.
    for flag in model.pharmacogenomic_safety_flags:
        if genotype_triggers(patient, flag.phenotype_dimension, flag.phenotype_value):
            tag = "AVOID" if flag.kind == "avoidance" else "AWARENESS"
            triggers = (" trigger(s): " + ", ".join(flag.trigger_agents)) if flag.trigger_agents else ""
            cite = f" [{flag.primary_citation}]" if flag.primary_citation else ""
            warnings.append(
                f"PGx SAFETY [{tag}, {flag.gene}]: {flag.action} — {flag.consequence}"
                f"{triggers}. This is a contraindication/awareness flag, NOT a dose change"
                f" (Tier {flag.evidence_tier} on the gene link){cite}")
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
    pd_model: Optional[Model] = None,
    patient: Optional[Dict[str, Any]] = None,
) -> None:
    """Populate the prediction-band fields on ``result`` in place.

    Honors the never-synthesize rule (spec §5): a model with no published BSV
    draws no band and instead records why. Bands are BSV-only by default (where is
    the individual's *true* curve?); ``residual`` adds Σ for observation-level bands
    (where will a *measured sample* land?).

    When ``pd_model`` is supplied, each virtual individual's *true* effect-site curve
    is pushed through the PD link to build an effect band (e.g. a BIS band). PD-parameter
    BSV (Ce50, γ) is not curated, so the effect band reflects only the PK between-subject
    variability and is a LOWER BOUND on true effect spread (spec §14).
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

    # PD effect band — push each individual's TRUE (BSV-only, pre-residual) effect-site
    # curve through the deterministic PD link. The Hill transform is non-linear and
    # monotone, so quantiles are taken on the effect draws directly (not mapped from
    # ce quantiles) and stay correct regardless of the effect's direction (BIS falls
    # as ce rises). _apply_pd is numpy-vectorized, so the whole (samples, n) array
    # transforms at once.
    effect_draws = None
    if pd_model is not None and patient is not None:
        effect_draws = _apply_pd(pd_model, ce_draws, patient)

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
    if effect_draws is not None:
        eff_q = np.percentile(effect_draws, qs, axis=0)
        result.effect_quantiles = {q: eff_q[k] for k, q in enumerate(qs)}
        result.warnings.append(
            "BAND: effect band propagates PK between-subject variability through the "
            "(fixed) PD link; PD-parameter BSV (Ce50, gamma) is not curated, so it is a "
            "LOWER BOUND on true effect spread (spec §14)"
        )
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
# Confidence band — seeded draws over the ESTIMATION uncertainty on θ (v0.3 E1)
# --------------------------------------------------------------------------- #
def _attach_confidence_band(
    result: SimulationResult,
    model: Model,
    params: MicroParams,
    dosing: Dosing,
    t: np.ndarray,
    *,
    percentile: Tuple[int, int],
    samples: int,
    seed: int,
) -> None:
    """Populate the confidence-band fields in place (v0.3 E1).

    The confidence band asks *how well does the fitting data pin down the typical curve?* —
    a reducible uncertainty, sampled from the curated per-θ estimation SE. It is a different
    object from the v0.2 prediction band (irreducible between-subject spread) and is typically
    far narrower. Honors the never-synthesize rule: a model with no published estimation
    uncertainty draws no confidence band and records why (v0.3 §5)."""
    lo, hi = int(percentile[0]), int(percentile[1])
    if not model.has_published_estimation:
        result.warnings.append(
            "CONFIDENCE BAND: no published estimation uncertainty (SE on θ) for this model — no "
            "confidence band drawn (never-synthesize rule; distinct from the prediction band)")
        return
    ses = model.estimation_ses()
    if not ses:
        result.warnings.append(
            "CONFIDENCE BAND: uncertainty_status declares estimation uncertainty but no parameter "
            "carries a usable SE — no band drawn")
        return

    rng = np.random.default_rng(seed)
    n = len(t)
    cp_draws = np.empty((samples, n))
    ce_draws = np.empty((samples, n))
    for i in range(samples):
        mp = sample_parameter_vector(params, ses, rng)
        traj = _simulate_ref(mp, dosing, t)
        cp_draws[i] = traj.cp
        ce_draws[i] = traj.ce
    result.cp_est_var = cp_draws.var(axis=0)
    result.ce_est_var = ce_draws.var(axis=0)
    qs = [lo, 50, hi]
    cp_q = np.percentile(cp_draws, qs, axis=0)
    ce_q = np.percentile(ce_draws, qs, axis=0)
    result.cp_confidence_quantiles = {q: cp_q[k] for k, q in enumerate(qs)}
    result.ce_confidence_quantiles = {q: ce_q[k] for k, q in enumerate(qs)}
    result.confidence_band_tier = worst_tier(
        [model.estimation_band_tier or model.tier, result.tier])
    result.band_percentile = result.band_percentile or (lo, hi)
    result.warnings.append(
        "CONFIDENCE BAND: seeded draws over the per-θ estimation SE — this is how well the data "
        "pin down the TYPICAL curve (reducible with more data), NOT the between-subject spread "
        f"(band-tier {result.confidence_band_tier}); the parameters without a curated SE are held fixed")


# --------------------------------------------------------------------------- #
# Covariate band — seeded draws over the covariate VALUE distribution (v0.7 C2 §7.2)
# --------------------------------------------------------------------------- #
def _equation_choice_variance(
    ds: Dataset, model: Model, kernel, point: Dict[str, Any], dosing: Dosing, t: np.ndarray,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Per-time variance across the model's admissible body-size equations (§7.3).

    The *equation-choice* half of the covariate uncertainty: at a fixed (point) covariate
    vector, how much does the prediction move purely because the body-size descriptor could
    be computed by a different equation? Reuses the C1 override mechanism. ``None`` when the
    model declares no covariate_model (covariate_sensitivity_status: none — no equation axis).
    """
    cm = model.covariate_model
    if cm is None:
        return None, None
    di = cm.derived_inputs[0]
    cand = [eid for eid, rec in ds.covariate_equations.items()
            if rec.get("quantity") in ("lbm", "ffm", "nfm", "ibw")]
    cand = [di.equation] + [c for c in cand if c != di.equation]
    cp_curves, ce_curves = [], []
    for eid in cand:
        val = _cov_eval(eid, point, ds=ds).value
        traj = _simulate_ref(kernel({**point, f"_{di.quantity}_override": val}), dosing, t)
        cp_curves.append(traj.cp)
        ce_curves.append(traj.ce)
    return np.vstack(cp_curves).var(axis=0), np.vstack(ce_curves).var(axis=0)


def _attach_covariate_band(
    result: SimulationResult,
    model: Model,
    kernel,
    dosing: Dosing,
    t: np.ndarray,
    *,
    patient: Dict[str, Any],
    point: Dict[str, Any],
    percentile: Tuple[int, int],
    samples: int,
    seed: int,
    pd_model: Optional[Model],
    ds: Dataset,
) -> None:
    """Populate the covariate-band fields in place (v0.7 §7.2/§7.3).

    Two distinct covariate uncertainties, kept separate (the two halves of §1):

    * **equation choice** — the per-time variance across admissible body-size equations at
      the fixed covariate vector (always computed for a covariate-scaled model);
    * **value uncertainty** — the seeded band from a caller-supplied covariate-value
      distribution (drawn only when one is supplied; the never-invent rule, §5).

    The administered dose schedule is held fixed at the point weight, so the band isolates
    the covariate -> PK effect (consistent with the v0.2 BSV band).
    """
    lo, hi = int(percentile[0]), int(percentile[1])
    n = len(t)

    # equation-choice variance (the "which equation" half) — always, if curated
    cp_eqn, ce_eqn = _equation_choice_variance(ds, model, kernel, point, dosing, t)
    result.cp_cov_equation_var = cp_eqn
    result.ce_cov_equation_var = ce_eqn

    # surface the model's own equation status at this patient (e.g. James inverted)
    if model.covariate_model is not None:
        di = model.covariate_model.derived_inputs[0]
        ev = _cov_eval(di.equation, point, ds=ds)
        result.covariate_warnings.extend(ev.warnings)

    # covariate-value band — only when a distribution is supplied (never invent)
    if not _has_cov_dist(patient):
        result.cp_cov_value_var = np.zeros(n)
        result.ce_cov_value_var = np.zeros(n)
        if model.covariate_model is not None:
            ltier = _cov_layer_tier(model, point, ds)
            result.covariate_band_tier = worst_tier([model.tier, ltier]) if ltier else model.tier
        result.warnings.append(
            "COVARIATE BAND: no covariate-value distribution supplied — equation-choice "
            "variance reported, but no value band drawn (never-invent rule, §5)"
        )
        return

    varying = _cov_dist_keys(patient)
    rng = np.random.default_rng(seed)
    cp_draws = np.empty((samples, n))
    ce_draws = np.empty((samples, n))
    for i in range(samples):
        pv = _sample_cov(patient, rng)
        traj = _simulate_ref(kernel(pv), dosing, t)
        cp_draws[i] = traj.cp
        ce_draws[i] = traj.ce
    result.cp_cov_value_var = cp_draws.var(axis=0)
    result.ce_cov_value_var = ce_draws.var(axis=0)

    qs = [lo, 50, hi]
    result.cp_covariate_band = {q: np.percentile(cp_draws, q, axis=0) for q in qs}
    result.ce_covariate_band = {q: np.percentile(ce_draws, q, axis=0) for q in qs}
    if pd_model is not None:
        eff_draws = _apply_pd(pd_model, ce_draws, point)
        result.effect_covariate_band = {q: np.percentile(eff_draws, q, axis=0) for q in qs}

    ltier = _cov_layer_tier(model, point, ds)
    result.covariate_band_tier = worst_tier([model.tier, ltier]) if ltier else model.tier
    result.band_percentile = result.band_percentile or (lo, hi)
    if len(varying) > 1:
        result.warnings.append(
            "COVARIATE BAND: >1 covariate varies; perturbed INDEPENDENTLY (no covariance "
            "assumed — recorded caveat, spec §14)"
        )
    result.warnings.append(
        f"COVARIATE BAND: propagates supplied covariate-value uncertainty ({', '.join(varying)}); "
        "the administered dose schedule is held fixed (covariate -> PK effect only)"
    )


# --------------------------------------------------------------------------- #
# Developmental extrapolation + pharmacogenomic modifiers (opt-in transforms)
# --------------------------------------------------------------------------- #
def _reference_adult_patient(sex: str, ref_wt: float) -> Dict[str, Any]:
    """A nominal reference adult at the allometric reference weight (v0.8 §5).

    Allometry scales the *reference (adult) disposition parameters*, so the adult model's
    typical values are taken at this reference individual and then size/maturation-scaled
    — never the adult covariate kernel run at neonatal covariates (which would break, e.g.
    Schnider's LBM at 3.4 kg)."""
    return {"weight": ref_wt, "height": 170.0, "age": 40.0, "sex": sex}


def _apply_developmental(model: Model, kernel, point: Dict[str, Any]):
    """Build the developmentally-extrapolated MicroParams + Tier-D warnings (v0.8 §5).

    Returns ``(MicroParams, warnings)``. Raises if the model carries no curated
    developmental block (never-invent: an extrapolation is applied only when curated)."""
    dev = model.developmental_model
    if dev is None:
        raise ValueError(
            f"{model.id}: developmental=True but no developmental_model is curated — "
            "Hypnos will not invent an allometric/maturation extrapolation (v0.8 §5).")
    ref_wt = dev.size.reference_weight_kg if dev.size else 70.0
    ref_params = kernel(_reference_adult_patient(str(point.get("sex", "M")), ref_wt))
    vc = ref_params.as_volumes_clearances()
    scaled = apply_developmental(vc, dev, point)
    mp = MicroParams.from_volumes_clearances(
        V1=scaled["V1"], Cl1=scaled["Cl1"], V2=scaled.get("V2", 0.0), Cl2=scaled.get("Cl2", 0.0),
        V3=scaled.get("V3", 0.0), Cl3=scaled.get("Cl3", 0.0), ke0=scaled.get("ke0", 0.0))

    cites = [c for c in [dev.primary_citation,
                         dev.size.primary_citation if dev.size else None,
                         dev.maturation.primary_citation if dev.maturation else None] if c]
    warns = [
        f"DEVELOPMENTAL EXTRAPOLATION ({dev.extrapolation_basis}): adult model carried below "
        f"its derivation age by allometry"
        + ("+maturation" if dev.has_maturation else "")
        + "; NOT validated in this patient -> Tier D"
        + (f" [{', '.join(cites)}]" if cites else "")]
    mf = maturation_value(dev, point)
    if mf is not None:
        warns.append(f"DEVELOPMENTAL: maturation factor MF(PMA={point.get('pma_weeks')}) = {mf:.3f} "
                     f"(TM50={dev.maturation.tm50_weeks} wk, Hill={dev.maturation.hill}) applied to "
                     f"{dev.maturation.affected_parameter}")
    if dev.extrapolation_basis == "allometry_only":
        warns.append(
            "DEVELOPMENTAL CAVEAT (allometry_only): maturation is UN-MODELED — neonatal "
            "clearance is OVER-stated and exposure UNDER-predicted, exactly where the margin "
            "is smallest (v0.8 §5). A missing maturation block is a true gap, not a mature patient.")
    if dev.caveat:
        warns.append(f"DEVELOPMENTAL: {dev.caveat}")
    return mp, warns


def _apply_pgx_modifiers(model: Model, params: "MicroParams", point: Dict[str, Any]):
    """Apply DECLARED, opt-in kinetic pharmacogenomic modifiers (v0.9 G1).

    Returns ``(MicroParams, warnings)``. A modifier fires only when (a) its phenotype is
    explicitly declared and triggered and (b) the model's drug is in the modifier's
    ``substrate_scope`` (the substrate guardrail — a genetic effect never leaks to a
    non-substrate drug, v0.9 Trap 4). Scaling the affected clearance prolongs the
    hydrolysis-dependent effect; the magnitude is illustrative (the caveat says so), and
    applying one forces Tier-D."""
    warns: List[str] = []
    vc = params.as_volumes_clearances()
    changed = False
    for mod in model.pharmacogenomic_modifiers:
        if not genotype_triggers(point, mod.phenotype_dimension, mod.phenotype_value):
            continue
        if model.drug_name not in mod.substrate_scope:
            continue  # substrate guardrail: not a substrate of this gene
        sf = mod.scale_factor
        if sf is None:
            warns.append(f"PGx MODIFIER [{mod.gene}]: no scale_factor curated — direction only "
                         f"(cited, Tier {mod.evidence_tier}); not applied numerically")
            continue
        affected = mod.affected_parameter
        targets = ("Cl1", "Cl2", "Cl3") if affected in ("Cl", "CL", "clearance") else (affected,)
        for sym in targets:
            if sym in vc and vc[sym]:
                vc[sym] = vc[sym] * sf
                changed = True
        warns.append(
            f"PGx MODIFIER [{mod.gene} {mod.phenotype_value}]: {affected} x{sf:g} (opt-in, "
            f"Tier D, forward-only) -> prolonged effect for a given dose"
            + (f". {mod.caveat}" if mod.caveat else "")
            + (f" [{mod.primary_citation}]" if mod.primary_citation else ""))
    if not changed:
        return params, warns
    mp = MicroParams.from_volumes_clearances(
        V1=vc["V1"], Cl1=vc["Cl1"], V2=vc.get("V2", 0.0), Cl2=vc.get("Cl2", 0.0),
        V3=vc.get("V3", 0.0), Cl3=vc.get("Cl3", 0.0), ke0=vc.get("ke0", 0.0))
    return mp, warns


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
    bands=False,
    percentile: Tuple[int, int] = (5, 95),
    samples: int = 2000,
    seed: Optional[int] = None,
    residual: bool = False,
    developmental: bool = False,
    pharmacogenomics: bool = False,
) -> SimulationResult:
    """Forward-simulate one PK (optionally + PD) model for one virtual patient.

    ``developmental=True`` (v0.8) opt-in carries an adult model into a child by the
    model's curated allometric size + PMA maturation block, forcing **Tier D** and
    warning. ``pharmacogenomics=True`` (v0.9) opt-in applies a declared genotype's
    cited, substrate-scoped kinetic modifier (also Tier-D). Both are off by default,
    never silently reparameterize a model, and never compute a dose (specs §5/§9).

    With ``bands=True`` (and a mandatory integer ``seed``), also draws a seeded
    Monte-Carlo prediction band from the model's curated between-subject variability
    (v0.2 §6). A model that publishes no BSV draws no band — the never-synthesize
    rule (§5): a missing band is honest; a borrowed one is a lie with error bars.

    ``bands`` may also be a list of kinds (v0.7 C2): ``["prediction", "covariate"]``
    adds the covariate band — covariate-equation choice always, plus a covariate-VALUE
    band when a covariate is supplied as a distribution (e.g. ``weight={"mean":70,"sd":6}``).
    A distribution-valued covariate is collapsed to its mean for every deterministic path.
    """
    kinds = _normalize_bands(bands)
    # A covariate may be supplied as a distribution dict; collapse to the point (mean)
    # vector for every deterministic path (kernel, envelope, PD). Scalar patients are
    # byte-identical, so this is backward-compatible.
    point = _point_patient(patient)
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
    if model.kernel_function not in KERNELS:
        raise NotImplementedError(
            f"{model_id} uses the '{model.kernel_function}' kernel, not an IV-disposition "
            "kernel that simulate() drives. Local anesthetics: use hypnos.la (site "
            "absorption); volatiles: use hypnos.mac/washin/washout."
        )
    kernel = KERNELS[model.kernel_function]

    # Opt-in forward transforms (off by default; never silently reparameterize). In
    # developmental mode the parameters are derived from the REFERENCE ADULT and then
    # size/maturation-scaled (v0.8 §5), so the adult covariate kernel is never run at
    # neonatal covariates (which would break, e.g. Schnider's LBM/height at 3.4 kg).
    extra_warnings: List[str] = []
    if developmental:
        params, dev_warns = _apply_developmental(model, kernel, point)
        extra_warnings.extend(dev_warns)
    else:
        params = kernel(point)
    if pharmacogenomics:
        params, pgx_warns = _apply_pgx_modifiers(model, params, point)
        extra_warnings.extend(pgx_warns)

    t = np.asarray(t, dtype=float)
    weight = float(point.get("weight", 70.0))
    dosing = build_dosing(schedule, weight)
    traj: Trajectory = _simulate_ref(params, dosing, t)

    drug_meta = ds.drug(model.drug_name) or {}
    tier_floor, warnings, excluded = evaluate_safety(model, point, drug_meta)
    warnings = list(warnings) + extra_warnings
    tier = worst_tier([model.tier, tier_floor])
    # An applied developmental/pgx-kinetic extrapolation forces Tier D by construction
    # (you cannot get an A-looking number out of an extrapolation; specs §5).
    if developmental or (pharmacogenomics and any("PGx MODIFIER" in w for w in extra_warnings)):
        tier = "D"

    result = SimulationResult(
        model_id=model_id, t=t, cp=traj.cp, ce=traj.ce, tier=tier,
        warnings=warnings, excluded=excluded, params=params, patient=dict(point),
        concentration_unit=drug_meta.get("concentration_unit", "ug/mL"),
    )
    # whether this model carries curated estimation uncertainty (drives the v0.3 §7
    # reducible/irreducible readout's `estimation_curated` flag)
    result.has_estimation = model.has_published_estimation

    pdm: Optional[Model] = None
    if pd_model is not None:
        pdm = ds[pd_model]
        if not pdm.kernel_implemented:
            raise NotImplementedError(f"{pd_model}: PD kernel not implemented")
        effect = _apply_pd(pdm, traj.ce, point)
        result.effect = effect
        result.effect_label = pdm.label
        result.pd_model_id = pd_model
        # composed simulation inherits the worst tier among PK, PD, and envelope floor
        result.tier = worst_tier([result.tier, pdm.tier])
        pd_floor, pd_warn, _ = evaluate_safety(pdm, point)
        result.tier = worst_tier([result.tier, pd_floor])
        result.warnings.extend(w for w in pd_warn if w not in result.warnings)

    if kinds:
        if seed is None:
            raise ValueError(
                "simulate(bands=...) requires an explicit integer seed — every "
                "band-producing call is seeded so quantiles are byte-reproducible (spec §6)."
            )
        if "prediction" in kinds:
            _attach_bands(result, model, params, dosing, t, percentile=percentile,
                          samples=samples, seed=seed, residual=residual,
                          pd_model=pdm, patient=point)
        if "confidence" in kinds:
            _attach_confidence_band(result, model, params, dosing, t, percentile=percentile,
                                    samples=samples, seed=seed)
        if "covariate" in kinds:
            _attach_covariate_band(result, model, kernel, dosing, t, patient=patient,
                                   point=point, percentile=percentile, samples=samples,
                                   seed=seed, pd_model=pdm, ds=ds)

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

    # v0.7 §7.3 — the fifth (covariate) component, present only when the covariate band
    # was computed (kinds include "covariate"); the legacy bands=True path leaves these
    # None, so the decomposition stays exactly 3-way (backward-compatible).
    def _cov_mean(suffix: str) -> Optional[float]:
        vals = [r.__dict__.get(f"{key}_{suffix}") for r in eligible]
        present = [v[t_star] for v in vals if v is not None]
        return float(np.mean(present)) if present else None

    eqn_var = _cov_mean("cov_equation_var")
    val_var = _cov_mean("cov_value_var")
    has_cov = eqn_var is not None or val_var is not None
    cov_var = (eqn_var or 0.0) + (val_var or 0.0)

    # v0.3 E2 — the estimation component, present only when a confidence band was drawn
    # (kinds include "confidence"); otherwise None and the decomposition is unchanged.
    est_var = _cov_mean("est_var")
    has_est = est_var is not None

    total = var_structural + bsv_var + resid_var + cov_var + (est_var or 0.0)
    if total > 0:
        struct_share = float(var_structural / total)
        bsv_share = float(bsv_var / total)
        resid_share = float(resid_var / total)
        cov_share = float(cov_var / total)
        est_share = float((est_var or 0.0) / total)
        vs = {
            "structural": round(struct_share, 4),
            "bsv": round(bsv_share, 4),
            "residual": round(resid_share, 4),
            "t_star_min": float(eligible[0].t[t_star]),
        }
        if has_est:
            vs["estimation"] = round(est_share, 4)
        if has_cov:
            vs["covariate"] = round(cov_share, 4)
        out["variance_share"] = vs
        # The v0.3 §7 reducible/irreducible decomposition — split the same variance
        # along the axis that decides what to DO about it. Structural (between-model) is
        # reducible by curating/validating more models; ESTIMATION is reducible with more
        # data per model; the v0.7 covariate component is reducible too (agree on the
        # equation / measure the covariate better); BSV and residual are irreducible (the
        # population is the limit; the assay is noisy).
        any_estimation = has_est or any(getattr(r, "has_estimation", False) for r in eligible)
        out["reducibility"] = {
            "reducible": round(struct_share + est_share + (cov_share if has_cov else 0.0), 4),
            "irreducible": round(bsv_share + resid_share, 4),
            "estimation_curated": bool(any_estimation),
            "note": ("reducible = between-model (more models/curation)"
                     + (" + estimation (more data per model)" if has_est else "")
                     + (" + covariate (agree on the equation / measure better)" if has_cov else "")
                     + ("" if has_est else
                        "; estimation component contributes 0 unless a confidence band is requested")),
        }
        if has_cov:
            # the two halves of §1, kept distinct: WHICH equation vs how well we know the value
            out["covariate_split"] = {
                "equation_choice": round(float((eqn_var or 0.0) / total), 4),
                "value_uncertainty": round(float((val_var or 0.0) / total), 4),
            }
            ctiers = [r.covariate_band_tier for r in eligible if r.covariate_band_tier]
            out["covariate_band_tier"] = worst_tier(ctiers) if ctiers else None

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
    bands=False,
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

    ``bands`` accepts the v0.7 C2 kinds too (``["prediction","covariate"]``): the variance
    decomposition then gains the **covariate** component, split into ``equation_choice`` and
    ``value_uncertainty`` (the two halves of §1), with a refined reducible/irreducible rollup.
    """
    from .filter import select

    kinds = _normalize_bands(bands)
    if kinds and seed is None:
        raise ValueError("compare(bands=...) requires an explicit integer seed (spec §6).")

    t = np.asarray(t, dtype=float)
    cmp = Comparison(drug=drug, purpose=purpose, t=t, bands=bool(kinds),
                     concentration_unit=(ds.drug(drug) or {}).get("concentration_unit", "ug/mL"))
    for m in select(ds, drug=drug, purpose=purpose):
        if not m.kernel_implemented:
            cmp.unavailable.append({"model_id": m.id, "tier": m.tier,
                                    "reason": "reference kernel pending verified transcription"})
            continue
        if m.kernel_function not in KERNELS:
            # not an IV-disposition model (e.g. a local anesthetic — use hypnos.la)
            cmp.unavailable.append({"model_id": m.id, "tier": m.tier,
                                    "reason": "not an IV-disposition model (see hypnos.la / volatiles)"})
            continue
        res = simulate(ds, m.id, patient=patient, schedule=schedule, t=t, pd_model=pd_model,
                       bands=bands, percentile=percentile, samples=samples, seed=seed,
                       residual=residual)
        if res.excluded:
            cmp.excluded.append({"model_id": m.id, "tier": res.tier, "reasons": res.warnings,
                                 "result": res})
        else:
            cmp.included.append(res)
            if "prediction" in kinds and not m.has_published_variability:
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
    if kinds:
        for k, rs in (("cp", cmp.included), ("ce", ce_results)):
            band_part = _band_divergence(rs, k)
            if band_part:
                cmp.divergence.setdefault(k, {}).update(band_part)
    return cmp


# --------------------------------------------------------------------------- #
# developmental_overlay — the v0.8 headline: fitted vs extrapolated vs visibly-wrong
# --------------------------------------------------------------------------- #
@dataclass
class DevelopmentalExtrapolation:
    """One adult model carried into a child by allometry (+ maturation), Tier-D (v0.8 §6)."""

    model_id: str
    basis: str                      # allometry_only | allometry_plus_maturation
    tier: str
    result: SimulationResult
    maturation_factor: Optional[float]
    caveat: Optional[str]


@dataclass
class DevelopmentalOverlay:
    """The developmental extrapolation overlay (v0.8 §6).

    ``fitted`` names every model with *actual* pediatric standing (in its band) or greyed
    (out of band, each reasoned); ``extrapolations`` are the Tier-D allometry(+maturation)
    bands for the adult models; ``linear_per_kg`` is the deliberately-labeled, visibly-wrong
    naive shortcut. The instrument's value is making the extrapolation honest, not
    resolving it — every curve here is Tier-D and labeled, and none is a dose."""

    drug: str
    patient: Dict[str, Any]
    t: np.ndarray
    fitted: List[Dict[str, Any]] = field(default_factory=list)
    extrapolations: List[DevelopmentalExtrapolation] = field(default_factory=list)
    linear_per_kg: Optional[SimulationResult] = None
    concentration_unit: str = "ug/mL"

    @property
    def conc_factor(self) -> float:
        return concentration_factor(self.concentration_unit)


def _linear_per_kg_result(model: Model, point: Dict[str, Any], schedule, t: np.ndarray,
                          cu: str) -> SimulationResult:
    """The naive linear-per-kg adult rescaling as a labeled, visibly-wrong reference (v0.8 §6)."""
    kernel = KERNELS[model.kernel_function]
    dev = model.developmental_model
    ref_wt = dev.size.reference_weight_kg if dev and dev.size else 70.0
    vc = kernel(_reference_adult_patient(str(point.get("sex", "M")), ref_wt)).as_volumes_clearances()
    weight = float(point.get("weight", 70.0))
    lin = linear_per_kg_scale(vc, weight, reference_weight_kg=ref_wt)
    mp = MicroParams.from_volumes_clearances(
        V1=lin["V1"], Cl1=lin["Cl1"], V2=lin.get("V2", 0.0), Cl2=lin.get("Cl2", 0.0),
        V3=lin.get("V3", 0.0), Cl3=lin.get("Cl3", 0.0), ke0=lin.get("ke0", 0.0))
    traj = _simulate_ref(mp, build_dosing(schedule, weight), t)
    return SimulationResult(
        model_id=f"{model.id} (linear-per-kg shortcut)", t=t, cp=traj.cp, ce=traj.ce, tier="D",
        warnings=["LINEAR-PER-KG REFERENCE: the naive adult rescaling (every disposition "
                  "parameter ∝ weight). Shown ONLY as a labeled, visibly-wrong reference: it "
                  "ignores the allometric ¾-power scaling of clearance AND the PMA maturation of "
                  "clearance — the two mechanisms (acting in opposite directions) that make a "
                  "neonate not a small adult. The divergence between this line and the mechanistic "
                  "extrapolation IS the teaching point; it is NEVER a dose (v0.8 §6)."],
        params=mp, patient=dict(point), concentration_unit=cu)


def developmental_overlay(
    ds: Dataset,
    *,
    drug: str,
    patient: Dict[str, Any],
    schedule: Optional[Schedule] = None,
    t: Optional[np.ndarray] = None,
    purpose: str = "pk",
) -> DevelopmentalOverlay:
    """The v0.8 headline: in this child, which models have fitted standing, how far does the
    allometry(+maturation) extrapolation of the adult models spread, and how badly does the
    per-kg shortcut miss (v0.8 §6)?

    The eligible-model set shrinks to those with actual pediatric standing (greyed + reasoned
    where out of band); every greyed adult model with a curated developmental block is fanned
    out as a Tier-D extrapolation; the naive linear-per-kg line is overlaid as a labeled,
    visibly-wrong reference. No curve is a dose."""
    from .filter import select
    from .presets import default_schedule_for

    point = _point_patient(patient)
    schedule = schedule if schedule is not None else default_schedule_for(drug)
    t = np.asarray(t if t is not None else np.linspace(0.0, 30.0, 300), dtype=float)
    cu = (ds.drug(drug) or {}).get("concentration_unit", "ug/mL")
    overlay = DevelopmentalOverlay(drug=drug, patient=dict(point), t=t, concentration_unit=cu)

    rep_for_linear: Optional[Model] = None
    for m in select(ds, drug=drug, purpose=purpose):
        if not m.kernel_implemented or m.kernel_function not in KERNELS:
            continue
        floor, warns, excluded = evaluate_safety(m, point, ds.drug(m.drug_name) or {})
        if m.is_fitted_pediatric:
            overlay.fitted.append({
                "model_id": m.id, "tier": worst_tier([m.tier, floor]),
                "in_standing": not excluded, "reasons": warns})
            continue
        if m.has_developmental_model:
            res = simulate(ds, m.id, patient=point, schedule=schedule, t=t, developmental=True)
            dev = m.developmental_model
            overlay.extrapolations.append(DevelopmentalExtrapolation(
                model_id=m.id, basis=dev.extrapolation_basis, tier=res.tier, result=res,
                maturation_factor=maturation_value(dev, point), caveat=dev.caveat))
            if rep_for_linear is None:
                rep_for_linear = m
        else:
            overlay.fitted.append({
                "model_id": m.id, "tier": worst_tier([m.tier, floor]), "in_standing": False,
                "reasons": warns + ["no developmental_model curated -> mechanistic extrapolation "
                                    "not available (an honest gap, never a fabricated curve)"]})

    if rep_for_linear is not None:
        overlay.linear_per_kg = _linear_per_kg_result(rep_for_linear, point, schedule, t, cu)
    return overlay


# --------------------------------------------------------------------------- #
# covariate_divergence — divergence WITHIN one model, over covariate equations (v0.7 C1)
# --------------------------------------------------------------------------- #
@dataclass
class EquationCurve:
    """One model curve computed with a substituted body-size equation (v0.7 §7.1)."""

    equation_id: str
    quantity: str                  # the model's derived-input quantity the value plugs into
    derived_value: float           # the equation's value (kg) for this patient
    verbatim: bool                 # True => the model's OWN equation; False => a substitution
    in_envelope: bool              # within the equation's own validity envelope
    inverted: bool                 # the James-style non-physical inversion
    tier: str
    cp: np.ndarray
    ce: np.ndarray
    status: str

    # alias so the shared `_divergence` machinery can name the driver by equation id
    @property
    def model_id(self) -> str:
        return self.equation_id

    @property
    def cp_peak(self) -> float:
        return float(np.max(self.cp))

    @property
    def ce_peak(self) -> float:
        return float(np.max(self.ce))


@dataclass
class CovariateDivergence:
    model_id: str
    derived_equation: str          # the model's own equation (verbatim)
    quantity: str
    patient: Dict[str, Any]
    t: np.ndarray
    key: str                       # "ce" (effect-site model) or "cp" (PK-only)
    by_equation: List[EquationCurve] = field(default_factory=list)
    divergence: Dict[str, Any] = field(default_factory=dict)
    concentration_unit: str = "ug/mL"

    @property
    def conc_factor(self) -> float:
        return concentration_factor(self.concentration_unit)

    @property
    def own(self) -> EquationCurve:
        return next(c for c in self.by_equation if c.verbatim)


# body-size descriptors that are mutually substitutable into one covariate slot — the
# documented pump-to-pump substitution (e.g. a TCI implementing "Schnider" with
# Janmahasatian FFM instead of James LBM) is exactly an LBM<->FFM swap (v0.7 §1).
_SUBSTITUTABLE_QUANTITIES = ("lbm", "ffm", "nfm", "ibw")


def covariate_divergence(
    ds: Dataset,
    model_id: str,
    *,
    patient: Dict[str, Any],
    schedule: Optional[Schedule] = None,
    t: Optional[np.ndarray] = None,
    candidates: Optional[List[str]] = None,
) -> CovariateDivergence:
    """Divergence *within* one model: overlay its predicted curve under each admissible
    covariate equation, greying any equation outside its own validity envelope (v0.7 §7.1).

    v0.1's :func:`compare` overlays *different models*; this overlays the *same model* under
    *different body-size equations* — the substitution some TCI pumps make silently. The
    model's own equation is marked ``verbatim``; substitutions are ``verbatim=False`` and
    shown only here, never as the model's own prediction (the never-invent rule, §10). When
    the model's own equation has inverted (James above BMI ≈ 37), this view shows the
    documented Schnider-in-obesity failure mode *at its covariate source*.
    """
    from .presets import default_schedule_for

    patient = _point_patient(patient)   # equation substitution is at a fixed covariate vector
    model = ds[model_id]
    cm = model.covariate_model
    if cm is None:
        raise ValueError(
            f"{model_id} declares no covariate_model — there is no covariate-equation "
            "divergence axis (it scales on raw covariates, or the binding is uncurated)"
        )
    if not model.kernel_implemented or model.kernel_function not in KERNELS:
        raise NotImplementedError(
            f"{model_id}: covariate_divergence needs an implemented IV-disposition kernel"
        )

    di = cm.derived_inputs[0]          # every curated covariate-scaled model has exactly one
    quantity = di.quantity
    own = di.equation

    lib = ds.covariate_equations
    if candidates is None:
        candidates = [eid for eid, rec in lib.items()
                      if rec.get("quantity") in _SUBSTITUTABLE_QUANTITIES]
    # the model's own equation always leads, substitutions follow (stable order)
    ordered = [own] + sorted(c for c in candidates if c != own)

    drug = model.drug_name
    schedule = schedule if schedule is not None else default_schedule_for(drug)
    t = np.asarray(t if t is not None else np.linspace(0.0, 60.0, 361), dtype=float)
    weight = float(patient.get("weight", 70.0))
    dosing = build_dosing(schedule, weight)
    kernel = KERNELS[model.kernel_function]

    curves: List[EquationCurve] = []
    for eid in ordered:
        if eid not in lib:
            continue
        r = _cov_eval(eid, patient, ds=ds)
        params = kernel({**patient, f"_{quantity}_override": r.value})
        traj = _simulate_ref(params, dosing, t)
        verbatim = eid == own
        if not r.out_of_envelope:
            status = "in-envelope" + ("" if verbatim else " (substitution; verbatim=false)")
        elif r.inverted:
            status = "OUTSIDE equation envelope — INVERSION (greyed)"
        else:
            status = "OUTSIDE equation envelope (greyed)"
        curves.append(EquationCurve(
            equation_id=eid, quantity=r.quantity, derived_value=r.value,
            verbatim=verbatim, in_envelope=not r.out_of_envelope, inverted=r.inverted,
            tier=r.tier, cp=traj.cp, ce=traj.ce, status=status,
        ))

    key = "ce" if (curves and curves[0].ce.max() > 0) else "cp"
    div = _divergence(curves, key) if len(curves) >= 2 else {}
    return CovariateDivergence(
        model_id=model_id, derived_equation=own, quantity=quantity, patient=dict(patient),
        t=t, key=key, by_equation=curves, divergence=div,
        concentration_unit=(ds.drug(drug) or {}).get("concentration_unit", "ug/mL"),
    )


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


# --------------------------------------------------------------------------- #
# simulate_reversal — antagonism as an interaction (v0.5 Part A)
# --------------------------------------------------------------------------- #
@dataclass
class ReversalResult:
    """Forward simulation of a reversal agent acting on an agonist already in the patient (v0.5 §A4)."""

    reversal_id: str
    agonist_id: str
    mechanism: str
    t: np.ndarray
    free_ce: np.ndarray                          # free (active) agonist effect-site conc
    tier: str
    warnings: List[str] = field(default_factory=list)
    effect: Optional[np.ndarray] = None          # composed effect (TOF / opioid), if a PD model given
    effect_label: Optional[str] = None
    free_ce_no_reversal: Optional[np.ndarray] = None   # the same agonist with NO reversal (the contrast)
    effect_no_reversal: Optional[np.ndarray] = None
    complex_umol: Optional[np.ndarray] = None    # encapsulation: inert bound complex
    ce50_multiplier: Optional[np.ndarray] = None  # competitive: apparent-Ce50 shift over time
    renarcotization: bool = False                # competitive: did the effect relapse?


def _antagonist_micro(model: Model) -> MicroParams:
    """Build the antagonist's own MicroParams from its curated disposition parameters."""
    pv = {p.symbol: p.central for p in model.parameters}
    return MicroParams.from_volumes_clearances(
        V1=pv.get("V1", 1.0), Cl1=pv.get("Cl1", 0.0),
        V2=pv.get("V2", 0.0), Cl2=pv.get("Cl2", 0.0),
        V3=pv.get("V3", 0.0), Cl3=pv.get("Cl3", 0.0), ke0=pv.get("ke0", 0.0))


def simulate_reversal(
    ds: Dataset,
    reversal_id: str,
    *,
    patient: Dict[str, Any],
    reversal_dose: str,
    agonist_dose: str,
    t: np.ndarray,
    agonist_id: Optional[str] = None,
    reversal_time_min: float = 10.0,
    pd_model: Optional[str] = None,
) -> ReversalResult:
    """Simulate a *given* reversal dose acting on an agonist already in the patient (v0.5 §A4).

    Dispatches on the reversal mechanism (encapsulation / competitive_receptor /
    enzyme_inhibition). Forward-only: it renders the trajectory of the dose you specify and the
    no-reversal contrast — it **never computes the reversal dose** ("how much sugammadex/
    naloxone?"), which is inverse control wearing an antidote's coat (§A6)."""
    from .reversal import competitive_shift, encapsulation, indirect_inhibition

    rev = ds[reversal_id]
    rm = rev.reversal_mechanism
    if rm is None:
        raise ValueError(f"{reversal_id} is not a reversal agent (no reversal_mechanism block)")
    point = _point_patient(patient)
    agonist_id = agonist_id or (rm.get("targets") or [None])[0]
    if agonist_id is None or agonist_id not in ds:
        raise ValueError(f"{reversal_id}: reversal target {agonist_id!r} does not resolve")
    agonist = ds[agonist_id]
    t = np.asarray(t, dtype=float)
    weight = float(point.get("weight", 70.0))
    rev_mg = parse_amount(reversal_dose, weight)
    ago_mg = parse_amount(agonist_dose, weight)
    warns = [f"REVERSAL ({rm['type']}): simulating a GIVEN dose ({reversal_dose}) of {rev.drug_name} "
             f"on {agonist.drug_name} — forward-only; Hypnos never computes a reversal dose (§A6)"]
    mech = rm["type"]

    if mech == "encapsulation":
        if not (agonist.kernel_implemented and agonist.kernel_function in KERNELS):
            raise NotImplementedError(f"{agonist_id}: agonist PK kernel not implemented")
        nmb = KERNELS[agonist.kernel_function](point)
        no = encapsulation(nmb, ago_mg, 0.0, t, sugammadex_time_min=reversal_time_min)
        rv = encapsulation(nmb, ago_mg, rev_mg, t, sugammadex_time_min=reversal_time_min)
        appr = (rm.get("binding") or {}).get("approximation")
        if appr:
            warns.append(f"ENCAPSULATION: 1:1 molar binding via the '{appr}' SIMPLIFICATION "
                         "(binding kinetics not curated) — Tier-C, labeled (§A5)")
        res = ReversalResult(
            reversal_id=reversal_id, agonist_id=agonist_id, mechanism=mech, t=t,
            free_ce=rv.free_ce, free_ce_no_reversal=no.free_ce, complex_umol=rv.complex_umol,
            tier=worst_tier([rev.tier, agonist.tier]), warnings=warns)
        if pd_model is not None:
            pdm = ds[pd_model]
            res.effect = _apply_pd(pdm, rv.free_ce, point)
            res.effect_no_reversal = _apply_pd(pdm, no.free_ce, point)
            res.effect_label = pdm.label
            res.tier = worst_tier([res.tier, pdm.tier])
        return res

    if mech == "competitive_receptor":
        # the antagonist's own effect-site concentration over time, then the apparent-Ce50 shift
        amic = _antagonist_micro(rev)
        atraj = _simulate_ref(amic, build_dosing([("bolus", reversal_time_min, reversal_dose)], weight), t)
        antag = rm.get("antagonism") or {}
        ki = float(antag.get("ki", 1.0))
        # internal concentrations are µg/mL (== mg/L); convert Ki to the same units so the
        # competitive ratio C/Ki is dimensionally correct (a ng/mL Ki is 1/1000 of a µg/mL).
        if str(antag.get("ki_units", "")).lower() == "ng/ml":
            ki = ki / 1000.0
        mult = competitive_shift(atraj.ce, ki)
        # renarcotization: the multiplier rises then FALLS back toward 1 as the antagonist clears
        i_peak = int(np.argmax(mult))
        relapsed = bool(i_peak < len(mult) - 1 and mult[-1] < 0.5 * (mult[i_peak] - 1.0) + 1.0)
        if relapsed:
            warns.append("RENARCOTIZATION: the apparent-Ce50 shift peaks then RELAPSES toward 1 as "
                         f"{rev.drug_name} clears — the antagonist is outlasted by the agonist; the "
                         "reversed effect returns (rendered forward, never a re-dosing schedule; §A5)")
        return ReversalResult(
            reversal_id=reversal_id, agonist_id=agonist_id, mechanism=mech, t=t,
            free_ce=atraj.ce, ce50_multiplier=mult, renarcotization=relapsed,
            tier=worst_tier([rev.tier, agonist.tier]), warnings=warns)

    if mech == "enzyme_inhibition":
        ind = rm.get("indirect") or {}
        # demonstrate the ceiling on a representative block-depth sweep (1=full block .. 0=recovered)
        block = np.linspace(1.0, 0.0, len(t))
        residual = indirect_inhibition(block, ceiling=float(ind.get("ceiling", 0.5)),
                                       depth_floor=float(ind.get("depth_floor", 0.9)))
        warns.append("CEILING: neostigmine cannot reverse a block deeper than the depth floor "
                     f"({ind.get('depth_floor')}); below it the residual block is unchanged (give "
                     "sugammadex). Muscarinic effects require an antimuscarinic co-agent (§A5).")
        return ReversalResult(
            reversal_id=reversal_id, agonist_id=agonist_id, mechanism=mech, t=t,
            free_ce=block, effect=residual, effect_label="residual block fraction",
            tier=worst_tier([rev.tier, agonist.tier]), warnings=warns)

    raise ValueError(f"unknown reversal mechanism {mech!r}")
