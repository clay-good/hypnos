"""Derived PK/PD characterizations (forward-only, safe by construction).

Two families live here:

* **Onset / offset** — **time to peak effect** (``tpeak``, the spec's
  ``effect_link`` §3) and the forward **decrement time** after a fixed-rate
  infusion. After a bolus the effect-site concentration peaks when its rate of
  change is zero, i.e. exactly when ``Ce == Cp``; ``tpeak`` characterizes a
  model's onset and is dose-independent (the system is linear).

* **External-validation metrics** (v0.4 §6) — **Varvel's** performance error and
  the four derived metrics (MDPE, MDAPE, wobble, divergence) that quantify how
  well a model's *predicted* concentration/effect matches *observed* data, plus a
  seeded population roll-up. ``validate_against_cohort`` drives the existing
  forward solver from each subject's recorded dose history and compares to the
  recorded observations. This is the engine v0.4 builds; the source-specific data
  adapters (Open-TCI, VitalDB) sit on top of it and are not part of this layer.

**Why no context-sensitive half-time (CSHT).** The classic CSHT is the time for
plasma to fall 50% after a target-controlled infusion that held plasma
*constant* for some duration. Computing it requires solving for the infusion that
maintains a target concentration, which is inverse control: exactly the step
Hypnos refuses to take (spec §10). Onset (``tpeak``) and the decline after a
*fixed* dose history are forward problems and in scope; constant-concentration
CSHT is not. The validation engine is likewise forward-only: it runs the
*recorded* dose history and never searches for a dose (v0.4 §10).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .load import Dataset
from .models import worst_tier
from .reference import Dosing, simulate as _simulate_ref


@dataclass
class DecrementTime:
    model_id: str
    infusion: str
    duration_min: float
    fraction: float
    conc_at_stop: float          # plasma conc at the moment infusion stops (internal ug/mL)
    decrement_min: float         # time after stop for plasma to fall by `fraction` (inf if not reached)
    tier: str
    warnings: List[str] = field(default_factory=list)


@dataclass
class PeakEffect:
    model_id: str
    tpeak_min: float                 # time of the effect-site peak after a bolus
    ke0: float                       # the effect-site equilibration rate constant
    ce_cp_ratio_at_peak: float       # ~1.0 by definition (Ce = Cp at the peak); a sanity check
    tier: str
    warnings: List[str] = field(default_factory=list)


def time_to_peak_effect(
    ds: Dataset, model_id: str, *, patient: Dict[str, Any], tmax: float = 20.0, n: int = 8001
) -> PeakEffect:
    """Time to peak effect-site concentration after a bolus (onset characterization).

    Forward-only and dose-independent. Requires a PK model with an effect
    compartment (ke0 > 0); raises for PK-only models (e.g. Kim remifentanil,
    Paedfusor) where effect-site onset is undefined.
    """
    from .covariates import point_patient
    from .export.registry import KERNELS
    from .simulate import evaluate_safety

    patient = point_patient(patient)   # collapse any covariate distribution to its mean (v0.7 C2)
    model = ds[model_id]
    if model.purpose != "pk":
        raise ValueError(f"time_to_peak_effect expects a PK model; {model_id} is '{model.purpose}'")
    if not model.kernel_implemented or model.kernel_function not in KERNELS:
        raise NotImplementedError(f"{model_id}: reference kernel not implemented")
    params = KERNELS[model.kernel_function](patient)
    if params.ke0 <= 0:
        raise ValueError(
            f"{model_id} has no effect compartment (ke0 = 0); time-to-peak-effect is undefined"
        )

    t = np.linspace(0.0, tmax, n)
    traj = _simulate_ref(params, Dosing(boluses=((0.0, 1.0),)), t)  # unit bolus; tpeak is dose-independent
    i = int(np.argmax(traj.ce))

    tier_floor, warnings, _ = evaluate_safety(model, patient)
    warnings = list(warnings)
    if i >= n - 1:
        warnings.append(f"WARNING: effect-site peak not reached within tmax={tmax:g} min; increase tmax")
    ratio = float(traj.ce[i] / traj.cp[i]) if traj.cp[i] > 0 else float("nan")

    return PeakEffect(
        model_id=model_id, tpeak_min=float(t[i]), ke0=params.ke0,
        ce_cp_ratio_at_peak=ratio, tier=worst_tier([model.tier, tier_floor]),
        warnings=warnings,
    )


def decrement_time(
    ds: Dataset, model_id: str, *, patient: Dict[str, Any], infusion: str,
    duration: float, fraction: float = 0.5,
) -> DecrementTime:
    """Plasma decrement time after a **constant-rate** infusion of given duration.

    Runs a fixed-rate infusion for ``duration`` minutes, stops it, and returns the
    time for the plasma concentration to fall by ``fraction`` (default 50%) from
    its value at the stop. Forward-only and a function of the *fixed dose history*.

    This is **not** the classic context-sensitive half-time, which is defined for a
    target-controlled infusion holding plasma *constant* — that requires inverse
    control, which Hypnos does not do (spec §10). This metric captures the same
    context-sensitivity qualitatively (it lengthens with infusion duration for
    accumulating drugs, and is near-flat for remifentanil) from a forward,
    constant-rate regimen.
    """
    from .covariates import point_patient
    from .export.registry import instantiate
    from .simulate import build_dosing, evaluate_safety

    if not 0.0 < fraction < 1.0:
        raise ValueError("fraction must be in (0, 1)")
    patient = point_patient(patient)   # collapse any covariate distribution to its mean (v0.7 C2)
    model = ds[model_id]
    if model.purpose != "pk":
        raise ValueError(f"decrement_time expects a PK model; {model_id} is '{model.purpose}'")
    params = instantiate(model, patient)  # raises NotImplementedError if kernel pending

    weight = float(patient.get("weight", 70.0))
    schedule = [("infusion", 0.0, infusion), ("infusion", float(duration), "0 mg/kg/h")]
    dosing = build_dosing(schedule, weight)

    horizon = duration + max(120.0, 6.0 * duration)
    grid = np.linspace(0.0, horizon, int(horizon * 4) + 1)
    t = np.union1d(grid, [float(duration)])  # ensure the stop instant is sampled
    traj = _simulate_ref(params, dosing, t)

    i_stop = int(np.searchsorted(t, duration))
    c_stop = float(traj.cp[i_stop])
    tier_floor, warnings, _ = evaluate_safety(model, patient)
    warnings = list(warnings)

    decrement = float("inf")
    if c_stop <= 0:
        warnings.append("WARNING: plasma concentration at stop is ~0; decrement undefined")
    else:
        threshold = (1.0 - fraction) * c_stop
        after = np.where((t >= duration) & (traj.cp <= threshold))[0]
        if len(after):
            decrement = float(t[after[0]] - duration)
        else:
            warnings.append(f"WARNING: {100*fraction:g}% decrement not reached within {horizon:g} min")

    return DecrementTime(
        model_id=model_id, infusion=infusion, duration_min=float(duration),
        fraction=fraction, conc_at_stop=c_stop, decrement_min=decrement,
        tier=worst_tier([model.tier, tier_floor]), warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# External-validation metrics — Varvel's framework (v0.4 §6)
# --------------------------------------------------------------------------- #
#
# Given a series of observed vs. predicted concentrations (or effects), Varvel's
# canonical anesthesia PK/PD validation methodology computes, per sample j of
# subject i:
#
#   performance error   PE_ij = 100 * (C_obs,ij - C_pred,ij) / C_pred,ij   (%)
#
# and, per subject, four summaries:
#
#   MDPE   = median_j(PE_ij)               — bias        (signed)
#   MDAPE  = median_j(|PE_ij|)             — inaccuracy  (magnitude)
#   wobble = median_j(|PE_ij - MDPE_i|)    — intra-individual variability
#   divergence = slope of |PE_ij| vs t_j   — drift of error with time (%/h)
#
# The population roll-up is the median (with a seeded nonparametric bootstrap CI)
# of each across subjects. C_pred is exactly what the reference kernels already
# produce; this module adds only the alignment + the metric math (pure NumPy).


@dataclass(frozen=True)
class VarvelResult:
    """The four Varvel metrics for one subject (or one pooled observation set)."""

    mdpe: float          # median performance error (%), signed — bias
    mdape: float         # median absolute performance error (%) — inaccuracy
    wobble: float        # median |PE - MDPE| (%) — intra-individual variability
    divergence: float    # slope of |PE| vs time (%/h) — drift; nan if < 2 valid samples
    n: int               # number of valid (finite-PE) samples that entered the metrics


@dataclass(frozen=True)
class PopulationPerformance:
    """Population roll-up of subject-level Varvel metrics with seeded bootstrap CIs."""

    mdpe: float
    mdape: float
    wobble: float
    divergence: float
    n_subjects: int
    ci95: Dict[str, Tuple[float, float]]   # metric name -> (low, high)
    seed: int


@dataclass(frozen=True)
class SubjectRecord:
    """One subject's recorded case: covariates, dose history, and observations.

    This is the common shape every source-specific adapter (Open-TCI, VitalDB)
    maps its native records to; the metric math downstream is adapter-agnostic.
    ``observations`` are ``(time_min, value, kind)`` with ``kind`` naming the
    measured quantity (``"cp"``, ``"ce"``, ``"bis"``, ``"tof"``).
    """

    covariates: Dict[str, Any]
    schedule: Sequence[Tuple[str, float, str]]
    observations: Sequence[Tuple[float, float, str]]
    subject_id: Optional[str] = None


@dataclass
class CohortValidation:
    """Hypnos-computed external validation of one model against a cohort (v0.4 §4.1)."""

    model_id: str
    dataset: str
    mode: str            # pk_concentration | pd_bis | pd_tof
    target: str          # cp | ce | bis | tof
    n_subjects: int
    population: PopulationPerformance
    per_subject: List[VarvelResult] = field(default_factory=list)
    in_envelope: Optional[bool] = None
    seed: int = 0

    def to_record(self) -> Dict[str, Any]:
        """Serialize to a schema ``external_validation[]`` entry (v0.4 §4.1).

        The provenance block is intentionally partial here: ``hypnos_version`` and
        ``git_commit`` are stamped by the caller that commits the artifact, never
        by wall-clock or environment reads inside the engine (v0.4 §3 determinism).
        """
        p = self.population
        metrics = [
            {"name": "MDPE", "value": _round(p.mdpe), "units": "%",
             "ci95": _ci(p.ci95.get("mdpe"))},
            {"name": "MDAPE", "value": _round(p.mdape), "units": "%",
             "ci95": _ci(p.ci95.get("mdape"))},
            {"name": "wobble", "value": _round(p.wobble), "units": "%",
             "ci95": _ci(p.ci95.get("wobble"))},
            {"name": "divergence", "value": _round(p.divergence), "units": "%/h",
             "ci95": _ci(p.ci95.get("divergence"))},
        ]
        return {
            "dataset": self.dataset,
            "mode": self.mode,
            "target": self.target,
            "cohort": {"n_subjects": self.n_subjects, "filter": "all",
                       "in_envelope": self.in_envelope},
            "metrics": metrics,
            "provenance": {"computed_by": "hypnos", "seed": self.seed},
            "reproducible": True,
        }


def _nanmedian(col: np.ndarray) -> float:
    """nan-aware median that returns nan (quietly) for an all-nan column."""
    if not np.isfinite(col).any():
        return float("nan")
    return float(np.nanmedian(col))


def _round(x: Optional[float], nd: int = 4) -> Optional[float]:
    if x is None or not np.isfinite(x):
        return None
    return float(round(x, nd))


def _ci(pair: Optional[Tuple[float, float]]) -> Optional[Dict[str, Optional[float]]]:
    if pair is None:
        return None
    return {"low": _round(pair[0]), "high": _round(pair[1])}


def performance_error(c_obs: Sequence[float], c_pred: Sequence[float]) -> np.ndarray:
    """Elementwise Varvel performance error ``PE = 100*(obs - pred)/pred`` (%).

    The prediction is the denominator, so where it is non-positive the PE is
    undefined; those samples become ``nan`` and are dropped by the metrics (the
    near-zero-prediction guard the metric requires — v0.4 §6). Returns an array
    the same shape as the inputs.
    """
    obs = np.asarray(c_obs, dtype=float)
    pred = np.asarray(c_pred, dtype=float)
    if obs.shape != pred.shape:
        raise ValueError(
            f"c_obs and c_pred must have the same shape ({obs.shape} != {pred.shape})"
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        pe = 100.0 * (obs - pred) / pred
    pe = np.where(pred > 0, pe, np.nan)
    return pe


def varvel_metrics(pe: Sequence[float], times: Optional[Sequence[float]] = None) -> VarvelResult:
    """MDPE, MDAPE, wobble, and divergence for one subject's performance errors.

    ``pe`` is the per-sample PE% (e.g. from :func:`performance_error`); ``nan``
    entries are ignored. ``times`` (minutes) are required for ``divergence`` (the
    least-squares slope of ``|PE|`` vs time, reported in %/h); without them, or
    with fewer than two distinct valid times, ``divergence`` is ``nan``.
    """
    pe = np.asarray(pe, dtype=float)
    valid = np.isfinite(pe)
    pev = pe[valid]
    n = int(pev.size)
    if n == 0:
        return VarvelResult(float("nan"), float("nan"), float("nan"), float("nan"), 0)
    mdpe = float(np.median(pev))
    mdape = float(np.median(np.abs(pev)))
    wobble = float(np.median(np.abs(pev - mdpe)))
    divergence = float("nan")
    if times is not None:
        t = np.asarray(times, dtype=float)
        if t.shape != valid.shape:
            raise ValueError("times must match the shape of pe")
        tv = t[valid]
        if tv.size >= 2 and float(np.ptp(tv)) > 0:
            # slope of |PE| (%) vs time (min) -> %/min, reported as %/h.
            slope_per_min = float(np.polyfit(tv, np.abs(pev), 1)[0])
            divergence = slope_per_min * 60.0
    return VarvelResult(mdpe, mdape, wobble, divergence, n)


def pooled_performance(
    per_subject: Sequence[VarvelResult], *, seed: int, n_boot: int = 2000
) -> PopulationPerformance:
    """Population median of each Varvel metric across subjects, with seeded CIs.

    The point estimate is the (nan-aware) median across subjects; the 95% CI is a
    seeded nonparametric bootstrap over subjects (resample subjects with
    replacement, recompute the median, take the 2.5/97.5 percentiles). Identical
    ``(per_subject, seed, n_boot)`` -> identical CIs (v0.4 §3 determinism). A
    subject whose divergence is ``nan`` (a single sample) is simply ignored by the
    nan-aware median for that one metric, never imputed.
    """
    subs = [r for r in per_subject if r.n > 0]
    k = len(subs)
    names = ("mdpe", "mdape", "wobble", "divergence")
    if k == 0:
        nanci = {name: (float("nan"), float("nan")) for name in names}
        return PopulationPerformance(
            float("nan"), float("nan"), float("nan"), float("nan"), 0, nanci, seed
        )
    cols = {name: np.array([getattr(r, name) for r in subs], dtype=float) for name in names}
    point = {name: _nanmedian(cols[name]) for name in names}

    rng = np.random.default_rng(seed)
    ci: Dict[str, Tuple[float, float]] = {}
    if k == 1:
        # one subject: the median is that subject; a bootstrap CI is degenerate.
        ci = {name: (point[name], point[name]) for name in names}
    else:
        idx = rng.integers(0, k, size=(n_boot, k))
        for name in names:
            col = cols[name]
            if not np.isfinite(col).any():   # all-nan metric (e.g. divergence, no times)
                ci[name] = (float("nan"), float("nan"))
                continue
            boot = np.nanmedian(col[idx], axis=1)
            boot = boot[np.isfinite(boot)]
            if boot.size:
                lo, hi = np.percentile(boot, [2.5, 97.5])
                ci[name] = (float(lo), float(hi))
            else:
                ci[name] = (float("nan"), float("nan"))
    return PopulationPerformance(
        point["mdpe"], point["mdape"], point["wobble"], point["divergence"], k, ci, seed
    )


_TARGET_MODE = {"cp": "pk_concentration", "ce": "pk_concentration",
                "bis": "pd_bis", "tof": "pd_tof"}


def validate_against_cohort(
    ds: Dataset,
    model_id: str,
    subjects: Sequence[SubjectRecord],
    *,
    target: str = "cp",
    pd_model: Optional[str] = None,
    dataset: str = "in_memory",
    seed: int = 0,
    grid_min_points: int = 200,
) -> CohortValidation:
    """Run a model forward against a cohort and compute its Varvel metrics (v0.4 §6).

    For each subject: simulate the model from the *recorded* dose history and
    covariates, interpolate the prediction to the observation timestamps, compute
    the performance error against the observed values, then the per-subject
    metrics. The population roll-up aggregates across subjects with a seeded
    bootstrap CI. Forward-only and deterministic — it never tunes a dose (v0.4 §10).

    ``target`` selects the predicted quantity (``cp``/``ce`` plasma/effect-site
    concentration, ``bis`` depth-of-anaesthesia effect — which needs ``pd_model``).
    Only observations whose ``kind`` matches ``target`` are scored. The
    source-specific adapters that produce :class:`SubjectRecord`\\ s live above this
    function; this is the adapter-agnostic engine.
    """
    from .simulate import simulate as _simulate

    if target not in _TARGET_MODE:
        raise ValueError(f"target must be one of {sorted(_TARGET_MODE)}; got {target!r}")
    if target in ("bis", "tof") and pd_model is None:
        raise ValueError(f"target={target!r} requires a pd_model (PK->effect->{target.upper()} stack)")

    per_subject: List[VarvelResult] = []
    for subj in subjects:
        obs = [(float(t), float(v)) for (t, v, kind) in subj.observations if kind == target]
        if not obs:
            continue
        obs_t = np.array([t for t, _ in obs], dtype=float)
        obs_v = np.array([v for _, v in obs], dtype=float)
        # grid spans the observation window and includes the observation instants.
        tmax = float(obs_t.max())
        grid = np.linspace(0.0, max(tmax, 1e-6), max(grid_min_points, 2))
        t = np.union1d(grid, obs_t)
        res = _simulate(ds, model_id, patient=subj.covariates,
                        schedule=list(subj.schedule), t=t, pd_model=pd_model)
        pred_curve = {"cp": res.cp, "ce": res.ce, "bis": res.effect,
                      "tof": res.effect}[target]
        if pred_curve is None:
            raise ValueError(
                f"{model_id}: no '{target}' prediction available "
                f"(did you pass pd_model for an effect target?)"
            )
        pred_at_obs = np.interp(obs_t, t, pred_curve)
        pe = performance_error(obs_v, pred_at_obs)
        per_subject.append(varvel_metrics(pe, obs_t))

    pop = pooled_performance(per_subject, seed=seed)
    return CohortValidation(
        model_id=model_id, dataset=dataset, mode=_TARGET_MODE[target], target=target,
        n_subjects=len(per_subject), population=pop, per_subject=per_subject, seed=seed,
    )


# --------------------------------------------------------------------------- #
# Envelope stratification + cross-model leaderboard (v0.4 VE2/VE3)
# --------------------------------------------------------------------------- #
def partition_by_envelope(
    ds: Dataset, model_id: str, subjects: Sequence[SubjectRecord]
) -> Tuple[List[SubjectRecord], List[SubjectRecord]]:
    """Split a cohort into (in-envelope, out-of-envelope) subjects for one model (v0.4 §4.2).

    A subject is *in envelope* when its covariates raise no
    :meth:`Envelope.check` violation for ``model_id`` — the exact demographic test
    `evaluate_safety` uses to grey a model. This is the stratification axis that turns
    v0.1's *asserted* failure modes ("Schnider misbehaves at high BMI") into *measured*
    ones: compute the metrics separately in each stratum and the envelope claim becomes
    empirical. The split is per-model (each model has its own envelope)."""
    from .covariates import point_patient

    model = ds[model_id]
    in_env: List[SubjectRecord] = []
    out_env: List[SubjectRecord] = []
    for subj in subjects:
        violations = model.applicability_envelope.check(point_patient(subj.covariates))
        (in_env if not violations else out_env).append(subj)
    return in_env, out_env


@dataclass
class LeaderboardEntry:
    """One model's place in the cross-model leaderboard, with its envelope strata."""

    model_id: str
    pd_model: Optional[str]
    declared_tier: str
    overall: CohortValidation
    in_envelope: Optional[CohortValidation] = None
    out_envelope: Optional[CohortValidation] = None

    @property
    def mdape(self) -> float:
        """Overall MDAPE (inaccuracy) — the ranking key. nan sorts last."""
        return self.overall.population.mdape


@dataclass
class Leaderboard:
    """A reproducible, envelope-stratified, cross-model leaderboard on ONE cohort (v0.4 §7.1).

    The point the published numbers cannot make: every model here is scored on the *same*
    subjects, by the *same* kernels, with a pinned seed — so the ranking compares models,
    not cohorts-and-methods. Entries are ordered best (lowest overall MDAPE) first."""

    dataset: str
    target: str
    n_subjects: int
    seed: int
    entries: List[LeaderboardEntry] = field(default_factory=list)
    manifest: Dict[str, Any] = field(default_factory=dict)

    def to_record(self) -> Dict[str, Any]:
        """Serialize the whole leaderboard to a committable, derived-metrics-only artifact."""
        def strat(cv: Optional[CohortValidation]) -> Optional[Dict[str, Any]]:
            return cv.to_record() if cv is not None else None
        return {
            "dataset": self.dataset, "target": self.target,
            "n_subjects": self.n_subjects, "seed": self.seed,
            "manifest": self.manifest,
            "leaderboard": [
                {
                    "rank": i + 1, "model_id": e.model_id, "pd_model": e.pd_model,
                    "declared_tier": e.declared_tier,
                    "overall": strat(e.overall),
                    "in_envelope": strat(e.in_envelope),
                    "out_envelope": strat(e.out_envelope),
                }
                for i, e in enumerate(self.entries)
            ],
        }


def cross_model_leaderboard(
    ds: Dataset,
    subjects: Sequence[SubjectRecord],
    candidates: Sequence[Tuple[str, Optional[str]]],
    *,
    target: str = "bis",
    stratify_by_envelope: bool = True,
    dataset: str = "cohort",
    seed: int = 0,
    manifest: Optional[Dict[str, Any]] = None,
) -> Leaderboard:
    """Score every candidate (PK[, PD]) stack on ONE cohort and rank them (v0.4 §7.1).

    ``candidates`` is a list of ``(pk_model_id, pd_model_id)`` pairs (``pd_model_id`` is
    ``None`` for a ``cp``/``ce`` target). Each is run through the *same*
    :func:`validate_against_cohort` engine on the *same* subjects with the *same* seed,
    so the ranking is apples-to-apples. With ``stratify_by_envelope`` each model is also
    scored on its in- and out-of-envelope sub-cohorts (v0.4 §4.2), turning the asserted
    failure modes into measured ones. A candidate that cannot be scored (no scorable
    observations, kernel pending) is skipped with no fabricated number. Deterministic:
    identical ``(subjects, candidates, seed)`` -> identical leaderboard."""
    entries: List[LeaderboardEntry] = []
    for pk_id, pd_id in candidates:
        try:
            overall = validate_against_cohort(
                ds, pk_id, subjects, target=target, pd_model=pd_id, dataset=dataset, seed=seed)
        except (ValueError, NotImplementedError):
            continue
        if overall.n_subjects == 0:
            continue
        in_cv = out_cv = None
        if stratify_by_envelope:
            in_subj, out_subj = partition_by_envelope(ds, pk_id, subjects)
            if in_subj:
                in_cv = validate_against_cohort(ds, pk_id, in_subj, target=target,
                                                pd_model=pd_id, dataset=dataset, seed=seed)
                in_cv.in_envelope = True
            if out_subj:
                out_cv = validate_against_cohort(ds, pk_id, out_subj, target=target,
                                                 pd_model=pd_id, dataset=dataset, seed=seed)
                out_cv.in_envelope = False
        entries.append(LeaderboardEntry(
            model_id=pk_id, pd_model=pd_id, declared_tier=ds[pk_id].tier,
            overall=overall, in_envelope=in_cv, out_envelope=out_cv))

    # rank by overall MDAPE (inaccuracy), nan last
    entries.sort(key=lambda e: (not np.isfinite(e.mdape), e.mdape))
    return Leaderboard(
        dataset=dataset, target=target, n_subjects=len(list(subjects)), seed=seed,
        entries=entries, manifest=manifest or {})


# Covariate columns recognized by the generic CSV adapter (numeric unless noted).
_CSV_COVARIATES = ("age", "weight", "height", "crcl_ml_min", "albumin_g_dl",
                   "ejection_fraction_pct")


def subjects_from_csv(rows: Sequence[Dict[str, str]]) -> List[SubjectRecord]:
    """Generic long-format cohort CSV -> :class:`SubjectRecord`\\ s (v0.4 §A adapter).

    Each row is **one observation**; rows are grouped by ``subject`` and the
    covariates + dose history are read (constant) from each subject's first row.
    This is the adapter-agnostic, source-neutral entry point — a researcher points
    it at their own observed concentrations (or an Open-TCI / VitalDB export mapped
    to these columns) and feeds the result to :func:`validate_against_cohort`. It
    never invents data: a row missing ``time_min``/``observed``/``kind`` is skipped.

    Recognized columns (others are ignored):

    * ``subject`` — subject id (defaults to ``"all"`` if absent, i.e. one cohort).
    * ``time_min``, ``observed``, ``kind`` — the observation (``kind`` ∈ cp/ce/bis/tof).
    * ``age``/``weight``/``height``/``crcl_ml_min``/``albumin_g_dl``/``ejection_fraction_pct``
      (numeric) and ``sex``/``child_pugh`` (string) — covariates.
    * ``bolus``/``infusion`` — dose specs (e.g. ``"2 mg/kg"`` / ``"6 mg/kg/h"``) applied at t=0.
    """
    grouped: Dict[str, Dict[str, Any]] = {}
    order: List[str] = []
    for row in rows:
        sid = (row.get("subject") or "all").strip()
        try:
            t = float(row["time_min"])
            v = float(row["observed"])
        except (KeyError, TypeError, ValueError):
            continue
        kind = (row.get("kind") or "cp").strip()
        if sid not in grouped:
            order.append(sid)
            cov: Dict[str, Any] = {}
            for c in _CSV_COVARIATES:
                if str(row.get(c, "")).strip():
                    cov[c] = float(row[c])
            for c in ("sex", "child_pugh"):
                if str(row.get(c, "")).strip():
                    cov[c] = row[c].strip()
            schedule: List[Tuple[str, float, str]] = []
            for kind_col in ("bolus", "infusion"):
                spec = str(row.get(kind_col, "")).strip()
                if spec:
                    schedule.append((kind_col, 0.0, spec))
            grouped[sid] = {"covariates": cov, "schedule": schedule, "obs": []}
        grouped[sid]["obs"].append((t, v, kind))
    return [SubjectRecord(covariates=grouped[s]["covariates"], schedule=grouped[s]["schedule"],
                          observations=grouped[s]["obs"], subject_id=s) for s in order]


def subjects_from_vitaldb(
    cases: Sequence[Dict[str, Any]], *, propofol_mg_per_ml: float = 20.0,
) -> List[SubjectRecord]:
    """Map fetched VitalDB cases -> :class:`SubjectRecord`\\ s for PD-BIS validation (v0.4 VE1).

    VitalDB (vitaldb.net, open) records the TCI pump's drug-delivery tracks and the
    *measured* bispectral index. The scientifically independent observation is the
    **measured BIS** (``BIS/BIS``) — not the pump's own predicted Ce, which is just
    another model's output — so this adapter validates a Hypnos PK→BIS stack against
    measured depth-of-anaesthesia. The propofol delivery (``Orchestra/PPF20_RATE``,
    mL/h, PPF20 = 20 mg/mL) is reconstructed as a step-infusion schedule (one event per
    rate change), and the measured BIS becomes the ``kind="bis"`` observations.

    Each ``case`` is a plain dict (the network fetch builds these; this transform is
    pure and offline-testable)::

        {"id": str, "age": float, "sex": "M"/"F", "height": float, "weight": float,
         "infusion_ml_h": [(t_sec, rate_ml_h), ...],   # Orchestra/PPF20_RATE
         "bis":           [(t_sec, bis_value), ...]}    # BIS/BIS

    It **never invents data**: a case with no usable BIS or infusion is dropped, and the
    clinical choices it bakes in (which track is the observation; PPF20 = 20 mg/mL; that
    a propofol-only stack ignores remifentanil's synergistic BIS deepening) are
    DOMAIN-REVIEW items, not verified facts — exactly the curation the community confirms.
    """
    out: List[SubjectRecord] = []
    for case in cases:
        bis = [(float(t) / 60.0, float(v), "bis")
               for t, v in case.get("bis", []) if v is not None and 0.0 < float(v) <= 100.0]
        rate = [(float(t), float(r)) for t, r in case.get("infusion_ml_h", []) if r is not None]
        if not bis or not rate:
            continue
        # reconstruct a step-infusion schedule: one ("infusion", t_min, "<mg/h>") event
        # per change in the delivered rate (mL/h * mg/mL = mg/h), starting from t=0.
        schedule: List[Tuple[str, float, str]] = []
        last: Optional[float] = None
        for t_sec, r_ml_h in rate:
            if last is None or abs(r_ml_h - last) > 1e-9:
                mg_h = r_ml_h * propofol_mg_per_ml
                schedule.append(("infusion", round(float(t_sec) / 60.0, 4), f"{mg_h:g} mg/h"))
                last = r_ml_h
        cov: Dict[str, Any] = {}
        for k in ("age", "height", "weight"):
            if case.get(k) is not None:
                cov[k] = float(case[k])
        if case.get("sex"):
            cov["sex"] = str(case["sex"]).upper()[:1]
        out.append(SubjectRecord(covariates=cov, schedule=schedule, observations=bis,
                                 subject_id=str(case.get("id", f"vdb{len(out) + 1}"))))
    return out


def subjects_from_cohort_self_consistency(
    ds: Dataset, model_id: str, *, target: str = "cp", pd_model: Optional[str] = None,
    offset_pct: float = 20.0, n_subjects: int = 3,
) -> List[SubjectRecord]:
    """A KNOWN-ANSWER cohort built from the model's OWN predictions (v0.4 §6 / §8).

    Observed = predicted·(1 + offset/100), so :func:`validate_against_cohort` must
    recover ``MDPE ≈ offset_pct`` and ``MDAPE ≈ |offset_pct|`` — a CI-runnable
    correctness check on the metric engine that needs **no external data**. This is
    explicitly *not* a clinical cohort; it is the self-consistency fixture the v0.4
    spec validates the engine with before any real-cohort run."""
    from .presets import default_schedule_for
    from .simulate import simulate as _simulate

    schedule = list(default_schedule_for(ds[model_id].drug_name))
    obs_t = np.array([5.0, 10.0, 20.0, 30.0, 45.0, 60.0])
    base = [dict(age=45, weight=70, height=170, sex="M"),
            dict(age=60, weight=85, height=178, sex="M"),
            dict(age=35, weight=60, height=162, sex="F")]
    grid = np.union1d(np.linspace(0.0, float(obs_t.max()), 200), obs_t)
    out: List[SubjectRecord] = []
    for i, cov in enumerate(base[:max(1, n_subjects)]):
        res = _simulate(ds, model_id, patient=cov, schedule=schedule, t=grid, pd_model=pd_model)
        curve = {"cp": res.cp, "ce": res.ce, "bis": res.effect, "tof": res.effect}.get(target)
        if curve is None:
            raise ValueError(f"target={target!r} unavailable (pass pd_model for an effect target)")
        pred = np.interp(obs_t, grid, curve)
        obs = [(float(t), float(p * (1.0 + offset_pct / 100.0)), target)
               for t, p in zip(obs_t, pred) if p > 0]
        out.append(SubjectRecord(covariates=cov, schedule=schedule,
                                 observations=obs, subject_id=f"sc{i + 1}"))
    return out
