"""Pure-NumPy/SciPy reference kernels.

Every Hypnos model binds to one of these kernels so exports can be round-trip
validated against it. The kernels implement exactly what is needed and nothing
more:

* an n-compartment mammillary linear-ODE PK solver (exact matrix exponential);
* the effect-compartment first-order link (the *k*e0 model);
* the sigmoid E_max / Hill PD transform.

Hypnos is **not** a TCI engine: these solve the *forward* problem (dose history
-> predicted concentration/effect). There is no inverse-control / target-seeking
code here, by design (see spec §10).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import expm


# --------------------------------------------------------------------------- #
# Covariate helpers
# --------------------------------------------------------------------------- #
def lbm_james(weight: float, height_cm: float, sex: str) -> float:
    """James (1976) lean body mass. ``sex`` in {'M','F'} (case-insensitive)."""
    ratio = (weight / height_cm) ** 2
    if str(sex).upper().startswith("M"):
        return 1.1 * weight - 128.0 * ratio
    return 1.07 * weight - 148.0 * ratio


def ffm_janmahasatian(weight: float, height_cm: float, sex: str) -> float:
    """Janmahasatian (2005) fat-free mass (kg). BMI-driven, sex-specific.

    The form Kim 2017's remifentanil disposition model scales on (and the size
    kernel of the Eleveld models is built from). ``sex`` in {'M','F'}.
    """
    b = bmi(weight, height_cm)
    if str(sex).upper().startswith("M"):
        return 9.27e3 * weight / (6.68e3 + 216.0 * b)
    return 9.27e3 * weight / (8.78e3 + 244.0 * b)


def ffm_al_sallami(weight: float, height_cm: float, age: float, sex: str) -> float:
    """Al-Sallami (2015) fat-free mass (kg) — the age-aware FFM the Eleveld models use.

    An age-dependent multiplier on the Janmahasatian (2005) FFM that extends validity
    into childhood (where Janmahasatian was not derived). ``sex`` in {'M','F'}.
    """
    jan = ffm_janmahasatian(weight, height_cm, sex)
    if str(sex).upper().startswith("M"):
        return (0.88 + (1.0 - 0.88) / (1.0 + (age / 13.4) ** (-12.7))) * jan
    return (1.11 + (1.0 - 1.11) / (1.0 + (age / 7.1) ** (-1.1))) * jan


def bmi(weight: float, height_cm: float) -> float:
    return weight / ((height_cm / 100.0) ** 2)


# --------------------------------------------------------------------------- #
# Parameter container
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MicroParams:
    """Micro-rate-constant parameterization of a (≤3)-compartment model + ke0.

    Any volumes/clearances model is reduced to this canonical form so the
    solver and exporters share a single representation.
    """

    V1: float
    k10: float
    k12: float = 0.0
    k21: float = 0.0
    k13: float = 0.0
    k31: float = 0.0
    ke0: float = 0.0
    n_compartments: int = 3
    derived: Dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_volumes_clearances(
        cls,
        V1: float,
        Cl1: float,
        V2: float = 0.0,
        Cl2: float = 0.0,
        V3: float = 0.0,
        Cl3: float = 0.0,
        ke0: float = 0.0,
    ) -> "MicroParams":
        k10 = Cl1 / V1
        k12 = Cl2 / V1 if V2 else 0.0
        k21 = Cl2 / V2 if V2 else 0.0
        k13 = Cl3 / V1 if V3 else 0.0
        k31 = Cl3 / V3 if V3 else 0.0
        n = 1 + (1 if V2 else 0) + (1 if V3 else 0)
        return cls(
            V1=V1, k10=k10, k12=k12, k21=k21, k13=k13, k31=k31, ke0=ke0,
            n_compartments=n,
            derived={"V2": V2, "V3": V3, "Cl1": Cl1, "Cl2": Cl2, "Cl3": Cl3},
        )

    def as_volumes_clearances(self) -> Dict[str, float]:
        """Present the model in volumes/clearances form (for export/reporting)."""
        V1 = self.V1
        Cl1 = self.k10 * V1
        V2 = (self.k12 * V1 / self.k21) if self.k21 else 0.0
        Cl2 = self.k12 * V1
        V3 = (self.k13 * V1 / self.k31) if self.k31 else 0.0
        Cl3 = self.k13 * V1
        return {"V1": V1, "V2": V2, "V3": V3, "Cl1": Cl1, "Cl2": Cl2, "Cl3": Cl3, "ke0": self.ke0}

    def system_matrix(self) -> np.ndarray:
        """4x4 state matrix for [A1, A2, A3, Ce]. Amounts in mg; Ce a concentration."""
        k10, k12, k21, k13, k31, ke0 = (
            self.k10, self.k12, self.k21, self.k13, self.k31, self.ke0
        )
        M = np.array(
            [
                [-(k10 + k12 + k13), k21, k31, 0.0],
                [k12, -k21, 0.0, 0.0],
                [k13, 0.0, -k31, 0.0],
                [ke0 / self.V1, 0.0, 0.0, -ke0],
            ],
            dtype=float,
        )
        return M


# --------------------------------------------------------------------------- #
# Dosing representation
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Dosing:
    """Normalized dosing: instantaneous boluses (mg) and step-changing infusions."""

    boluses: Tuple[Tuple[float, float], ...] = ()           # (time_min, amount_mg)
    infusions: Tuple[Tuple[float, float], ...] = ()         # (start_min, rate_mg_per_min)

    def active_rate(self, t: float) -> float:
        """Step-function infusion rate: the most recent infusion event at/<= t."""
        rate = 0.0
        latest = -math.inf
        for start, r in self.infusions:
            if start <= t and start >= latest:
                latest = start
                rate = r
        return rate


# --------------------------------------------------------------------------- #
# Trajectory result
# --------------------------------------------------------------------------- #
@dataclass
class Trajectory:
    t: np.ndarray
    cp: np.ndarray          # plasma concentration (ug/mL == mg/L)
    ce: np.ndarray          # effect-site concentration (ug/mL)
    amounts: np.ndarray     # (len(t), 4) state amounts/conc
    units: str = "ug/mL"


# --------------------------------------------------------------------------- #
# Exact solver (matrix exponential)
# --------------------------------------------------------------------------- #
def simulate(params: MicroParams, dosing: Dosing, t: np.ndarray) -> Trajectory:
    """Exact forward simulation of Cp(t) and Ce(t) via the matrix exponential.

    The system is linear time-invariant with piecewise-constant infusion and
    instantaneous boluses, so each segment has a closed form:

        x(t0+dt) = expm(M·dt)·x(t0) + M⁻¹·(expm(M·dt) − I)·b

    where ``b`` injects the constant infusion rate into the central compartment.
    """
    t = np.asarray(t, dtype=float)
    M = params.system_matrix()
    n = M.shape[0]

    # Breakpoints: every place the dynamics or sampling changes.
    bps = set(t.tolist())
    bps.add(0.0)
    for bt, _ in dosing.boluses:
        bps.add(bt)
    for st, _ in dosing.infusions:
        bps.add(st)
    ordered = sorted(b for b in bps if b >= t.min() or b <= 0.0)
    # ensure we start at the earliest of 0 and the first sample
    start = min(0.0, float(t.min()))
    ordered = sorted(set([start] + [b for b in ordered if b >= start]))

    # Cache state at each breakpoint.
    state = np.zeros(n)
    # Memoize the augmented matrix exponential by (dt, rate): on a uniform time
    # grid with piecewise-constant infusion the same propagator recurs at almost
    # every step, so this turns O(n) expm calls into ~1 (identical results).
    phi_cache: Dict[Any, np.ndarray] = {}
    state_at: Dict[float, np.ndarray] = {}
    for i, bp in enumerate(ordered):
        # apply boluses at this instant
        for bt, amt in dosing.boluses:
            if bt == bp:
                state = state.copy()
                state[0] += amt
        state_at[bp] = state.copy()
        if i + 1 < len(ordered):
            dt = ordered[i + 1] - bp
            if dt > 0:
                rate = dosing.active_rate(bp + dt / 2.0)
                key = (round(dt, 12), rate)
                phi_aug = phi_cache.get(key)
                if phi_aug is None:
                    # Augmented matrix exponential solves x' = M x + b exactly for
                    # constant b, even when M is singular (e.g. ke0 = 0):
                    #   [x'; 0] = [[M, b], [0, 0]] [x; 1]
                    aug = np.zeros((n + 1, n + 1))
                    aug[:n, :n] = M
                    aug[0, n] = rate  # infusion enters the central compartment
                    phi_aug = expm(aug * dt)
                    phi_cache[key] = phi_aug
                state = phi_aug[:n, :n] @ state + phi_aug[:n, n]

    # Sample at requested times (each t is a breakpoint, so look it up).
    out = np.zeros((len(t), n))
    for j, tj in enumerate(t):
        out[j] = state_at[float(tj)]

    cp = out[:, 0] / params.V1
    ce = out[:, 3]
    return Trajectory(t=t, cp=cp, ce=ce, amounts=out)


# --------------------------------------------------------------------------- #
# Population variability — virtual-individual sampling (v0.2 spec §6)
# --------------------------------------------------------------------------- #
# Symbols carried in the volumes/clearances view a curated BSV block may key on.
_VC_SYMBOLS = ("V1", "V2", "V3", "Cl1", "Cl2", "Cl3", "ke0")


def sample_individual(typical: MicroParams, omegas: Dict[str, float], rng) -> MicroParams:
    """Draw one *virtual individual* by perturbing the typical-value parameters.

    Exponential (log-normal) between-subject variability, the near-universal
    convention for disposition parameters:

        P_i = P_typical · exp(η_i),   η_i ~ N(0, ω²)

    ``omegas`` maps a volumes/clearances symbol (V1, V2, V3, Cl1, Cl2, Cl3, ke0) to
    its η-scale variance ω² (``omega2``). Parameters absent from ``omegas`` are held
    at the typical value (no fabricated variability — spec §5). Draws are taken from
    the caller-supplied seeded ``rng`` so a band is byte-reproducible (spec §6).

    Off-diagonal Ω (correlated η's) is curated separately and not yet sampled here;
    each η is drawn independently, which is the honest default a ``diagonal``
    ``variability_status`` records (spec §3.3 / §5).
    """
    vc = typical.as_volumes_clearances()
    for sym in _VC_SYMBOLS:
        om2 = omegas.get(sym)
        if om2 and om2 > 0 and vc.get(sym, 0.0) > 0:
            vc[sym] = vc[sym] * math.exp(rng.normal(0.0, math.sqrt(om2)))
    return MicroParams.from_volumes_clearances(
        V1=vc["V1"], Cl1=vc["Cl1"], V2=vc["V2"], Cl2=vc["Cl2"],
        V3=vc["V3"], Cl3=vc["Cl3"], ke0=vc["ke0"],
    )


def sample_parameter_vector(typical: MicroParams, ses: Dict[str, float], rng) -> MicroParams:
    """Draw one *typical-value* parameter vector perturbed by its ESTIMATION uncertainty (v0.3 E1).

    This is the confidence-band analog of :func:`sample_individual`. Where ``sample_individual``
    scatters a *virtual individual* around the typical patient by the between-subject Ω
    (irreducible population spread), this scatters the *typical curve itself* by the SE on each
    estimated θ — how well the fitting data pinned the population mean down (reducible with more
    data). The two are conceptually distinct and never conflated (v0.3 §1).

        θ_b = θ_hat + N(0, SE_θ)     (natural scale; floored at a small positive)

    ``ses`` maps a volumes/clearances symbol (V1/V2/V3/Cl1/Cl2/Cl3/ke0) to its natural-scale SE.
    Parameters absent from ``ses`` are held at the point estimate (no fabricated uncertainty —
    the never-synthesize rule, v0.3 §5). Draws come from the seeded ``rng`` so a confidence band
    is byte-reproducible."""
    vc = typical.as_volumes_clearances()
    for sym in _VC_SYMBOLS:
        se = ses.get(sym)
        if se and se > 0 and vc.get(sym, 0.0) > 0:
            vc[sym] = max(vc[sym] + rng.normal(0.0, se), 1e-9)
    return MicroParams.from_volumes_clearances(
        V1=vc["V1"], Cl1=vc["Cl1"], V2=vc["V2"], Cl2=vc["Cl2"],
        V3=vc["V3"], Cl3=vc["Cl3"], ke0=vc["ke0"],
    )


def residual_std(
    conc: np.ndarray,
    model: str,
    *,
    prop_var: Optional[float] = None,
    add_sd: Optional[float] = None,
    log_sd: Optional[float] = None,
) -> np.ndarray:
    """Pointwise residual standard deviation of an *observation* given the predicted
    concentration ``conc`` and a residual-error model (spec §2.2).

    * ``proportional`` — SD = conc · sqrt(prop_var)
    * ``additive``     — SD = add_sd (constant)
    * ``combined``     — SD = sqrt((conc·sqrt(prop_var))² + add_sd²)
    * ``log``          — additive on log scale; SD ≈ conc · sqrt(exp(log_sd²) − 1)
    """
    conc = np.asarray(conc, dtype=float)
    if model == "proportional":
        return conc * math.sqrt(prop_var or 0.0)
    if model == "additive":
        return np.full_like(conc, float(add_sd or 0.0))
    if model == "combined":
        return np.sqrt((conc * math.sqrt(prop_var or 0.0)) ** 2 + float(add_sd or 0.0) ** 2)
    if model == "log":
        return conc * math.sqrt(math.exp(float(log_sd or 0.0) ** 2) - 1.0)
    raise ValueError(f"unknown residual-error model {model!r}")


def apply_residual(
    conc: np.ndarray,
    model: str,
    rng,
    *,
    prop_var: Optional[float] = None,
    add_sd: Optional[float] = None,
    log_sd: Optional[float] = None,
) -> np.ndarray:
    """Draw observation-level samples by applying residual error to a trajectory."""
    conc = np.asarray(conc, dtype=float)
    if model == "log":
        return conc * np.exp(rng.normal(0.0, float(log_sd or 0.0), size=conc.shape))
    if model == "additive":
        return conc + rng.normal(0.0, float(add_sd or 0.0), size=conc.shape)
    if model == "proportional":
        return conc * (1.0 + rng.normal(0.0, math.sqrt(prop_var or 0.0), size=conc.shape))
    if model == "combined":
        return conc * (1.0 + rng.normal(0.0, math.sqrt(prop_var or 0.0), size=conc.shape)) \
            + rng.normal(0.0, float(add_sd or 0.0), size=conc.shape)
    raise ValueError(f"unknown residual-error model {model!r}")


# --------------------------------------------------------------------------- #
# Independent numerical solver (for cross-validation only)
# --------------------------------------------------------------------------- #
def simulate_numeric(params: MicroParams, dosing: Dosing, t: np.ndarray) -> Trajectory:
    """Independent scipy ``solve_ivp`` integration, used to cross-check :func:`simulate`."""
    from scipy.integrate import solve_ivp

    t = np.asarray(t, dtype=float)
    M = params.system_matrix()
    n = M.shape[0]

    # event breakpoints for boluses
    bolus_times = sorted({bt for bt, _ in dosing.boluses})

    def rhs(tt, x):
        b = np.zeros(n)
        b[0] = dosing.active_rate(tt)
        return M @ x + b

    # integrate piecewise so boluses are applied as jumps
    seg_edges = sorted(set([min(0.0, float(t.min())), float(t.max())] + bolus_times
                           + [st for st, _ in dosing.infusions]))
    seg_edges = [e for e in seg_edges if e <= t.max()]
    state = np.zeros(n)
    # storage
    ts_all: List[float] = []
    xs_all: List[np.ndarray] = []
    for i in range(len(seg_edges) - 1):
        a, b_ = seg_edges[i], seg_edges[i + 1]
        for bt, amt in dosing.boluses:
            if bt == a:
                state = state.copy()
                state[0] += amt
        if b_ <= a:
            continue
        sample = t[(t >= a) & (t <= b_)]
        sol = solve_ivp(rhs, (a, b_), state, t_eval=np.unique(np.concatenate([[a, b_], sample])),
                        rtol=1e-9, atol=1e-12, method="LSODA")
        for k, tk in enumerate(sol.t):
            ts_all.append(tk)
            xs_all.append(sol.y[:, k])
        state = sol.y[:, -1]
    # apply a bolus exactly at t.max() edge if present
    # interpolate to requested t
    ts_arr = np.array(ts_all)
    xs_arr = np.array(xs_all)
    out = np.zeros((len(t), n))
    for ci in range(n):
        out[:, ci] = np.interp(t, ts_arr, xs_arr[:, ci])
    cp = out[:, 0] / params.V1
    ce = out[:, 3]
    return Trajectory(t=t, cp=cp, ce=ce, amounts=out)


# --------------------------------------------------------------------------- #
# PD transform
# --------------------------------------------------------------------------- #
def sigmoid_emax(
    ce: np.ndarray, E0: float, Emax: float, Ce50: float, gamma: float
) -> np.ndarray:
    """Sigmoid E_max (Hill) effect: E = E0 − Emax·Ce^γ / (Ce50^γ + Ce^γ)."""
    ce = np.asarray(ce, dtype=float)
    num = np.power(np.clip(ce, 0.0, None), gamma)
    den = Ce50 ** gamma + num
    return E0 - Emax * np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def mac_age_corrected(mac40: float, age: float) -> float:
    """Age-corrected MAC via the Mapleson/Nickalls relation (BJA 2003).

        MAC(age) = MAC40 · 10^(−0.00269·(age − 40))

    Universal across inhaled anaesthetics (~6% decrease per decade), anchored at
    age 40. Valid for age > 1 year.
    """
    return mac40 * 10.0 ** (-0.00269 * (age - 40.0))


def mac_fraction(end_tidal_pct: float, mac40: float, age: float) -> float:
    """MAC fraction (depth surrogate) = end-tidal concentration / age-corrected MAC."""
    return end_tidal_pct / mac_age_corrected(mac40, age)


def alveolar_washin(
    blood_gas: float,
    t: np.ndarray,
    *,
    alveolar_ventilation: float = 4.0,
    frc: float = 2.5,
    cardiac_output: float = 5.0,
):
    """Single-compartment alveolar wash-in of an inhaled agent: FA/FI(t).

    Mass balance on the lung (the functional residual capacity, FRC), with the
    inspired fraction FI held constant and the early, pre-recirculation assumption
    that the mixed-venous partial pressure is still ≈ 0 (so uptake into blood is
    ``λ·Q̇·FA``, where ``λ`` is the blood:gas partition coefficient and ``Q̇`` is
    cardiac output):

        FRC · dFA/dt = V̇_A·(FI − FA) − λ·Q̇·FA

    which integrates (FA(0)=0) to

        FA/FI(t) = plateau · (1 − e^(−t/τ)),
        plateau  = V̇_A / (V̇_A + λ·Q̇),     τ = FRC / (V̇_A + λ·Q̇)

    The discriminating, solubility-driven quantity is the **early plateau** — the
    FA/FI "knee" reached in the first few minutes before tissue uptake dominates.
    A less soluble agent (smaller ``λ``) has a *higher* plateau and rises toward it
    faster, i.e. washes in faster — the canonical ordering (desflurane / nitrous
    oxide fast, isoflurane slow). This is a comparative, education-grade kernel
    with stated standard 70-kg-adult ventilation constants (V̇_A ≈ 4 L/min, FRC ≈
    2.5 L, Q̇ ≈ 5 L/min, all overridable); it is **not** a per-patient FA/FI
    predictor, and full multi-compartment tissue uptake (Mapleson) is out of scope.

    Returns ``(fa_fi, plateau, tau)`` — the FA/FI array over ``t`` (minutes), the
    early plateau (dimensionless), and the time constant ``τ`` (minutes).
    """
    t = np.asarray(t, dtype=float)
    denom = alveolar_ventilation + blood_gas * cardiac_output
    plateau = alveolar_ventilation / denom
    tau = frc / denom
    fa_fi = plateau * (1.0 - np.exp(-t / tau))
    return fa_fi, plateau, tau


def alveolar_washout(
    blood_gas: float,
    t: np.ndarray,
    *,
    alveolar_ventilation: float = 4.0,
    frc: float = 2.5,
    cardiac_output: float = 5.0,
):
    """Single-compartment alveolar wash-out (emergence) of an inhaled agent: FA/FA₀(t).

    The offset mirror of :func:`alveolar_washin`. The inspired fraction is now zero
    (``FI = 0``) and the agent leaves the alveoli through ventilation; the symmetric
    early-phase idealization is that the mixed-venous fraction is still ≈ the initial
    alveolar fraction ``FA₀`` (tissues are saturated and have not yet given up their
    load), so blood *redelivers* agent at ``λ·Q̇·FA₀``:

        FRC · dFA/dt = −V̇_A·FA + λ·Q̇·(FA₀ − FA)

    which integrates (FA(0)=FA₀) to

        FA/FA₀(t) = floor + (1 − floor)·e^(−t/τ),
        floor     = λ·Q̇ / (V̇_A + λ·Q̇) = 1 − plateau,   τ = FRC / (V̇_A + λ·Q̇)

    The discriminating, solubility-driven quantity is the **early elimination floor**
    ``floor = 1 − plateau`` — the residual alveolar fraction the single-compartment
    model settles toward before tissue release dominates. A *less* soluble agent
    (smaller ``λ``) has a *lower* floor and so washes out more completely and faster —
    the canonical ordering (desflurane / nitrous oxide fast emergence, isoflurane
    slow). This is the exact symmetry of the wash-in plateau, and it reproduces the
    clinical reason desflurane is chosen for long cases.

    Same honesty boundary as wash-in: a comparative, education-grade kernel with
    stated standard 70-kg-adult ventilation constants, all overridable; it is **not**
    a per-patient emergence predictor, and the full multi-compartment tissue release
    (Mapleson) that governs the long offset tail is out of scope — this captures only
    the early, lung-dominated phase.

    Returns ``(fa_fa0, floor, tau)`` — the FA/FA₀ array over ``t`` (minutes), the
    early elimination floor (dimensionless), and the time constant ``τ`` (minutes).
    """
    t = np.asarray(t, dtype=float)
    denom = alveolar_ventilation + blood_gas * cardiac_output
    floor = blood_gas * cardiac_output / denom
    tau = frc / denom
    fa_fa0 = floor + (1.0 - floor) * np.exp(-t / tau)
    return fa_fa0, floor, tau


def sigmoid_emax_twoslope(
    ce: np.ndarray, E0: float, Emax: float, Ce50: float, gamma_low: float, gamma_high: float
) -> np.ndarray:
    """Two-slope sigmoid E_max (the Eleveld BIS form): a Hill curve whose steepness
    differs below vs above Ce50.

        γ = γ_low  if Ce ≤ Ce50  else  γ_high
        E = E0 − Emax · Ce^γ / (Ce50^γ + Ce^γ)

    Continuous at Ce = Ce50 (both branches give the half-maximal effect there),
    but asymmetric: propofol's concentration–BIS relationship is shallower below
    the half-effect point and steeper above it.
    """
    ce = np.clip(np.asarray(ce, dtype=float), 0.0, None)
    gam = np.where(ce <= Ce50, gamma_low, gamma_high)
    num = np.power(ce, gam)
    den = np.power(Ce50, gam) + num
    return E0 - Emax * np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def greco_response_surface(
    ce_a: np.ndarray,
    ce_b: np.ndarray,
    E0: float,
    Emax: float,
    Ce50_a: float,
    Ce50_b: float,
    alpha: float,
    gamma: float,
) -> np.ndarray:
    """Greco-type two-drug response surface (the canonical interaction model).

    With normalized concentrations ``uA = ce_a/Ce50_a`` and ``uB = ce_b/Ce50_b``::

        U = uA + uB + alpha * uA * uB
        E = E0 − Emax · U^γ / (1 + U^γ)

    ``alpha`` is the interaction term: ``alpha > 0`` is **synergy** (the combination
    is more potent than additivity), ``alpha = 0`` is exact additivity, ``alpha < 0``
    is antagonism. At ``ce_b = 0`` this collapses to the single-drug sigmoid in
    drug A, so an interaction surface composes consistently with the PK kernels.

    Concentrations and Ce50s must share one unit system (Hypnos uses µg/mL = mg/L
    throughout; for remifentanil that is ng/mL ÷ 1000).
    """
    ua = np.clip(np.asarray(ce_a, dtype=float), 0.0, None) / Ce50_a
    ub = np.clip(np.asarray(ce_b, dtype=float), 0.0, None) / Ce50_b
    U = ua + ub + alpha * ua * ub
    num = np.power(U, gamma)
    den = 1.0 + num
    return E0 - Emax * np.divide(num, den, out=np.zeros_like(num), where=den > 0)
