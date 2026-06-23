"""Reversal-agent kernels — the v0.5 Part A subsystem.

A reversal agent has no interesting trajectory *alone*; its purpose is to act on a drug
already in the patient. Antagonism is one of three mechanistically distinct *interactions*
(spec §A1), each a small forward-only kernel:

* :func:`encapsulation` — **sugammadex** binds free rocuronium ~1:1 in plasma, removing it
  from the active pool; the free (effect-driving) NMB concentration drops and train-of-four
  recovers. A PK-level coupling: a binding reaction layered onto the NMB's linear PK.
* :func:`competitive_shift` — **naloxone** / **flumazenil** compete at the receptor, shifting
  the agonist's apparent potency (Ce50 → Ce50·(1 + C_antag/Ki)). A PD-level coupling; the
  agonist PK is untouched. The kernel makes **renarcotization** visible: a short-lived
  antagonist washes out before a long-acting agonist, so the effect transiently recovers then
  **relapses** — rendered as honest forward output, never as a re-dosing schedule (§A5).
* :func:`indirect_inhibition` — **neostigmine** raises synaptic ACh to out-compete the NMB,
  with a pronounced **ceiling** at deep block (useless below a depth floor).

There is **no inverse control** (§A6): Hypnos simulates the trajectory of a *given* reversal
dose; it never computes the reversal dose ("how much sugammadex / naloxone?"). That is the
single most acute dosing temptation in the project, and the line is absolute.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .reference import MicroParams

# Molar masses (g/mol) for the molar 1:1 encapsulation bookkeeping. Sugammadex binds
# rocuronium mole-for-mole, so dosing must be compared in moles, not milligrams.
MW_ROCURONIUM = 529.8
MW_SUGAMMADEX = 2178.0


@dataclass
class ReversalTrajectory:
    t: np.ndarray
    free_cp: np.ndarray          # free (active) NMB plasma concentration (ug/mL)
    free_ce: np.ndarray          # free NMB effect-site concentration (drives TOF)
    complex_umol: np.ndarray     # bound (inert) NMB–sugammadex complex (µmol)
    antagonist: Optional[np.ndarray] = None   # antagonist concentration, where relevant


def encapsulation(
    nmb: MicroParams,
    nmb_dose_mg: float,
    sugammadex_dose_mg: float,
    t: np.ndarray,
    *,
    sugammadex_time_min: float = 3.0,
    stoichiometry: float = 1.0,
    k_on: float = 50.0,
    nmb_mw: float = MW_ROCURONIUM,
    sug_mw: float = MW_SUGAMMADEX,
) -> ReversalTrajectory:
    """Forward simulation of sugammadex encapsulating rocuronium (spec §A4).

    Augments the NMB's linear three-compartment PK with a plasma binding reaction that
    depletes the *free* central NMB: ``d/dt = k_on · [NMB_free] · [Sug_free]`` (1:1, the
    bound complex is inert and is not re-released). With a fast ``k_on`` this approximates
    the commonly-modeled instantaneous-1:1 encapsulation — *labeled a simplification* until
    binding kinetics are curated (§A5). The free effect-site concentration (the quantity the
    TOF PD consumes, unchanged) falls after the sugammadex dose, and the block recovers.

    All amounts are tracked in micromoles; doses are converted from mg via the molar masses,
    because the binding is mole-for-mole. Forward-only; computes no dose."""
    from scipy.integrate import solve_ivp

    t = np.asarray(t, dtype=float)
    V1 = nmb.V1
    k10, k12, k21, k13, k31, ke0 = nmb.k10, nmb.k12, nmb.k21, nmb.k13, nmb.k31, nmb.ke0
    nmb_umol = 1e3 * nmb_dose_mg / nmb_mw          # mg / (g/mol) * 1e3 -> µmol
    sug_umol = 1e3 * sugammadex_dose_mg / sug_mw

    def rhs(tt, x):
        A1, A2, A3, S, C, Ce = x
        cr = A1 / V1                                # free central NMB conc (µmol/L)
        cs = max(S, 0.0) / V1                       # free sugammadex conc (plasma)
        bind = k_on * cr * cs * V1                  # µmol/min bound in the central volume
        dA1 = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3 - bind
        dA2 = k12 * A1 - k21 * A2
        dA3 = k13 * A1 - k31 * A3
        dS = -bind / stoichiometry
        dC = bind
        dCe = ke0 * (cr - Ce)
        return [dA1, dA2, dA3, dS, dC, dCe]

    # segment 1: rocuronium bolus at t=0, sugammadex not yet given
    edges = sorted({0.0, float(sugammadex_time_min), float(t.max())})
    x0 = [nmb_umol, 0.0, 0.0, 0.0, 0.0, 0.0]
    ts: list = []
    xs: list = []
    state = x0
    for a, b in zip(edges[:-1], edges[1:]):
        if abs(a - sugammadex_time_min) < 1e-9:
            state = list(state)
            state[3] += sug_umol                    # administer sugammadex (free pool jumps)
        sample = t[(t >= a) & (t <= b)]
        teval = np.unique(np.concatenate([[a, b], sample]))
        sol = solve_ivp(rhs, (a, b), state, t_eval=teval, rtol=1e-7, atol=1e-10, method="LSODA")
        ts.extend(sol.t.tolist())
        xs.extend(sol.y.T.tolist())
        state = sol.y[:, -1]
    ts_arr = np.array(ts)
    xs_arr = np.array(xs)
    # interpolate each state to the requested grid
    out = np.vstack([np.interp(t, ts_arr, xs_arr[:, i]) for i in range(6)]).T
    # convert free central amount (µmol) back to a mass concentration (ug/mL == mg/L):
    # conc_ug_ml = (A1 µmol / V1 L) * MW (g/mol) / 1000
    free_cp = (out[:, 0] / V1) * nmb_mw / 1000.0
    free_ce = (out[:, 5]) * nmb_mw / 1000.0
    return ReversalTrajectory(t=t, free_cp=free_cp, free_ce=free_ce, complex_umol=out[:, 4])


def competitive_shift(antagonist_conc: np.ndarray, ki: float) -> np.ndarray:
    """Apparent-Ce50 multiplier for competitive receptor antagonism (spec §A4).

    Returns ``1 + C_antag / Ki`` elementwise — the factor by which a competitive antagonist
    (naloxone at the µ-opioid receptor, flumazenil at the benzodiazepine site) raises the
    agonist's apparent Ce50, shifting its Hill curve rightward. Fed into
    :func:`hypnos.reference.sigmoid_emax` as ``Ce50·multiplier``; the agonist PK is untouched,
    only its PD potency shifts. Where the antagonist clears faster than the agonist, the
    multiplier falls back toward 1 over time — the mechanistic origin of renarcotization."""
    c = np.clip(np.asarray(antagonist_conc, dtype=float), 0.0, None)
    return 1.0 + c / float(ki)


def indirect_inhibition(
    block_fraction: np.ndarray,
    *,
    ceiling: float = 0.5,
    depth_floor: float = 0.9,
) -> np.ndarray:
    """Neostigmine indirect (anticholinesterase) reversal with its ceiling (spec §A4/§A5).

    ``block_fraction`` is the current neuromuscular block (1 = full block, 0 = recovered).
    Neostigmine raises ACh to out-compete the NMB, but the effect **ceilings**: it cannot
    reverse a block deeper than ``depth_floor`` (the documented deep-block failure — give
    sugammadex instead), and even where it works it removes at most a ``ceiling`` fraction of
    the block. Returns the *residual* block fraction after reversal — a labeled, bounded
    effect, never a dose. The ceiling/depth-floor is the safety content, not a free parameter."""
    b = np.clip(np.asarray(block_fraction, dtype=float), 0.0, 1.0)
    reversible = np.where(b >= depth_floor, 0.0, ceiling)   # no reversal below the depth floor
    return b * (1.0 - reversible)
