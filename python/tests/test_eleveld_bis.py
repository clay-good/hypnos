"""Eleveld two-slope BIS PD model — validated against the published form
(tci ``emax_eleveld``): BIS = E0 - Emax*Ce^g/(Ce50^g+Ce^g), g split at Ce50."""
import numpy as np
import pytest

import hypnos
from hypnos.reference import sigmoid_emax_twoslope

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
ELEVELD_BIS = "pd_effect.propofol.eleveld_bis"
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
T = np.linspace(0, 60, 361)
SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_twoslope_endpoints_and_halfeffect():
    b = sigmoid_emax_twoslope(np.array([0.0, 3.08, 1e6]), 93.0, 93.0, 3.08, 1.47, 1.89)
    assert abs(b[0] - 93.0) < 1e-9          # no drug -> baseline
    assert abs(b[1] - 46.5) < 1e-9          # at Ce50 -> half-maximal effect exactly
    assert b[2] < 0.5                        # saturating -> approaches 0


def test_twoslope_is_continuous_at_ce50():
    eps = 1e-6
    lo = sigmoid_emax_twoslope(np.array([3.08 - eps]), 93, 93, 3.08, 1.47, 1.89)[0]
    hi = sigmoid_emax_twoslope(np.array([3.08 + eps]), 93, 93, 3.08, 1.47, 1.89)[0]
    assert abs(lo - hi) < 1e-3              # continuous across the slope switch


def test_steeper_above_ce50():
    # with gamma_high > gamma_low the curve falls faster just above Ce50 than it
    # rises (in BIS-drop terms) just below it
    drop_below = 46.5 - sigmoid_emax_twoslope(np.array([3.08 * 1.2]), 93, 93, 3.08, 1.47, 1.89)[0]
    drop_below_low_gamma = 46.5 - sigmoid_emax_twoslope(np.array([3.08 * 1.2]), 93, 93, 3.08, 1.47, 1.47)[0]
    assert drop_below > drop_below_low_gamma  # the steeper high-side gamma deepens effect faster


def test_eleveld_pk_pd_composition(ds):
    patient = dict(age=72, weight=60, height=162, sex="F")
    res = hypnos.simulate(ds, ELEVELD, patient=patient, schedule=SCHED, t=T, pd_model=ELEVELD_BIS)
    assert res.effect is not None
    assert 0.0 <= res.effect.min() and res.effect.max() <= 93.0 + 1e-9
    assert res.tier == "B"  # worst of Eleveld PK (A) and Eleveld BIS (B)


def test_age_lowers_ce50_deepens_bis(ds):
    # Eleveld Ce50 = 3.08*exp(-0.00635*(age-35)) -> older patients are more sensitive
    young = hypnos.simulate(ds, ELEVELD, patient=dict(age=25, weight=70, height=175, sex="M"),
                            schedule=SCHED, t=T, pd_model=ELEVELD_BIS)
    old = hypnos.simulate(ds, ELEVELD, patient=dict(age=85, weight=70, height=175, sex="M"),
                          schedule=SCHED, t=T, pd_model=ELEVELD_BIS)
    assert old.effect.min() < young.effect.min()


def test_single_slope_pd_still_works(ds):
    # regression: the original single-gamma BIS sigmoid is unaffected by the dispatch
    res = hypnos.simulate(ds, SCHNIDER, patient=dict(age=50, weight=77, height=177, sex="M"),
                          schedule=SCHED, t=T, pd_model="pd_effect.propofol.bis_sigmoid")
    assert res.effect is not None and res.tier == "C"
