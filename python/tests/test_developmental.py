"""Tests for the v0.8 developmental-pharmacology subsystem.

The two mechanistic kernels (allometric size, PMA maturation) must be verbatim and
composable; an adult model carried into a child must be Tier-D, opt-in, and warned;
and the validator must enforce the five developmental traps — above all that maturation
is driven by PMA, never chronological age.
"""

import numpy as np
import pytest

import hypnos
from hypnos.developmental import (
    allometric_scale,
    apply_developmental,
    linear_per_kg_scale,
    maturation_factor,
)
from hypnos.models import Model
from hypnos.validate import _check_developmental

MARSH = "hypnotics_iv.propofol.marsh_1991"          # allometry_plus_maturation
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"    # allometry_only
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"      # fitted-pediatric (no dev block)


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --- kernels ---------------------------------------------------------------- #
def test_maturation_factor_is_anderson_holford():
    # verbatim sigmoid: MF = PMA^Hill / (TM50^Hill + PMA^Hill)
    pma, tm50, hill = 40.0, 47.7, 3.4
    expect = pma ** hill / (tm50 ** hill + pma ** hill)
    assert maturation_factor(pma, tm50, hill) == pytest.approx(expect)


def test_maturation_monotone_and_bounded():
    vals = [maturation_factor(p, 47.7, 3.4) for p in (28, 40, 60, 120, 600)]
    assert all(0.0 < v <= 1.0 for v in vals)
    assert vals == sorted(vals)                      # increasing in PMA
    assert maturation_factor(47.7, 47.7, 3.4) == pytest.approx(0.5)   # half-max at TM50
    assert maturation_factor(2000, 47.7, 3.4) > 0.99                  # mature at high PMA


def test_allometric_exponents_on_the_right_parameters():
    ref = {"V1": 4.0, "V2": 10.0, "Cl1": 1.8, "Cl2": 1.0, "ke0": 0.46}
    out = allometric_scale(ref, 35.0, cl_exponent=0.75, v_exponent=1.0, reference_weight_kg=70.0)
    r = 35.0 / 70.0
    assert out["V1"] == pytest.approx(4.0 * r ** 1.0)       # volume ∝ weight^1
    assert out["Cl1"] == pytest.approx(1.8 * r ** 0.75)     # clearance ∝ weight^¾
    assert out["ke0"] == pytest.approx(0.46)                # rate constant untouched


def test_allometric_identity_at_reference_weight():
    ref = {"V1": 4.0, "Cl1": 1.8}
    out = allometric_scale(ref, 70.0, reference_weight_kg=70.0)
    assert out == pytest.approx(ref)


def test_linear_per_kg_differs_from_allometry_on_clearance():
    ref = {"V1": 4.0, "Cl1": 1.8}
    w = 3.5
    allo = allometric_scale(ref, w, reference_weight_kg=70.0)
    lin = linear_per_kg_scale(ref, w, reference_weight_kg=70.0)
    # volume exponent is 1.0, so volume coincides; clearance diverges (¾ vs 1)
    assert lin["V1"] == pytest.approx(allo["V1"])
    assert lin["Cl1"] != pytest.approx(allo["Cl1"])
    assert allo["Cl1"] > lin["Cl1"]                 # ¾-power: small bodies clear FASTER per kg


def test_apply_developmental_composes_size_then_maturation(ds):
    m = ds[MARSH]
    dev = m.developmental_model
    ref = {"V1": 4.0, "V2": 10.0, "V3": 100.0, "Cl1": 1.8, "Cl2": 1.0, "Cl3": 0.8, "ke0": 0.0}
    patient = {"weight": 3.4, "pma_weeks": 40}
    out = apply_developmental(ref, dev, patient)
    r = 3.4 / 70.0
    mf = maturation_factor(40, dev.maturation.tm50_weeks, dev.maturation.hill)
    assert out["Cl1"] == pytest.approx(1.8 * r ** 0.75 * mf)   # size then maturation on clearance
    assert out["V1"] == pytest.approx(4.0 * r ** 1.0)          # volume size-only


def test_maturation_requires_pma_not_chronological_age(ds):
    dev = ds[MARSH].developmental_model
    with pytest.raises(ValueError, match="pma_weeks"):
        apply_developmental({"V1": 4.0, "Cl1": 1.8}, dev, {"weight": 3.4})  # no PMA


# --- model layer ------------------------------------------------------------ #
def test_developmental_blocks_curated(ds):
    assert ds[MARSH].developmental_model.extrapolation_basis == "allometry_plus_maturation"
    assert ds[MARSH].developmental_model.has_maturation
    assert ds[SCHNIDER].developmental_model.extrapolation_basis == "allometry_only"
    assert not ds[SCHNIDER].developmental_model.has_maturation
    assert ds[ELEVELD].developmental_model is None     # fitted-pediatric: no extrapolation block
    assert ds[ELEVELD].is_fitted_pediatric             # but it has fitted standing (neonate→elderly)


# --- simulate opt-in -------------------------------------------------------- #
def test_developmental_simulate_forces_tier_d_and_warns(ds):
    neo = dict(pma_weeks=40, weight=3.4, sex="M")
    res = hypnos.simulate(ds, MARSH, patient=neo, schedule=[("bolus", 0.0, "2 mg/kg")],
                          t=np.linspace(0, 30, 200), developmental=True)
    assert res.tier == "D"
    assert any("DEVELOPMENTAL EXTRAPOLATION" in w for w in res.warnings)
    assert any("maturation factor" in w.lower() for w in res.warnings)


def test_allometry_only_surfaces_over_dose_caveat(ds):
    neo = dict(pma_weeks=40, weight=3.4, sex="M")
    res = hypnos.simulate(ds, SCHNIDER, patient=neo, schedule=[("bolus", 0.0, "2 mg/kg")],
                          t=np.linspace(0, 30, 200), developmental=True)
    assert res.tier == "D"
    assert any("allometry_only" in w and "OVER-stated" in w for w in res.warnings)


def test_developmental_never_invents(ds):
    # opt-in on a model with no curated developmental block must refuse, not fabricate
    with pytest.raises(ValueError, match="no developmental_model"):
        hypnos.simulate(ds, ELEVELD, patient=dict(pma_weeks=40, weight=3.4, sex="M"),
                        schedule=[("bolus", 0.0, "2 mg/kg")], t=np.linspace(0, 30, 50),
                        developmental=True)


# --- the headline overlay --------------------------------------------------- #
def test_developmental_overlay_neonate(ds):
    neo = dict(pma_weeks=40, postnatal_age_days=3, weight=3.4, age=0.06, sex="M")
    ov = hypnos.developmental_overlay(ds, drug="propofol", patient=neo,
                                      schedule=[("bolus", 0.0, "2 mg/kg")], t=np.linspace(0, 30, 200))
    # fitted-pediatric models are all greyed at 3.4 kg (no model has standing here)
    assert ov.fitted and all(not f["in_standing"] for f in ov.fitted)
    # the adult models fan out as Tier-D extrapolations, both bases represented
    bases = {e.basis for e in ov.extrapolations}
    assert {"allometry_only", "allometry_plus_maturation"} <= bases
    assert all(e.tier == "D" for e in ov.extrapolations)
    # the maturation extrapolation carries a maturation factor; allometry_only does not
    by_basis = {e.basis: e for e in ov.extrapolations}
    assert by_basis["allometry_plus_maturation"].maturation_factor is not None
    assert by_basis["allometry_only"].maturation_factor is None
    # the visibly-wrong linear-per-kg reference line is present and labeled
    assert ov.linear_per_kg is not None
    assert any("visibly-wrong" in w for w in ov.linear_per_kg.warnings)


# --- validation ------------------------------------------------------------- #
def test_validate_accepts_dataset(ds):
    assert hypnos.validate_dataset(ds) == []


def _dev_model(dev_block):
    return Model({"id": "x.y.z", "drug": {"name": "propofol"}, "purpose": "pk",
                  "tier": "D", "primary_citation": "marsh-1991-propofol-pk",
                  "extraction": {"review_status": "unverified"},
                  "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                  "parameters": [], "developmental_model": dev_block})


def test_validate_rejects_chronological_age_maturation_clock():
    bad = _dev_model({
        "extrapolation_basis": "allometry_plus_maturation", "evidence_tier": "D",
        "applied_by_default": False, "primary_citation": "anderson-holford-2008",
        "maturation": {"function": "sigmoidal_pma", "tm50_weeks": 47.7, "hill": 3.4,
                       "driver": "age_years", "affected_parameter": "Cl"}})
    problems = _check_developmental(bad, {"anderson-holford-2008"})
    assert any("PMA-class clock" in p for p in problems)


def test_validate_rejects_default_applied_extrapolation():
    bad = _dev_model({"extrapolation_basis": "allometry_only", "evidence_tier": "D",
                      "applied_by_default": True, "primary_citation": "anderson-holford-2008"})
    assert any("applied_by_default" in p for p in _check_developmental(bad, {"anderson-holford-2008"}))


def test_validate_rejects_non_d_extrapolation_tier():
    bad = _dev_model({"extrapolation_basis": "allometry_only", "evidence_tier": "B",
                      "applied_by_default": False, "primary_citation": "anderson-holford-2008"})
    assert any("must be D" in p for p in _check_developmental(bad, {"anderson-holford-2008"}))
