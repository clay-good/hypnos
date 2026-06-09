"""Eleveld 2018 propofol kernel — validated against the published reference
individual and the source equations (tci R package ``pkmod_eleveld_ppf``)."""
import numpy as np
import pytest

import hypnos
from hypnos.export.registry import instantiate

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
T = np.linspace(0, 60, 361)
SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_reference_individual_exact(ds):
    # 35 y, 70 kg, 170 cm, male, no concomitant anaesthetics, arterial sampling.
    vc = instantiate(ds[ELEVELD], dict(age=35, weight=70, height=170, sex="M")).as_volumes_clearances()
    assert abs(vc["V1"] - 6.28) < 1e-6
    assert abs(vc["V2"] - 25.5) < 1e-6
    assert abs(vc["V3"] - 273.0) < 1e-6
    assert abs(vc["Cl1"] - 1.79) < 1e-6
    assert abs(vc["Cl3"] - 1.11) < 1e-6
    assert abs(vc["ke0"] - 0.146) < 1e-6
    # Q2 carries the (1 + theta16*(1-Q3mat)) term, so the computed reference is ~1.83, not 1.75
    assert abs(vc["Cl2"] - 1.8306) < 1e-3


def test_allometric_clearance_scaling(ds):
    # CL scales ~ (WGT/70)^0.75 for same-age adults (maturation ~ flat in adults)
    light = instantiate(ds[ELEVELD], dict(age=35, weight=50, height=170, sex="M")).as_volumes_clearances()
    heavy = instantiate(ds[ELEVELD], dict(age=35, weight=100, height=170, sex="M")).as_volumes_clearances()
    ratio = heavy["Cl1"] / light["Cl1"]
    assert abs(ratio - (100 / 50) ** 0.75) < 0.02


def test_sex_effect_on_clearance(ds):
    male = instantiate(ds[ELEVELD], dict(age=40, weight=70, height=175, sex="M")).as_volumes_clearances()
    female = instantiate(ds[ELEVELD], dict(age=40, weight=70, height=175, sex="F")).as_volumes_clearances()
    assert female["Cl1"] != male["Cl1"]  # theta4 (male) vs theta15 (female)


def test_opiate_coadministration_lowers_v3(ds):
    base = instantiate(ds[ELEVELD], dict(age=40, weight=70, height=175, sex="M")).as_volumes_clearances()
    op = instantiate(ds[ELEVELD],
                     dict(age=40, weight=70, height=175, sex="M", opiate_coadministration=True)).as_volumes_clearances()
    assert op["V3"] < base["V3"]  # fopiate(theta13<0) shrinks V3


def test_simulate_and_covers_broad_envelope(ds):
    # in-envelope for neonate-to-obese; a child and an obese adult both simulate
    for patient in (dict(age=6, weight=20, height=115, sex="M"),
                    dict(age=40, weight=140, height=172, sex="M")):
        res = hypnos.simulate(ds, ELEVELD, patient=patient, schedule=SCHED, t=T)
        assert res.tier == "A"  # broad envelope -> not tiered down
        assert res.cp_peak > 0 and np.all(np.isfinite(res.cp))


def test_appears_in_three_way_divergence(ds):
    patient = dict(age=72, weight=60, height=162, sex="F")
    cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=SCHED, t=T)
    included = {r.model_id.split(".")[-1] for r in cmp.included}
    assert {"marsh_1991", "schnider_1998", "eleveld_2018"} <= included
