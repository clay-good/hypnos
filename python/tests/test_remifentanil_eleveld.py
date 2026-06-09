"""Remifentanil Eleveld 2017 kernel — validated against the published reference
individual; completes the spec's named Minto + Eleveld remifentanil pair."""
import numpy as np
import pytest

import hypnos
from hypnos.export.registry import instantiate

ELEVELD = "opioids.remifentanil.eleveld_2017"
MINTO = "opioids.remifentanil.minto_1997"
T = np.linspace(0, 30, 181)
REMI = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_reference_individual_exact(ds):
    vc = instantiate(ds[ELEVELD], dict(age=35, weight=70, height=170, sex="M")).as_volumes_clearances()
    assert abs(vc["V1"] - 5.81) < 1e-6
    assert abs(vc["V2"] - 8.82) < 1e-6
    assert abs(vc["V3"] - 5.03) < 1e-6   # V3ref (not the tci V2ref typo) -> matches the paper
    assert abs(vc["Cl1"] - 2.58) < 1e-6
    assert abs(vc["Cl2"] - 1.72) < 1e-6
    assert abs(vc["Cl3"] - 0.124) < 1e-6
    assert abs(vc["ke0"] - 1.09) < 1e-6
    assert ds[ELEVELD].review_status == "unverified"


def test_faster_ke0_than_minto(ds):
    # Eleveld remi ke0 (~1.09) is markedly faster than Minto's (~0.6) for a 35y adult
    e = instantiate(ds[ELEVELD], dict(age=35, weight=70, height=170, sex="M")).ke0
    m = instantiate(ds[MINTO], dict(age=35, weight=70, height=170, sex="M")).ke0
    assert e > m


def test_broad_envelope_covers_obese_and_child(ds):
    # Eleveld (FFM, broad envelope) stays in-envelope where Minto (James LBM, adult) is greyed.
    KIM = "opioids.remifentanil.kim_2017"
    # obese adult: Eleveld + Kim (obesity model) both cover; Minto greyed
    obese = hypnos.compare(ds, drug="remifentanil",
                           patient=dict(age=40, weight=140, height=172, sex="M"), schedule=REMI, t=T)
    obese_inc = {r.model_id for r in obese.included}
    assert ELEVELD in obese_inc and KIM in obese_inc
    assert any(e["model_id"] == MINTO for e in obese.excluded)
    # child: only Eleveld covers (Minto and Kim are adult-only)
    child = hypnos.compare(ds, drug="remifentanil",
                           patient=dict(age=6, weight=20, height=115, sex="M"), schedule=REMI, t=T)
    assert {r.model_id for r in child.included} == {ELEVELD}
    greyed = {e["model_id"] for e in child.excluded}
    assert MINTO in greyed and KIM in greyed


def test_adults_cross_validate(ds):
    # for a standard adult the two models should agree reasonably (a sanity cross-check)
    patient = dict(age=40, weight=75, height=178, sex="M")
    cmp = hypnos.compare(ds, drug="remifentanil", patient=patient, schedule=REMI, t=T)
    assert {MINTO, ELEVELD} <= {r.model_id for r in cmp.included}
    assert cmp.divergence["ce"]["max_rel"] < 0.5  # within 50% — closer than the propofol trio
