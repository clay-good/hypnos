"""Remifentanil Kim 2017 (obesity) kernel + the effect-site divergence fix.

Completes the spec's named remifentanil trio (Minto, Eleveld, Kim). Kim is
PK-only (no published ke0) and uses the Janmahasatian fat-free mass."""
import numpy as np
import pytest

import hypnos
from hypnos.export.registry import instantiate

KIM = "opioids.remifentanil.kim_2017"
MINTO = "opioids.remifentanil.minto_1997"
ELEVELD = "opioids.remifentanil.eleveld_2017"
PAEDFUSOR = "hypnotics_iv.propofol.paedfusor_2005"
ELEVELD_PPF = "hypnotics_iv.propofol.eleveld_2018"
T = np.linspace(0, 30, 181)
REMI = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_reference_individual(ds):
    # weight/age-only parameters are exact at (37 y, 74.5 kg); V2 depends on FFM(height)
    vc = instantiate(ds[KIM], dict(age=37, weight=74.5, height=170, sex="M")).as_volumes_clearances()
    assert abs(vc["V1"] - 4.76) < 1e-6
    assert abs(vc["V3"] - 4.0) < 1e-6
    assert abs(vc["Cl1"] - 2.77) < 1e-6
    assert abs(vc["Cl2"] - 1.94) < 1e-6
    assert abs(vc["Cl3"] - 0.197) < 1e-6
    assert vc["ke0"] == 0.0  # PK-only: no published effect-site link


def test_kim_is_pk_only_no_effect_site(ds):
    res = hypnos.simulate(ds, KIM, patient=dict(age=40, weight=75, height=178, sex="M"),
                          schedule=REMI, t=T)
    assert np.allclose(res.ce, 0.0)  # no ke0 -> no effect-site curve


def test_three_way_remifentanil_compare(ds):
    cmp = hypnos.compare(ds, drug="remifentanil",
                         patient=dict(age=40, weight=75, height=178, sex="M"), schedule=REMI, t=T)
    assert {MINTO, ELEVELD, KIM} <= {r.model_id for r in cmp.included}


def test_kim_covers_morbid_obesity_where_others_greyed(ds):
    cmp = hypnos.compare(ds, drug="remifentanil",
                         patient=dict(age=45, weight=160, height=170, sex="M"), schedule=REMI, t=T)
    assert {r.model_id for r in cmp.included} == {KIM}
    greyed = {e["model_id"] for e in cmp.excluded}
    assert MINTO in greyed and ELEVELD in greyed


def test_ce_divergence_excludes_pk_only_models(ds):
    # the fix: a PK-only model (ke0=0) must not pollute the effect-site spread.
    # For a child, only Eleveld carries an effect compartment among the in-envelope
    # models (Paedfusor is PK-only) -> ce divergence is empty (needs >=2 ke0 models).
    cmp = hypnos.compare(ds, drug="propofol", patient=dict(age=6, weight=20, height=115, sex="M"),
                         schedule=[("bolus", 0.0, "2 mg/kg")], t=T)
    included = {r.model_id for r in cmp.included}
    assert PAEDFUSOR in included and ELEVELD_PPF in included
    assert cmp.divergence["ce"] == {}  # not a spurious huge spread vs Paedfusor's zero ce
    # plasma divergence still spans both
    assert cmp.divergence["cp"] != {}
