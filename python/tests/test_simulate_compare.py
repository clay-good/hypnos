import numpy as np
import pytest

import hypnos

T = np.linspace(0, 60, 361)
SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
MARSH = "hypnotics_iv.propofol.marsh_1991"
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
BIS = "pd_effect.propofol.bis_sigmoid"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_bolus_sets_central_concentration(ds):
    # 2 mg/kg into a 70 kg patient = 140 mg; Schnider V1 = 4.27 L -> Cp(0)=140/4.27.
    patient = dict(age=50, weight=70, height=175, sex="M")
    res = hypnos.simulate(ds, SCHNIDER, patient=patient, schedule=[("bolus", 0.0, "2 mg/kg")], t=T)
    assert abs(res.cp[0] - 140.0 / 4.27) < 1e-6


def test_within_envelope_no_warning(ds):
    patient = dict(age=50, weight=77, height=177, sex="M")
    res = hypnos.simulate(ds, SCHNIDER, patient=patient, schedule=SCHED, t=T)
    assert res.tier == "B"
    assert res.warnings == []
    assert not res.excluded


def test_out_of_envelope_tiers_down_to_D(ds):
    obese = dict(age=40, weight=140, height=172, sex="M")
    res = hypnos.simulate(ds, SCHNIDER, patient=obese, schedule=SCHED, t=T)
    assert res.tier == "D"
    assert res.excluded
    assert any("ENVELOPE" in w for w in res.warnings)
    assert any("FAILURE MODE" in w for w in res.warnings)


def test_eleveld_kernel_reproduces_reference(ds):
    # Eleveld now has an executable kernel; it must reproduce the published
    # reference individual (35 y, 70 kg, 170 cm, male): V1=6.28, CL=1.79, ke0=0.146.
    from hypnos.export.registry import instantiate

    p = instantiate(ds[ELEVELD], dict(age=35, weight=70, height=170, sex="M"))
    vc = p.as_volumes_clearances()
    assert abs(vc["V1"] - 6.28) < 1e-6
    assert abs(vc["V2"] - 25.5) < 1e-6
    assert abs(vc["V3"] - 273.0) < 1e-6
    assert abs(vc["Cl1"] - 1.79) < 1e-6
    assert abs(vc["Cl3"] - 1.11) < 1e-6
    assert abs(vc["ke0"] - 0.146) < 1e-6
    # the record stays unverified: an LLM transcription is not a human PDF check
    assert ds[ELEVELD].review_status == "unverified"


def test_pd_tier_propagation_worst_wins(ds):
    # PK is Tier B, PD is Tier C -> composed simulation is Tier C.
    patient = dict(age=50, weight=77, height=177, sex="M")
    res = hypnos.simulate(ds, SCHNIDER, patient=patient, schedule=SCHED, t=T, pd_model=BIS)
    assert res.tier == "C"
    assert res.effect is not None
    assert res.effect.min() >= 0.0
    assert res.effect.max() <= 93.0 + 1e-9


def test_compare_groups_and_divergence(ds):
    patient = dict(age=72, weight=60, height=162, sex="F")
    cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=SCHED, t=T)
    ids = {r.model_id for r in cmp.included}
    # all three adult propofol models now have kernels and are in-envelope here
    assert {MARSH, SCHNIDER, ELEVELD} <= ids
    assert cmp.divergence["ce"]["max_abs"] > 0.5  # the models genuinely disagree


def test_divergence_names_the_driver_pair(ds):
    # the spread says how much; the driver says which model is the outlier.
    patient = dict(age=72, weight=60, height=162, sex="F")
    cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=SCHED, t=T)
    drv = cmp.divergence["ce"]["driver"]
    # Schnider's small fixed V1 + fast ke0 make it the effect-site outlier here
    assert SCHNIDER in (drv["high"], drv["low"])
    assert drv["high"] != drv["low"]
    # the named gap equals the reported peak absolute spread
    assert abs(drv["gap"] - cmp.divergence["ce"]["max_abs"]) < 1e-9


def test_plasma_divergence_spans_pk_only_models(ds):
    # a pediatric compare: Kataria & Paedfusor are PK-only (no ce), so the effect-site
    # spread can't see them — but the plasma spread must, and name the pair.
    child = dict(age=6, weight=20, height=115, sex="M")
    cmp = hypnos.compare(ds, drug="propofol", patient=child, schedule=SCHED, t=T)
    cp = cmp.divergence["cp"]
    assert cp and cp["max_abs"] > 0 and "driver" in cp
    # at least the two pediatric models contribute to the plasma comparison
    ids = {r.model_id for r in cmp.included}
    assert "hypnotics_iv.propofol.kataria_1994" in ids
    assert "hypnotics_iv.propofol.paedfusor_2005" in ids


def test_compare_greys_out_envelope_violator(ds):
    obese = dict(age=40, weight=140, height=172, sex="M")
    cmp = hypnos.compare(ds, drug="propofol", patient=obese, schedule=SCHED, t=T)
    excluded_ids = {e["model_id"] for e in cmp.excluded}
    assert SCHNIDER in excluded_ids
    assert MARSH in {r.model_id for r in cmp.included}


def test_predicate_evaluator_is_sandboxed(ds):
    from hypnos.simulate import _eval_predicate

    assert _eval_predicate("bmi > 42", {"bmi": 47.0}) is True
    assert _eval_predicate("age > 65", {"age": 50}) is False
    with pytest.raises(ValueError):
        _eval_predicate("__import__('os')", {})
