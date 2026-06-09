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


def test_eleveld_kernel_refuses(ds):
    patient = dict(age=50, weight=77, height=177, sex="M")
    with pytest.raises(NotImplementedError):
        hypnos.simulate(ds, ELEVELD, patient=patient, schedule=SCHED, t=T)


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
    assert MARSH in ids and SCHNIDER in ids
    assert any(u["model_id"] == ELEVELD for u in cmp.unavailable)
    assert cmp.divergence["ce"]["max_abs"] > 0.5  # the models genuinely disagree


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
