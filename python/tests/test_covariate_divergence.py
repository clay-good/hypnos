"""Tests for the v0.7 C1 covariate-equation divergence view.

Divergence *within* one model: overlay its predicted curve under each admissible
body-size equation, greying any outside its own validity envelope. The key
invariant — the model's OWN (verbatim) equation must reproduce the model's normal
prediction exactly — guards the backward-compatibility of the kernel override hook.
"""
import numpy as np
import pytest

import hypnos
from hypnos.simulate import covariate_divergence

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
KIM = "opioids.remifentanil.kim_2017"
MARSH = "hypnotics_iv.propofol.marsh_1991"

OBESE = dict(age=50, weight=130, height=170, sex="M")     # BMI ~45
NORMAL = dict(age=50, weight=75, height=178, sex="M")     # BMI ~24


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_own_equation_leads_and_is_verbatim(ds):
    cd = covariate_divergence(ds, SCHNIDER, patient=NORMAL)
    assert cd.derived_equation == "james_1976"
    assert cd.quantity == "lbm"
    assert cd.by_equation[0].verbatim is True
    assert cd.by_equation[0].equation_id == "james_1976"
    assert sum(c.verbatim for c in cd.by_equation) == 1


def test_own_curve_equals_plain_simulate(ds):
    """The verbatim curve must reproduce the model's normal prediction exactly —
    the override hook is invisible when the own equation is used (backward compat)."""
    t = np.linspace(0, 60, 361)
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    plain = hypnos.simulate(ds, SCHNIDER, patient=OBESE, schedule=sched, t=t)
    cd = covariate_divergence(ds, SCHNIDER, patient=OBESE, schedule=sched, t=t)
    np.testing.assert_allclose(cd.own.ce, plain.ce, rtol=0, atol=0)
    np.testing.assert_allclose(cd.own.cp, plain.cp, rtol=0, atol=0)


def test_obese_greys_own_james_and_substitutions_stay(ds):
    cd = covariate_divergence(ds, SCHNIDER, patient=OBESE)
    own = cd.own
    assert own.in_envelope is False
    assert own.inverted is True
    subs = [c for c in cd.by_equation if not c.verbatim]
    assert subs and all(c.in_envelope for c in subs)        # FFM equations valid at BMI 45


def test_normal_patient_all_in_envelope_small_spread(ds):
    cd = covariate_divergence(ds, SCHNIDER, patient=NORMAL)
    assert all(c.in_envelope for c in cd.by_equation)
    assert cd.own.inverted is False
    assert cd.divergence["max_rel"] < 0.05                  # equations agree closely in-range


def test_substitution_actually_moves_the_curve_in_obesity(ds):
    cd = covariate_divergence(ds, SCHNIDER, patient=OBESE)
    jan = next(c for c in cd.by_equation if c.equation_id == "janmahasatian_2005")
    # the FFM substitution gives a materially different prediction than the inverted James
    assert abs(jan.ce_peak - cd.own.ce_peak) > 0.1
    assert cd.divergence["max_rel"] > 0.03
    assert {cd.divergence["driver"]["high"], cd.divergence["driver"]["low"]} >= {"james_1976"}


def test_pk_only_model_uses_cp_key(ds):
    cd = covariate_divergence(ds, KIM, patient=dict(age=40, weight=120, height=165, sex="M"))
    assert cd.key == "cp"                                    # Kim has no ke0 -> Ce is zero
    assert cd.quantity == "ffm"
    assert cd.derived_equation == "janmahasatian_2005"


def test_model_without_covariate_model_raises(ds):
    with pytest.raises(ValueError, match="no covariate_model"):
        covariate_divergence(ds, MARSH, patient=NORMAL)


def test_cli_covariate_divergence(ds, capsys):
    from hypnos.cli import main
    rc = main(["covariate-divergence", "--model", "propofol.schnider_1998",
               "--age", "50", "--weight", "130", "--height", "170", "--sex", "M"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "covariate-equation divergence" in out
    assert "INVERTED" in out and "james_1976" in out
    assert "failure mode" in out
