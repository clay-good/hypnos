
import numpy as np
import pytest

import hypnos
from hypnos.export import export_model
from hypnos.export.registry import instantiate
from hypnos.reference import Dosing, MicroParams, greco_response_surface, sigmoid_emax, simulate

REMI = "opioids.remifentanil.minto_1997"
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
SURFACE = "interactions.propofol_remifentanil.greco_bis"

PROP_SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
REMI_SCHED = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]
T = np.linspace(0, 30, 181)


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_remifentanil_simulation_sane(ds):
    patient = dict(age=40, weight=70, height=170, sex="M")
    res = hypnos.simulate(ds, REMI, patient=patient, schedule=REMI_SCHED, t=T)
    assert res.tier == "B"
    # 1 mcg/kg into ~5 L central -> ~10-15 ng/mL (0.010-0.015 ug/mL) plasma peak
    assert 0.008 < res.cp_peak < 0.020
    assert res.ce_peak < res.cp_peak  # effect site lags/attenuates the spike


def test_remifentanil_obese_envelope_violation(ds):
    obese = dict(age=40, weight=150, height=170, sex="M")  # BMI ~52
    res = hypnos.simulate(ds, REMI, patient=obese, schedule=REMI_SCHED, t=T)
    assert res.tier == "D"
    assert res.excluded


def test_greco_collapses_to_single_drug_at_zero_opioid():
    ce_prop = np.array([0.0, 4.0, 8.0])
    zero = np.zeros_like(ce_prop)
    surf = greco_response_surface(ce_prop, zero, E0=93, Emax=93, Ce50_a=4.0,
                                  Ce50_b=0.04, alpha=4.0, gamma=1.4)
    sig = sigmoid_emax(ce_prop, E0=93, Emax=93, Ce50=4.0, gamma=1.4)
    assert np.allclose(surf, sig, atol=1e-12)


def test_greco_synergy_lowers_effect():
    ce_prop = np.array([4.0])
    no_opioid = greco_response_surface(ce_prop, np.array([0.0]), 93, 93, 4.0, 0.04, 4.0, 1.4)
    with_opioid = greco_response_surface(ce_prop, np.array([0.012]), 93, 93, 4.0, 0.04, 4.0, 1.4)
    assert with_opioid[0] < no_opioid[0]  # adding the opioid deepens hypnosis


def test_simulate_interaction_shows_synergy(ds):
    patient = dict(age=40, weight=70, height=170, sex="M")
    ir = hypnos.simulate_interaction(
        ds, SURFACE, pk_a=SCHNIDER, pk_b=REMI, patient=patient,
        schedule_a=PROP_SCHED, schedule_b=REMI_SCHED, t=T,
    )
    alone = hypnos.simulate(ds, SCHNIDER, patient=patient, schedule=PROP_SCHED, t=T,
                            pd_model="pd_effect.propofol.bis_sigmoid")
    assert ir.tier == "C"  # surface is Tier C -> worst input wins
    assert ir.effect_min < float(np.min(alone.effect))  # combination is deeper


def test_interaction_inherits_worst_envelope_tier(ds):
    obese = dict(age=40, weight=150, height=170, sex="M")  # propofol+remi both out of envelope
    ir = hypnos.simulate_interaction(
        ds, SURFACE, pk_a=SCHNIDER, pk_b=REMI, patient=obese,
        schedule_a=PROP_SCHED, schedule_b=REMI_SCHED, t=T,
    )
    assert ir.tier == "D"
    assert any("ENVELOPE" in w for w in ir.warnings)


def test_interaction_rejects_non_surface(ds):
    with pytest.raises(ValueError):
        hypnos.simulate_interaction(ds, SCHNIDER, pk_a=SCHNIDER, pk_b=REMI,
                                    patient=dict(age=40, weight=70, height=170, sex="M"),
                                    schedule_a=PROP_SCHED, schedule_b=REMI_SCHED, t=T)


def _parse_params(text):
    line = next(ln for ln in text.splitlines() if "hypnos.params:" in ln)
    kv = dict(tok.split("=") for tok in line.split("hypnos.params:")[1].split())
    return MicroParams(**{k: float(v) for k, v in kv.items()})


@pytest.mark.parametrize("fmt", ["rxode2", "pumas"])
def test_r_julia_exports_round_trip(ds, fmt):
    m = ds[REMI]
    patient = {"age": 50, "weight": 77, "height": 177, "sex": "M"}
    _, text = export_model(fmt, m, ds, patient)
    assert "PROHIBITED" in text
    recovered = _parse_params(text)
    direct = instantiate(m, patient)
    dosing = Dosing(boluses=((0.0, 0.1),), infusions=((0.0, 0.02),))
    t = np.linspace(0, 30, 120)
    assert np.allclose(simulate(direct, dosing, t).cp,
                       simulate(recovered, dosing, t).cp, rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("fmt", ["rxode2", "pumas"])
def test_r_julia_structural_tokens(ds, fmt):
    _, text = export_model(fmt, ds[SCHNIDER], ds)
    if fmt == "rxode2":
        assert "rxode2({" in text and "d/dt(A1)" in text
    else:
        assert "@model begin" in text and "@dynamics begin" in text
