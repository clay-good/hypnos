"""Fentanyl Shafer 1990: the kernel wired to its sourced constants (pending_human_review),
the Shibutani 2004 obesity failure mode, and naloxone renarcotization computed from the
agonist's own kinetics."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.export.registry import KERNELS
from hypnos.simulate import simulate_reversal

FENTANYL = "opioids.fentanyl.shafer_1990"
REMI = "opioids.remifentanil.minto_1997"
NALOXONE = "reversal.naloxone.competitive"
ADULT = dict(age=40, weight=70, height=175, sex="M")
BOLUS = [("bolus", 0.0, "2 mcg/kg")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_kernel_reproduces_curated_parameters_exactly(ds):
    # the kernel computes exactly the sourced record values (never a second copy that can drift)
    m = ds[FENTANYL]
    assert m.kernel_implemented and m.kernel_function == "fentanyl_shafer_1990"
    vc = KERNELS["fentanyl_shafer_1990"](ADULT).as_volumes_clearances()
    for p in m.parameters:
        if p.symbol == "ke0":
            assert KERNELS["fentanyl_shafer_1990"](ADULT).ke0 == pytest.approx(p.central)
        else:
            assert vc[p.symbol] == pytest.approx(p.central, rel=1e-9), p.symbol


def test_opentiva_micro_constants_reproduced():
    # the cross-check source: opentiva k10 .0827, k12 .471, k13 .225, k21 .102, k31 .006
    mp = KERNELS["fentanyl_shafer_1990"](ADULT)
    for k, ref in [("k10", 0.0827), ("k12", 0.471), ("k13", 0.225), ("k21", 0.102), ("k31", 0.006)]:
        assert getattr(mp, k) == pytest.approx(ref, abs=6e-4), k


def test_still_pending_human_review_never_verified(ds):
    assert ds[FENTANYL].raw["extraction"]["review_status"] == "pending_human_review"
    assert ds[FENTANYL].raw["extraction"]["source_review"]["human_verified"] is False


def test_bolus_time_to_peak_and_units(ds):
    t = np.linspace(0, 30, 3001)
    r = hypnos.simulate(ds, FENTANYL, patient=ADULT, schedule=BOLUS, t=t)
    assert r.tier == "C" and r.concentration_unit == "ng/mL" and not r.warnings
    tp = t[int(np.argmax(r.ce))]
    assert 2.5 < tp < 4.0   # ke0 0.147/min -> ~3.2 min (literature ~3.6 min)
    # slower onset than remifentanil (Minto ke0 ~0.6/min)
    rr = hypnos.simulate(ds, REMI, patient=ADULT, schedule=BOLUS, t=t)
    assert t[int(np.argmax(rr.ce))] < tp


def test_context_sensitive_decrement_unlike_remifentanil(ds):
    # fentanyl's 50% decrement lengthens steeply with infusion duration; remifentanil's does not
    f60 = hypnos.decrement_time(ds, FENTANYL, patient=ADULT, infusion="2 mcg/kg/h", duration=60).decrement_min
    f480 = hypnos.decrement_time(ds, FENTANYL, patient=ADULT, infusion="2 mcg/kg/h", duration=480).decrement_min
    r60 = hypnos.decrement_time(ds, REMI, patient=ADULT, infusion="0.1 mcg/kg/min", duration=60).decrement_min
    r480 = hypnos.decrement_time(ds, REMI, patient=ADULT, infusion="0.1 mcg/kg/min", duration=480).decrement_min
    assert f480 > 5 * f60
    assert (f480 - f60) > 20 * abs(r480 - r60)


def test_obesity_failure_mode_tiers_down(ds):
    t = np.linspace(0, 60, 121)
    heavy = dict(ADULT, weight=88)
    r = hypnos.simulate(ds, FENTANYL, patient=heavy, schedule=BOLUS, t=t)
    assert r.tier == "D"
    assert any("Shibutani" in w and "OVERESTIMATES" in w for w in r.warnings)
    fm = ds[FENTANYL].raw["known_failure_modes"][0]
    assert fm["citation"] == "shibutani-2004-fentanyl-obese" and fm["citation"] in ds.citations


def test_naloxone_renarcotization_from_agonist_kinetics(ds):
    t = np.linspace(0, 240, 481)
    r = simulate_reversal(ds, NALOXONE, patient=ADULT, reversal_dose="0.4 mg",
                          agonist_dose="2 mcg/kg", t=t, reversal_time_min=10.0)
    assert r.agonist_id == FENTANYL and r.agonist_unit == "ng/mL"
    ag, eff = r.agonist_ce, r.agonist_ce_effective
    # identical before naloxone arrives; suppressed afterwards
    assert eff[10] == pytest.approx(ag[10])
    assert eff[40] < 0.5 * ag[40]
    # the adjusted curve reaches a nadir and then rises again as naloxone clears
    after = eff[t >= 10]
    assert after[-1] > 1.05 * after.min()
    assert r.renarcotization
    assert any("RENARCOTIZATION (agonist kinetics)" in w for w in r.warnings)
    np.testing.assert_allclose(r.free_ce_no_reversal, ag)


def test_short_horizon_shows_no_agonist_relapse(ds):
    # within 60 min naloxone is still near peak: no relapse of the adjusted agonist curve
    t = np.linspace(0, 60, 121)
    r = simulate_reversal(ds, NALOXONE, patient=ADULT, reversal_dose="0.4 mg",
                          agonist_dose="2 mcg/kg", t=t, reversal_time_min=10.0)
    assert not any("agonist kinetics" in w for w in r.warnings)


def test_cli_reverse_prints_agonist_table(capsys):
    assert main(["reverse", "--agent", NALOXONE, "--dose", "0.4 mg", "--agonist-dose", "2 mcg/kg",
                 "--tmax", "240", "--age", "40", "--weight", "70", "--height", "175", "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "antagonist-adjusted" in out and "ng/mL" in out and "renarcotization" in out


def test_cli_simulate_fentanyl(capsys):
    assert main(["simulate", FENTANYL, "--age", "40", "--weight", "70", "--height", "175",
                 "--sex", "M"]) == 0
    assert "fentanyl" in capsys.readouterr().out.lower()
