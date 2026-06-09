"""Conventional concentration units + drug-appropriate CLI default dosing."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.presets import default_schedule_for as _default_schedule_for
from hypnos.models import concentration_factor

T = np.linspace(0, 30, 181)
REMI = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]
PROP = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_concentration_factor():
    assert concentration_factor("ng/mL") == 1000.0
    assert concentration_factor("ug/mL") == 1.0
    assert concentration_factor("mg/L") == 1.0
    assert concentration_factor(None) == 1.0
    assert concentration_factor("weird") == 1.0  # unknown -> identity


def test_remifentanil_reports_ng_per_ml(ds):
    res = hypnos.simulate(ds, "opioids.remifentanil.minto_1997",
                          patient=dict(age=40, weight=70, height=175, sex="M"), schedule=REMI, t=T)
    assert res.concentration_unit == "ng/mL"
    assert res.conc_factor == 1000.0
    assert abs(res.cp_peak_display - res.cp_peak * 1000.0) < 1e-9
    assert 5.0 < res.cp_peak_display < 50.0   # clinically sensible ng/mL (not a µg/mL artifact)


def test_propofol_reports_ug_per_ml(ds):
    res = hypnos.simulate(ds, "hypnotics_iv.propofol.schnider_1998",
                          patient=dict(age=50, weight=77, height=177, sex="M"), schedule=PROP, t=T)
    assert res.concentration_unit == "ug/mL"
    assert res.cp_peak_display == res.cp_peak   # factor 1


def test_compare_carries_unit(ds):
    cmp = hypnos.compare(ds, drug="remifentanil",
                         patient=dict(age=40, weight=75, height=178, sex="M"), schedule=REMI, t=T)
    assert cmp.concentration_unit == "ng/mL" and cmp.conc_factor == 1000.0


def test_cli_default_schedule_is_drug_appropriate():
    # propofol keeps mg/kg; remifentanil uses mcg/kg (a 2 mg/kg remi dose would be a ~1000x overdose)
    assert _default_schedule_for("propofol")[0][2] == "2 mg/kg"
    assert "mcg/kg" in _default_schedule_for("remifentanil")[0][2]
    assert _default_schedule_for("rocuronium")[0][2] == "0.6 mg/kg"


def test_cli_remifentanil_simulate_is_sensible(capsys):
    rc = main(["simulate", "opioids.remifentanil.minto_1997",
               "--age", "40", "--weight", "70", "--height", "175", "--sex", "M"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "ng/mL" in out
    cp = float(out.split("Cp peak:")[1].split("ng/mL")[0].strip())
    assert cp < 100.0   # not the old 26874 ng/mL propofol-dose artifact
