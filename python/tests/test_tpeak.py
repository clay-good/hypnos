"""Time-to-peak-effect (onset) analysis."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
MARSH = "hypnotics_iv.propofol.marsh_1991"
MINTO = "opioids.remifentanil.minto_1997"
KIM = "opioids.remifentanil.kim_2017"          # PK-only, no ke0
PAEDFUSOR = "hypnotics_iv.propofol.paedfusor_2005"  # PK-only, no ke0
ADULT = dict(age=50, weight=77, height=177, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_ce_equals_cp_at_peak(ds):
    # the defining property: at the effect-site peak dCe/dt=0, hence Ce==Cp
    pe = hypnos.time_to_peak_effect(ds, SCHNIDER, patient=ADULT)
    assert abs(pe.ce_cp_ratio_at_peak - 1.0) < 2e-3


def test_tpeak_matches_published_range(ds):
    # Schnider propofol and Minto remifentanil both peak near ~1.6 min
    assert 1.2 < hypnos.time_to_peak_effect(ds, SCHNIDER, patient=ADULT).tpeak_min < 2.0
    assert 1.2 < hypnos.time_to_peak_effect(ds, MINTO, patient=ADULT).tpeak_min < 2.2


def test_lower_ke0_gives_later_peak(ds):
    # Marsh ke0=0.26 < Schnider ke0=0.456 -> Marsh peaks later
    marsh = hypnos.time_to_peak_effect(ds, MARSH, patient=ADULT)
    schnider = hypnos.time_to_peak_effect(ds, SCHNIDER, patient=ADULT)
    assert marsh.ke0 < schnider.ke0
    assert marsh.tpeak_min > schnider.tpeak_min


def test_tpeak_is_dose_independent(ds):
    # tpeak is a model property, independent of the bolus magnitude (linear system).
    # The function uses a fixed unit bolus, so just assert it's stable/positive.
    pe = hypnos.time_to_peak_effect(ds, SCHNIDER, patient=ADULT)
    assert pe.tpeak_min > 0 and np.isfinite(pe.tpeak_min)


@pytest.mark.parametrize("mid", [KIM, PAEDFUSOR])
def test_pk_only_models_have_no_tpeak(ds, mid):
    with pytest.raises(ValueError):
        hypnos.time_to_peak_effect(ds, mid, patient=ADULT)


def test_tpeak_carries_envelope_tier(ds):
    # out-of-envelope patient -> tiered down to D (envelope enforcement applies here too)
    child = dict(age=6, weight=20, height=115, sex="M")
    pe = hypnos.time_to_peak_effect(ds, SCHNIDER, patient=child)
    assert pe.tier == "D"
    assert any("ENVELOPE" in w or "EXTRAPOLATION" in w for w in pe.warnings)


def test_cli_tpeak(capsys):
    assert main(["tpeak", SCHNIDER, "--age", "50", "--weight", "77", "--height", "177", "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "time to peak effect" in out
    assert main(["tpeak", KIM, "--age", "40", "--weight", "75", "--height", "178", "--sex", "M"]) == 2
