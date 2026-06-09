"""Plasma decrement-time (offset) analysis + the expm-cache correctness check."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.export.registry import instantiate
from hypnos.reference import Dosing, simulate

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
MINTO = "opioids.remifentanil.minto_1997"
FENTANYL = "opioids.fentanyl.shafer_1990"
ADULT = dict(age=50, weight=77, height=177, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def _decr(ds, model, infusion, dur):
    return hypnos.decrement_time(ds, model, patient=ADULT, infusion=infusion, duration=dur).decrement_min


def test_propofol_decrement_lengthens_with_duration(ds):
    # context-sensitivity: longer infusion -> slower plasma decline (accumulation)
    d10 = _decr(ds, SCHNIDER, "10 mg/kg/h", 10)
    d60 = _decr(ds, SCHNIDER, "10 mg/kg/h", 60)
    d600 = _decr(ds, SCHNIDER, "10 mg/kg/h", 600)
    assert d10 < d60 < d600
    assert all(np.isfinite([d10, d60, d600]))


def test_remifentanil_decrement_is_near_flat(ds):
    # remifentanil's celebrated context-insensitivity: the decrement plateaus
    d60 = _decr(ds, MINTO, "0.25 mcg/kg/min", 60)
    d600 = _decr(ds, MINTO, "0.25 mcg/kg/min", 600)
    assert abs(d600 - d60) < 0.6
    # and it is flatter than propofol over the same duration span
    p60 = _decr(ds, SCHNIDER, "10 mg/kg/h", 60)
    p600 = _decr(ds, SCHNIDER, "10 mg/kg/h", 600)
    assert (d600 - d60) < (p600 - p60)


def test_decrement_validation_and_kernel_pending(ds):
    with pytest.raises(ValueError):
        hypnos.decrement_time(ds, SCHNIDER, patient=ADULT, infusion="10 mg/kg/h", duration=60, fraction=1.0)
    with pytest.raises(NotImplementedError):  # fentanyl kernel pending
        hypnos.decrement_time(ds, FENTANYL, patient=ADULT, infusion="2 mcg/kg/min", duration=60)


def test_conc_at_stop_positive(ds):
    dt = hypnos.decrement_time(ds, SCHNIDER, patient=ADULT, infusion="10 mg/kg/h", duration=60)
    assert dt.conc_at_stop > 0


def test_cli_decrement(capsys):
    rc = main(["decrement", SCHNIDER, "--infusion", "10 mg/kg/h", "--duration", "240",
               "--age", "50", "--weight", "77", "--height", "177", "--sex", "M"])
    assert rc == 0
    assert "decrement time" in capsys.readouterr().out


# --- expm-cache correctness (the solver optimization must not change results) ---
def test_expm_cache_grid_independent(ds):
    # the memoized propagator must give identical concentrations regardless of grid
    # density and irregular sample points (different dt values share/refresh the cache)
    p = instantiate(ds[SCHNIDER], ADULT)
    dosing = Dosing(boluses=((0.0, 150.0),), infusions=((0.0, 6.0), (20.0, 0.0)))
    coarse = simulate(p, dosing, np.array([0.0, 5.0, 20.0, 45.123, 60.0]))
    fine = simulate(p, dosing, np.linspace(0.0, 60.0, 6001))
    # compare at the shared/known times
    for tj, cp_coarse in zip(coarse.t, coarse.cp):
        k = int(np.argmin(np.abs(fine.t - tj)))
        if abs(fine.t[k] - tj) < 1e-6:
            assert abs(fine.cp[k] - cp_coarse) < 1e-9
