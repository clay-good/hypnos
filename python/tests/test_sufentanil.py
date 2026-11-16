"""Sufentanil Gepts 1995: fentanyl family, second member. Printed volumes/clearances (abstract
verbatim), ke0 solved for the 5.6-min time to peak effect, and the rational-opioid-selection
contrast with fentanyl rendered forward."""
import numpy as np
import pytest

import hypnos
from hypnos.export.registry import KERNELS

GEPTS = "opioids.sufentanil.gepts_1995"
FENT = "opioids.fentanyl.shafer_1990"
ADULT = dict(age=40, weight=70, height=175, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_kernel_reproduces_record_and_abstract_cis(ds):
    m = ds[GEPTS]
    mp = KERNELS[m.kernel_function](ADULT)
    vc = mp.as_volumes_clearances()
    for p in m.parameters:
        got = mp.ke0 if p.symbol == "ke0" else vc[p.symbol]
        assert got == pytest.approx(p.central, rel=1e-9), p.symbol
        if p.symbol != "ke0":   # each printed value sits inside its own published 95% CI
            assert p.value["low"] <= p.central <= p.value["high"]


def test_v3_is_the_printed_value_not_the_rate_constant_artifact():
    # back-calculating V3 from k31 rounded to 0.0013 gives 251.9 L, outside the CI
    assert KERNELS["sufentanil_gepts_1995"](ADULT).as_volumes_clearances()["V3"] == 261.6


def test_ke0_reproduces_published_time_to_peak(ds):
    t = np.linspace(0, 15, 15001)
    r = hypnos.simulate(ds, GEPTS, patient=ADULT, schedule=[("bolus", 0.0, "0.2 mcg/kg")], t=t)
    assert t[int(np.argmax(r.ce))] == pytest.approx(5.6, abs=0.05)


def test_worst_input_wins_tier(ds):
    m = ds[GEPTS]
    assert {p.symbol: p.tier for p in m.parameters}["ke0"] == "C" and m.tier == "C"


def test_sufentanil_offset_far_less_context_sensitive_than_fentanyl(ds):
    def dec(mid, dur):
        return hypnos.decrement_time(ds, mid, patient=ADULT, infusion="1 mcg/kg/h",
                                     duration=dur).decrement_min
    # Shafer & Varvel 1991: after long infusions sufentanil's offset is much faster
    assert dec(GEPTS, 600) < 0.3 * dec(FENT, 600)
    assert dec(GEPTS, 60) < dec(FENT, 60)
