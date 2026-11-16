"""Etomidate Kaneda 2011: the last IV hypnotic on the v0.1 list. Abstract-verbatim PK."""
import math

import numpy as np
import pytest

import hypnos
from hypnos.export.registry import KERNELS
from hypnos.presets import DEFAULT_SCHEDULES

KANEDA = "hypnotics_iv.etomidate.kaneda_2011"
ADULT = dict(age=40, weight=70, height=175, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_kernel_reproduces_abstract_values(ds):
    mp = KERNELS["etomidate_kaneda_2011"](ADULT)
    vc = mp.as_volumes_clearances()
    assert (vc["V1"], vc["V2"], vc["Cl1"], vc["Cl2"]) == (4.45, 74.9, 0.63, 3.16)
    assert mp.n_compartments == 2
    assert mp.ke0 == pytest.approx(math.log(2) / 1.55, abs=1e-4)
    for p in ds[KANEDA].parameters:
        got = mp.ke0 if p.symbol == "ke0" else vc[p.symbol]
        assert got == pytest.approx(p.central, rel=1e-9)


def test_bolus_reaches_hypnotic_effect_site_quickly(ds):
    t = np.linspace(0, 10, 1001)
    r = hypnos.simulate(ds, KANEDA, patient=ADULT, schedule=DEFAULT_SCHEDULES["etomidate"], t=t)
    tp = t[int(np.argmax(r.ce))]
    assert 1.0 < tp < 2.5
    # the published BIS EC50 is 0.526 ug/mL: an induction bolus clearly exceeds it
    assert r.ce_peak > 0.526


def test_envelope_and_short_infusion_caveat(ds):
    m = ds[KANEDA]
    assert m.raw["applicability_envelope"]["derivation_n"] == 18
    fm = m.raw["known_failure_modes"][0]
    assert fm["predicate"] is None and "10 minutes" in fm["behavior"]
    r = hypnos.simulate(ds, KANEDA, patient=dict(ADULT, age=75), schedule=DEFAULT_SCHEDULES["etomidate"],
                        t=np.linspace(0, 10, 11))
    assert r.tier == "D"   # outside the 18-55 y envelope -> named, tiered down


def test_no_invented_pd_record(ds):
    # the BIS EC50/gamma are published but E0/Emax are not in the abstract: no PD record
    assert not [m for m in ds if m.purpose == "pd" and m.drug_name == "etomidate"]
