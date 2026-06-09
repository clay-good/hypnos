import numpy as np
import pytest

import hypnos
from hypnos.reference import mac_age_corrected, mac_fraction, sigmoid_emax

SEVO = "volatiles.sevoflurane.mac"
DES = "volatiles.desflurane.mac"
ISO = "volatiles.isoflurane.mac"
N2O = "volatiles.nitrous_oxide.mac"
ROC = "nmb_agents.rocuronium.wierda_1991"
TOF = "pd_effect.rocuronium.tof_sigmoid"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_mac_age_correction_formula():
    # at age 40, MAC == MAC40
    assert abs(mac_age_corrected(2.0, 40) - 2.0) < 1e-12
    # ~6% decrease per decade (Mapleson)
    ratio = mac_age_corrected(2.0, 50) / mac_age_corrected(2.0, 40)
    assert abs(ratio - 10 ** (-0.00269 * 10)) < 1e-12
    assert 0.93 < ratio < 0.95


def test_mac_decreases_with_age(ds):
    young = hypnos.mac(ds, SEVO, age=20)
    old = hypnos.mac(ds, SEVO, age=80)
    assert young.mac_age > old.mac_age
    assert abs(young.mac_age - mac_age_corrected(1.8, 20)) < 1e-9


def test_mac_fraction_at_one_mac(ds):
    # end-tidal == age-corrected MAC -> fraction 1.0
    r = hypnos.mac(ds, SEVO, age=40, end_tidal_pct=1.8)
    assert abs(r.mac_fraction - 1.0) < 1e-9


def test_nitrous_oxide_additivity(ds):
    # 0.5 MAC sevo (0.9 vol% at age 40) + 0.5 MAC N2O (52 vol%) -> 1.0 combined
    r = hypnos.mac(ds, SEVO, age=40, end_tidal_pct=0.9, n2o_end_tidal_pct=52.0)
    assert abs(r.mac_fraction - 0.5) < 1e-6
    assert abs(r.combined_mac_fraction - 1.0) < 1e-3


def test_mac_age_extrapolation_tiers_down(ds):
    # MAC age-correction is only valid for age > 1 y
    r = hypnos.mac(ds, SEVO, age=0.5)
    assert r.tier == "D"
    assert any("ENVELOPE" in w for w in r.warnings)


def test_mac_rejects_non_physicochemical(ds):
    with pytest.raises(ValueError):
        hypnos.mac(ds, "hypnotics_iv.propofol.schnider_1998", age=40)


def test_solubility_ordering(ds):
    # blood:gas governs speed; des < sevo < iso
    des = hypnos.mac(ds, DES, age=40).blood_gas
    sevo = hypnos.mac(ds, SEVO, age=40).blood_gas
    iso = hypnos.mac(ds, ISO, age=40).blood_gas
    assert des < sevo < iso


def test_train_of_four_sigmoid_shape():
    # T1 (twitch height) drops from 100 -> 50 at Ce50 -> ~0 at high Ce; steep gamma
    tof = sigmoid_emax(np.array([0.0, 0.823, 3.0]), 100.0, 100.0, 0.823, 4.5)
    assert abs(tof[0] - 100.0) < 1e-9
    assert abs(tof[1] - 50.0) < 1e-9
    assert tof[2] < 2.0


def test_rocuronium_kernel_pending_refuses(ds):
    with pytest.raises(NotImplementedError):
        hypnos.simulate(ds, ROC, patient=dict(age=40, weight=70, height=175, sex="M"),
                        schedule=[("bolus", 0.0, "0.6 mg/kg")], t=np.linspace(0, 60, 10))


def test_phase_d_subsystems_present(ds):
    s = hypnos.summary(ds)
    assert s["by_subsystem"].get("volatiles") == 4
    assert s["by_subsystem"].get("nmb_agents") == 1
    assert s["n_models"] >= 15
    assert s["n_drugs"] >= 9
