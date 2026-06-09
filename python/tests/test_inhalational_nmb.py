import numpy as np
import pytest

import hypnos
from hypnos.reference import (
    alveolar_washin,
    alveolar_washout,
    mac_age_corrected,
    sigmoid_emax,
)

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


def test_alveolar_washin_kernel_math():
    # FA/FI(0) == 0; monotonic increasing; bounded by the plateau = V̇_A/(V̇_A+λQ̇)
    t = np.linspace(0, 10, 101)
    fa_fi, plateau, tau = alveolar_washin(0.65, t, alveolar_ventilation=4.0, frc=2.5, cardiac_output=5.0)
    assert abs(fa_fi[0]) < 1e-12
    assert np.all(np.diff(fa_fi) >= -1e-15)            # monotonic non-decreasing
    assert np.all(fa_fi <= plateau + 1e-12)            # never exceeds the plateau
    assert abs(plateau - 4.0 / (4.0 + 0.65 * 5.0)) < 1e-12
    assert abs(tau - 2.5 / (4.0 + 0.65 * 5.0)) < 1e-12
    # one time constant -> ~63.2% of the plateau
    one_tau, p, _ = alveolar_washin(0.65, np.array([tau]))
    assert abs(one_tau[0] / p - (1 - np.exp(-1))) < 1e-9


def test_washin_lower_solubility_is_faster(ds):
    # less soluble -> higher early FA/FI plateau (the discriminating, monotonic quantity)
    des = hypnos.washin(ds, DES).plateau
    n2o = hypnos.washin(ds, N2O).plateau
    sevo = hypnos.washin(ds, SEVO).plateau
    iso = hypnos.washin(ds, ISO).plateau
    assert des > n2o > sevo > iso          # exactly the blood:gas ordering, inverted


def test_washin_comparison_sorted_fastest_first(ds):
    rows = hypnos.washin_comparison(ds)
    assert [r.agent_id.split(".")[1] for r in rows][0] == "desflurane"
    assert rows[-1].agent_id.split(".")[1] == "isoflurane"
    plateaus = [r.plateau for r in rows]
    assert plateaus == sorted(plateaus, reverse=True)
    assert all(0.0 < r.plateau < 1.0 for r in rows)


def test_washin_rejects_non_physicochemical(ds):
    with pytest.raises(ValueError):
        hypnos.washin(ds, "hypnotics_iv.propofol.schnider_1998")


def test_cli_washin(capsys):
    from hypnos.cli import main
    assert main(["washin"]) == 0
    out = capsys.readouterr().out
    assert "desflurane" in out and "isoflurane" in out
    assert main(["washin", "--agent", "sevoflurane"]) == 0
    assert "plateau" in capsys.readouterr().out


def test_alveolar_washout_kernel_math():
    # FA/FA₀(0) == 1; monotonic decreasing; bounded below by the floor = λQ̇/(V̇_A+λQ̇)
    t = np.linspace(0, 10, 101)
    fa, floor, tau = alveolar_washout(0.65, t, alveolar_ventilation=4.0, frc=2.5, cardiac_output=5.0)
    assert abs(fa[0] - 1.0) < 1e-12
    assert np.all(np.diff(fa) <= 1e-15)                 # monotonic non-increasing
    assert np.all(fa >= floor - 1e-12)                  # never falls below the floor
    assert abs(floor - 0.65 * 5.0 / (4.0 + 0.65 * 5.0)) < 1e-12
    assert abs(tau - 2.5 / (4.0 + 0.65 * 5.0)) < 1e-12
    # wash-out floor and wash-in plateau are exact complements; same time constant
    _, plateau, tau_in = alveolar_washin(0.65, t)
    assert abs(floor - (1.0 - plateau)) < 1e-12
    assert abs(tau - tau_in) < 1e-12
    # one time constant -> decayed ~63.2% of the way from 1 toward the floor
    one_tau, f, _ = alveolar_washout(0.65, np.array([tau]))
    assert abs((1.0 - one_tau[0]) / (1.0 - f) - (1 - np.exp(-1))) < 1e-9


def test_washout_lower_solubility_is_faster(ds):
    # less soluble -> lower elimination floor -> more complete, faster wash-out
    des = hypnos.washout(ds, DES).floor
    n2o = hypnos.washout(ds, N2O).floor
    sevo = hypnos.washout(ds, SEVO).floor
    iso = hypnos.washout(ds, ISO).floor
    assert des < n2o < sevo < iso          # exactly the blood:gas ordering


def test_washout_comparison_sorted_fastest_first(ds):
    rows = hypnos.washout_comparison(ds)
    assert rows[0].agent_id.split(".")[1] == "desflurane"
    assert rows[-1].agent_id.split(".")[1] == "isoflurane"
    floors = [r.floor for r in rows]
    assert floors == sorted(floors)
    assert all(0.0 < r.floor < 1.0 for r in rows)


def test_washout_rejects_non_physicochemical(ds):
    with pytest.raises(ValueError):
        hypnos.washout(ds, "hypnotics_iv.propofol.schnider_1998")


def test_cli_washout(capsys):
    from hypnos.cli import main
    assert main(["washout"]) == 0
    out = capsys.readouterr().out
    assert "desflurane" in out and "isoflurane" in out
    assert main(["washout", "--agent", "sevoflurane"]) == 0
    assert "floor" in capsys.readouterr().out


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
