import numpy as np

from hypnos.reference import (
    Dosing,
    MicroParams,
    lbm_james,
    sigmoid_emax,
    simulate,
    simulate_numeric,
)


def test_james_lbm_known_values():
    # 77 kg, 177 cm male
    lbm = lbm_james(77, 177, "M")
    assert abs(lbm - (1.1 * 77 - 128 * (77 / 177) ** 2)) < 1e-9
    # female uses different constants
    assert lbm_james(60, 165, "F") != lbm_james(60, 165, "M")


def test_one_compartment_bolus_closed_form():
    # Single compartment, instantaneous bolus -> C(t) = (D/V) e^{-k10 t}
    p = MicroParams(V1=10.0, k10=0.2, n_compartments=1)
    dosing = Dosing(boluses=((0.0, 100.0),))
    t = np.linspace(0, 30, 200)
    traj = simulate(p, dosing, t)
    expected = (100.0 / 10.0) * np.exp(-0.2 * t)
    assert np.allclose(traj.cp, expected, rtol=1e-7, atol=1e-9)


def test_analytic_matches_numeric_three_comp():
    # Schnider-like params; exact matrix-exp solver vs independent scipy integration.
    p = MicroParams.from_volumes_clearances(
        V1=4.27, Cl1=1.79, V2=18.9, Cl2=1.29, V3=238.0, Cl3=0.836, ke0=0.456
    )
    dosing = Dosing(boluses=((0.0, 150.0),), infusions=((0.0, 6.0), (20.0, 0.0)))
    t = np.linspace(0, 60, 121)
    a = simulate(p, dosing, t)
    n = simulate_numeric(p, dosing, t)
    assert np.allclose(a.cp, n.cp, rtol=1e-4, atol=1e-4)
    assert np.allclose(a.ce, n.ce, rtol=1e-4, atol=1e-4)


def test_effect_site_lags_plasma():
    p = MicroParams.from_volumes_clearances(
        V1=4.27, Cl1=1.79, V2=18.9, Cl2=1.29, V3=238.0, Cl3=0.836, ke0=0.456
    )
    dosing = Dosing(boluses=((0.0, 150.0),))
    t = np.linspace(0, 10, 200)
    traj = simulate(p, dosing, t)
    # Cp peaks at t=0 (bolus); Ce peaks later.
    assert np.argmax(traj.cp) == 0
    assert np.argmax(traj.ce) > 0


def test_sigmoid_emax_bounds():
    ce = np.array([0.0, 4.7, 1e6])
    e = sigmoid_emax(ce, E0=93, Emax=93, Ce50=4.7, gamma=1.43)
    assert abs(e[0] - 93.0) < 1e-9          # no drug -> baseline
    assert abs(e[1] - 93.0 / 2) < 1e-9       # at Ce50 -> half maximal effect
    assert e[2] < 0.5                         # saturating -> approaches 0
