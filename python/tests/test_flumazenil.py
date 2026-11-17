"""Flumazenil: competitive benzodiazepine reversal of midazolam (Fiset 1995), the v0.5 B1 agent
that was blocked until a benzodiazepine kernel existed; plus agonist infusions in competitive
reversal and the resedation computed from both drugs' kinetics."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.simulate import _antagonist_micro, simulate_reversal
from hypnos.validate import validate_dataset

FLU = "reversal.flumazenil.competitive"
MID = "hypnotics_iv.midazolam.buhrer_1990"
MAN = dict(age=40, weight=80, height=177, sex="M")
T = np.linspace(0, 300, 601)


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_targets_resolve_and_validate(ds):
    rm = ds[FLU].reversal_mechanism
    assert rm["targets"] == [MID] and ds[MID].kernel_implemented
    assert not [p for p in validate_dataset(ds) if FLU in p]


def test_fiset_kinetics_reproduce_table_half_lives(ds):
    mp = _antagonist_micro(ds[FLU])
    a = np.array([[-(mp.k10 + mp.k12 + mp.k13), mp.k21, mp.k31],
                  [mp.k12, -mp.k21, 0], [mp.k13, 0, -mp.k31]])
    half = sorted(np.log(2) / -np.linalg.eigvals(a).real)
    for got, ref in zip(half, [0.9, 11.4, 63.2]):
        assert got == pytest.approx(ref, rel=0.06)
    vc = mp.as_volumes_clearances()
    assert vc["V1"] + vc["V2"] + vc["V3"] == pytest.approx(61, abs=0.1)


def test_resedation_on_a_midazolam_infusion(ds):
    r = simulate_reversal(ds, FLU, patient=MAN, reversal_dose="0.5 mg",
                          agonist_dose="0.1 mg/kg/h", t=T, reversal_time_min=60)
    share = r.agonist_ce_effective / np.where(r.agonist_ce > 0, r.agonist_ce, 1)
    i60 = int(np.argmin(np.abs(T - 60)))
    assert share[i60] == pytest.approx(1.0)
    assert share.min() < 0.6           # flumazenil roughly halves the effective midazolam
    assert share[-1] > 0.95            # ...and the reversal is gone by 5 h: resedation
    assert r.renarcotization
    assert any("of the un-reversed Ce" in w for w in r.warnings)


def test_bolus_agonist_still_supported(ds):
    r = simulate_reversal(ds, FLU, patient=MAN, reversal_dose="0.5 mg", agonist_dose="0.1 mg/kg",
                          t=T, reversal_time_min=30)
    assert r.agonist_ce is not None and r.agonist_ce[0] == 0.0


def test_infusion_rejected_for_non_competitive_mechanisms(ds):
    with pytest.raises(ValueError, match="only for competitive"):
        simulate_reversal(ds, "reversal.sugammadex.encapsulation", patient=MAN,
                          reversal_dose="2 mg/kg", agonist_dose="0.6 mg/kg/h", t=T)


def test_ki_scale_is_the_published_reversal_ec50(ds):
    ant = ds[FLU].reversal_mechanism["antagonism"]
    assert ant["ki"] == 13.7 and ant["ki_units"] == "ng/mL"


def test_cli_reverse_flumazenil_on_infusion(capsys):
    assert main(["reverse", "--agent", FLU, "--dose", "0.5 mg", "--agonist-dose", "0.1 mg/kg/h",
                 "--at", "60", "--tmax", "300", "--age", "40", "--weight", "80", "--height", "177",
                 "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "midazolam" in out and "renarcotization" in out and "True" in out
