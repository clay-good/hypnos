"""Midazolam: the healthy-volunteer (Buhrer 1990) vs post-CABG ICU (Zomorodi 1998) divergence
pair, both printed in Zomorodi 1998 Table 3; plus the single-in-envelope-model compare
regression (it crashed the CLI with KeyError 'max_abs')."""
import math

import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.export.registry import KERNELS
from hypnos.presets import DEFAULT_SCHEDULES
from hypnos.filter import pk_drugs

BUHRER = "hypnotics_iv.midazolam.buhrer_1990"
ZOMORODI = "hypnotics_iv.midazolam.zomorodi_1998"
ADULT = dict(age=40, weight=70, height=175, sex="M")
ICU = dict(age=62, weight=85, height=175, sex="M")
T = np.linspace(0, 240, 241)


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


@pytest.mark.parametrize("mid", [BUHRER, ZOMORODI])
def test_kernel_reproduces_record_exactly(ds, mid):
    m = ds[mid]
    mp = KERNELS[m.kernel_function](ADULT)
    vc = mp.as_volumes_clearances()
    for p in m.parameters:
        got = mp.ke0 if p.symbol == "ke0" else vc[p.symbol]
        assert got == pytest.approx(p.central, rel=1e-9), p.symbol
    assert m.raw["extraction"]["review_status"] == "pending_human_review"
    assert m.raw["extraction"]["source_review"]["human_verified"] is False


def test_buhrer_table3_exponents_and_ke0():
    # Zomorodi Table 3 prints the Buhrer exponents 1.097 / 0.047 / 0.0031 /min
    mp = KERNELS["midazolam_buhrer_1990"](ADULT)
    k10, k12, k13, k21, k31 = mp.k10, mp.k12, mp.k13, mp.k21, mp.k31
    a = np.array([[-(k10 + k12 + k13), k21, k31], [k12, -k21, 0], [k13, 0, -k31]])
    lam = sorted(-np.linalg.eigvals(a).real, reverse=True)
    for got, ref in zip(lam, [1.097, 0.047, 0.0031]):
        assert got == pytest.approx(ref, rel=0.02)
    assert mp.ke0 == pytest.approx(math.log(2) / 5.6, abs=1e-4)


def test_steady_state_volumes_match_table3():
    for fn, vss in [("midazolam_buhrer_1990", 117.7), ("midazolam_zomorodi_1998", 330.2)]:
        vc = KERNELS[fn](ADULT).as_volumes_clearances()
        assert vc["V1"] + vc["V2"] + vc["V3"] == pytest.approx(vss, abs=0.1)


def test_zomorodi_omega_is_exact_lognormal_cv(ds):
    om = ds[ZOMORODI].bsv_omegas()
    assert set(om) == {"V1", "V2", "V3", "Cl1", "Cl2", "Cl3"}
    assert om["V3"] == pytest.approx(math.log(1 + 0.71 ** 2), abs=1e-6)


def test_midazolam_is_simulatable_with_its_own_default_dose(ds):
    assert "midazolam" in pk_drugs(ds)
    # never the propofol fallback (2 mg/kg of midazolam would be a ~40x overdose)
    assert DEFAULT_SCHEDULES["midazolam"][0][2] == "0.05 mg/kg"


def test_every_simulatable_drug_has_an_explicit_default(ds):
    # hardening: the propofol fallback must never silently dose another drug
    assert set(pk_drugs(ds)) <= set(DEFAULT_SCHEDULES)


def test_healthy_vs_icu_divergence(ds):
    cmp = hypnos.compare(ds, drug="midazolam", patient=ADULT,
                         schedule=DEFAULT_SCHEDULES["midazolam"], t=T)
    assert {r.model_id for r in cmp.included} == {BUHRER, ZOMORODI}
    d = cmp.divergence["cp"]
    assert d["max_rel"] > 1.0     # the two populations disagree by >100% somewhere


def test_elderly_icu_patient_greys_the_volunteer_model(ds):
    cmp = hypnos.compare(ds, drug="midazolam", patient=ICU,
                         schedule=DEFAULT_SCHEDULES["midazolam"], t=T, bands=True, seed=7)
    assert [r.model_id for r in cmp.included] == [ZOMORODI]
    assert {e["model_id"] for e in cmp.excluded} == {BUHRER}
    z = cmp.included[0]
    assert z.cp_quantiles is not None   # first non-propofol/remifentanil BSV band


def test_cli_compare_single_model_with_bands_does_not_crash(capsys):
    rc = main(["compare", "--drug", "midazolam", "--age", "62", "--weight", "85",
               "--height", "175", "--sex", "M", "--bands"])
    assert rc == 0
    assert "only one model in envelope" in capsys.readouterr().out


def test_ejection_fraction_standing_only_for_icu_model(ds):
    env = ds[ZOMORODI].raw["applicability_envelope"]
    assert env["ejection_fraction_pct"]["min"] == 25
    assert "ejection_fraction_pct" not in ds[BUHRER].raw["applicability_envelope"]
