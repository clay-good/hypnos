"""Ketamine Kamp 2020: the first dissociative. Pooled arterial population model, PK-only,
linear weight scaling as in the one open reproduction (pending_human_review, partial_match)."""
import numpy as np
import pytest

import hypnos
from hypnos.cli import main
from hypnos.export.registry import KERNELS
from hypnos.filter import pk_drugs
from hypnos.presets import DEFAULT_SCHEDULES

KAMP = "dissociatives.ketamine.kamp_2020"
ADULT = dict(age=40, weight=70, height=175, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_kernel_reproduces_record_at_70kg(ds):
    m = ds[KAMP]
    vc = KERNELS[m.kernel_function](ADULT).as_volumes_clearances()
    for p in m.parameters:
        assert vc[p.symbol] == pytest.approx(p.central, rel=1e-3), p.symbol


def test_inside_the_abstract_meta_analysis_cis():
    # Kamp 2020 abstract: Vd 252 L/70 kg (95% CI 200-304); CL 79 L/h (95% CI 69-90)
    vc = KERNELS["ketamine_kamp_2020"](ADULT).as_volumes_clearances()
    assert 200 <= vc["V1"] + vc["V2"] + vc["V3"] <= 304
    assert 69 <= vc["Cl1"] * 60 <= 90


def test_linear_weight_scaling():
    a = KERNELS["ketamine_kamp_2020"](ADULT).as_volumes_clearances()
    b = KERNELS["ketamine_kamp_2020"](dict(ADULT, weight=140)).as_volumes_clearances()
    for k in ("V1", "V3", "Cl1", "Cl3"):
        assert b[k] == pytest.approx(2 * a[k])


def test_pk_only_no_invented_ke0(ds):
    assert KERNELS["ketamine_kamp_2020"](ADULT).ke0 == 0.0
    assert "ke0" not in {p.symbol for p in ds[KAMP].parameters}


def test_provenance_is_honest(ds):
    ex = ds[KAMP].raw["extraction"]
    assert ex["review_status"] == "pending_human_review"
    assert ex["source_review"]["outcome"] == "partial_match"
    assert ex["source_review"]["human_verified"] is False
    # the venous-vs-arterial caveat is carried as a cited, prose (non-predicate) failure mode
    fm = ds[KAMP].raw["known_failure_modes"][0]
    assert fm["predicate"] is None and "ARTERIAL" in fm["behavior"]


def test_simulatable_with_own_default(ds):
    assert "ketamine" in pk_drugs(ds) and "ketamine" in DEFAULT_SCHEDULES
    r = hypnos.simulate(ds, KAMP, patient=ADULT, schedule=DEFAULT_SCHEDULES["ketamine"],
                        t=np.linspace(0, 120, 121))
    assert r.tier == "C" and r.concentration_unit == "ng/mL" and r.cp_peak > 0


def test_cli_labels_pk_only_effect_site(capsys):
    assert main(["simulate", KAMP, "--age", "40", "--weight", "70", "--height", "175",
                 "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "Ce peak: n/a (PK-only" in out and "Ce peak: 0.000" not in out
