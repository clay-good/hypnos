"""Alfentanil: Maitre 1987 (population; prospectively MDPE +53%) vs Scott 1985 (two-stage;
MDPE +1%) — the pair Raemer 1990 ran in one pump study. Pedigree is not performance."""
import numpy as np
import pytest

import hypnos
from hypnos.export import BUILDERS, export_model
from hypnos.export.registry import KERNELS
from hypnos.presets import DEFAULT_SCHEDULES

MAITRE = "opioids.alfentanil.maitre_1987"
SCOTT = "opioids.alfentanil.scott_1985"
MAN = dict(age=35, weight=70, height=175, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_maitre_abstract_reference_individual():
    # abstract: CL 0.356 L/min and Vc 7.77 L for a 70-kg male < 40 y (so 0.111 L/kg, not 0.11)
    vc = KERNELS["alfentanil_maitre_1987"](MAN).as_volumes_clearances()
    assert vc["V1"] == pytest.approx(7.77) and vc["Cl1"] == pytest.approx(0.356)


def test_maitre_covariates_sex_and_age():
    f = KERNELS["alfentanil_maitre_1987"](dict(MAN, sex="F"))
    assert f.V1 == pytest.approx(7.77 * 1.15)
    old = KERNELS["alfentanil_maitre_1987"](dict(MAN, age=60))
    assert old.k10 * old.V1 == pytest.approx(0.356 - 0.00269 * 20)
    assert old.k31 == pytest.approx(0.0126 - 0.000113 * 20)
    with pytest.raises(ValueError, match="sex"):
        KERNELS["alfentanil_maitre_1987"](dict(MAN, sex=None))


def test_scott_matches_stanpump_volume_form():
    vc = KERNELS["alfentanil_scott_1985"](MAN).as_volumes_clearances()
    for k, ref in [("V1", 2.1853), ("V2", 6.698864), ("V3", 14.52582), ("Cl1", 0.1988623),
                   ("Cl2", 1.433557), ("Cl3", 0.2469389)]:
        assert vc[k] == pytest.approx(ref, rel=2e-3), k


def test_kernels_reproduce_record_values(ds):
    for mid in (MAITRE, SCOTT):
        mp = KERNELS[ds[mid].kernel_function](MAN)
        vc = mp.as_volumes_clearances()
        for p in ds[mid].parameters:
            got = vc.get(p.symbol, getattr(mp, p.symbol, None))
            assert got == pytest.approx(p.central, rel=1e-3), (mid, p.symbol)


def test_divergence_direction_matches_the_prospective_bias(ds):
    # Maitre's measured concentrations exceeded target (MDPE +53%): it predicts LOWER
    # than the well-calibrated Scott set for the same infusion.
    t = np.linspace(0, 60, 121)
    sched = DEFAULT_SCHEDULES["alfentanil"]
    m = hypnos.simulate(ds, MAITRE, patient=MAN, schedule=sched, t=t)
    s = hypnos.simulate(ds, SCOTT, patient=MAN, schedule=sched, t=t)
    assert np.all(s.cp[t >= 10] > m.cp[t >= 10])
    mdpe = {p["metric"]: p["value"] for p in ds[MAITRE].raw["predictive_performance"]}
    assert mdpe["MDPE"] == 53.0


def test_only_cl_bsv_is_curated_for_maitre(ds):
    # the 33% Vc CV is deliberately not sampled (micro-constant eta on Vc != VC-space eta)
    assert set(ds[MAITRE].bsv_omegas()) == {"Cl1"}


@pytest.mark.parametrize("fmt", sorted(BUILDERS))
def test_micro_rate_constant_records_export(ds, fmt):
    for mid in (MAITRE, SCOTT):
        _, text = export_model(fmt, ds[mid], ds)
        assert text
