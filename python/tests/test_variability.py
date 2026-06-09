"""Tests for the v0.2 population-variability layer.

Covers V0 (schema + curation + validate consistency), V1 (sample_individual +
seeded Monte-Carlo bands + the never-synthesize rule), V2 (separation index +
variance decomposition), and V3 (NONMEM $OMEGA/$SIGMA + tci_json passthrough).
"""
import json
import math

import numpy as np
import pytest

import hypnos
from hypnos.export import export_model
from hypnos.models import Model
from hypnos.reference import (
    MicroParams,
    apply_residual,
    residual_std,
    sample_individual,
)
from hypnos.simulate import SimulationResult, _band_divergence
from hypnos.validate import _check_variability

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
MARSH = "hypnotics_iv.propofol.marsh_1991"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --------------------------------------------------------------------------- #
# V0 — curation + typed views + validate consistency
# --------------------------------------------------------------------------- #
def test_dataset_valid(ds):
    assert hypnos.validate_dataset(ds) == []


def test_eleveld_carries_diagonal_variability(ds):
    m = ds[ELEVELD]
    assert m.variability_status == "diagonal"
    assert m.has_published_variability
    omegas = m.bsv_omegas()
    assert set(omegas) == {"V1", "V2", "V3", "Cl1", "Cl2", "Cl3", "ke0"}
    assert omegas["V1"] == pytest.approx(0.610)
    assert m.residual_error.model == "log"
    assert m.residual_error.log["sd"] == pytest.approx(0.191)


def test_band_tier_is_at_or_below_point_tier(ds):
    m = ds[ELEVELD]
    # the median line keeps its A tier; the band around it is labeled lower (B)
    assert m.tier == "A"
    assert m.band_tier == "B"
    assert m.variability_tier == "B"


def test_cv_percent_recomputes_from_omega2(ds):
    for p in ds[ELEVELD].parameters:
        v = p.variability
        if v and v.omega2 is not None:
            assert v.cv_percent == pytest.approx(v.cv_from_omega2, abs=1.0)


def test_marsh_has_no_variability(ds):
    m = ds[MARSH]
    assert m.variability_status == "none"
    assert not m.has_published_variability
    assert m.band_tier is None
    assert m.bsv_omegas() == {}


def _raw_with_variability(**root):
    """Minimal model raw dict for exercising _check_variability in isolation."""
    base = {
        "id": "x.y.z", "subsystem": "x", "drug": {"name": "y"}, "purpose": "pk",
        "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
        "parameters": [{
            "symbol": "V1", "value": {"central": 1.0, "units": "L"}, "tier": "A",
            "extraction": {"review_status": "unverified"},
        }],
        "tier": "A", "extraction": {"review_status": "unverified"},
        "primary_citation": "c",
    }
    base.update(root)
    return base


def test_validate_flags_cv_omega2_mismatch():
    raw = _raw_with_variability(variability_status="diagonal")
    raw["parameters"][0]["variability"] = {
        "bsv": {"omega2": 0.25, "cv_percent": 99.0},  # wrong: sqrt(exp(.25)-1)~53%
        "tier": "B", "extraction": {"review_status": "unverified"},
    }
    problems = _check_variability(Model(raw), known_citations={"c"})
    assert any("cv_percent" in p for p in problems)


def test_validate_flags_full_without_block():
    raw = _raw_with_variability(variability_status="full")
    raw["parameters"][0]["variability"] = {
        "bsv": {"omega2": 0.25}, "tier": "B", "extraction": {"review_status": "unverified"},
    }
    problems = _check_variability(Model(raw), known_citations={"c"})
    assert any("requires an omega_block" in p for p in problems)


def test_validate_flags_none_with_curated_bsv():
    raw = _raw_with_variability(variability_status="none")
    raw["parameters"][0]["variability"] = {
        "bsv": {"omega2": 0.25}, "tier": "B", "extraction": {"review_status": "unverified"},
    }
    problems = _check_variability(Model(raw), known_citations={"c"})
    assert any("'none'" in p for p in problems)


def test_validate_flags_unknown_variability_citation():
    raw = _raw_with_variability(variability_status="diagonal")
    raw["parameters"][0]["variability"] = {
        "bsv": {"omega2": 0.25}, "tier": "B", "primary_citation": "ghost",
        "extraction": {"review_status": "unverified"},
    }
    problems = _check_variability(Model(raw), known_citations={"c"})
    assert any("ghost" in p for p in problems)


# --------------------------------------------------------------------------- #
# V1 — sampling kernel + residual helpers
# --------------------------------------------------------------------------- #
def test_sample_individual_no_omega_is_identity():
    typ = MicroParams.from_volumes_clearances(V1=4.0, Cl1=2.0, V2=20.0, Cl2=1.0, V3=200.0, Cl3=0.8, ke0=0.4)
    rng = np.random.default_rng(0)
    drawn = sample_individual(typ, {}, rng)
    assert drawn.as_volumes_clearances() == pytest.approx(typ.as_volumes_clearances())


def test_sample_individual_is_seed_reproducible():
    typ = MicroParams.from_volumes_clearances(V1=4.0, Cl1=2.0, V2=20.0, Cl2=1.0, V3=200.0, Cl3=0.8, ke0=0.4)
    omegas = {"V1": 0.1, "Cl1": 0.2}
    a = sample_individual(typ, omegas, np.random.default_rng(42)).as_volumes_clearances()
    b = sample_individual(typ, omegas, np.random.default_rng(42)).as_volumes_clearances()
    assert a == pytest.approx(b)
    # a perturbed draw differs from the typical value
    assert a["Cl1"] != pytest.approx(2.0)


def test_sample_individual_lognormal_mean_recovers_typical():
    typ = MicroParams.from_volumes_clearances(V1=4.0, Cl1=2.0, V2=20.0, Cl2=1.0, V3=200.0, Cl3=0.8, ke0=0.4)
    rng = np.random.default_rng(1)
    cls = [sample_individual(typ, {"Cl1": 0.09}, rng).as_volumes_clearances()["Cl1"] for _ in range(4000)]
    # E[exp(eta)] = exp(omega2/2); median ~ typical
    assert np.median(cls) == pytest.approx(2.0, rel=0.05)


def test_residual_std_models():
    c = np.array([0.0, 2.0, 4.0])
    assert residual_std(c, "proportional", prop_var=0.04) == pytest.approx(c * 0.2)
    assert residual_std(c, "additive", add_sd=0.5) == pytest.approx(np.full(3, 0.5))
    log = residual_std(c, "log", log_sd=0.191)
    assert log[0] == pytest.approx(0.0)
    assert log[2] == pytest.approx(4.0 * math.sqrt(math.exp(0.191 ** 2) - 1.0))


def test_apply_residual_log_is_multiplicative():
    c = np.full(5000, 3.0)
    out = apply_residual(c, "log", np.random.default_rng(3), log_sd=0.2)
    assert np.all(out > 0)
    assert np.median(out) == pytest.approx(3.0, rel=0.05)


# --------------------------------------------------------------------------- #
# V1 — bands in simulate
# --------------------------------------------------------------------------- #
PATIENT = dict(age=50, weight=77, height=177, sex="M")
SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]


def test_bands_require_seed(ds):
    t = np.linspace(0, 30, 50)
    with pytest.raises(ValueError):
        hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True)


def test_bands_are_ordered_and_reproducible(ds):
    t = np.linspace(0, 30, 60)
    r1 = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=11, samples=300)
    r2 = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=11, samples=300)
    lo, hi = r1.band_percentile
    assert np.array_equal(r1.ce_quantiles[lo], r2.ce_quantiles[lo])  # byte-reproducible
    assert np.all(r1.ce_quantiles[lo] <= r1.ce_quantiles[50] + 1e-9)
    assert np.all(r1.ce_quantiles[50] <= r1.ce_quantiles[hi] + 1e-9)
    assert r1.band_tier == "B"


def test_never_synthesize_band(ds):
    t = np.linspace(0, 30, 40)
    r = hypnos.simulate(ds, MARSH, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=1, samples=20)
    assert r.ce_quantiles is None and r.cp_quantiles is None
    assert any("never-synthesize" in w for w in r.warnings)


def test_residual_band_is_wider(ds):
    t = np.linspace(0, 30, 50)
    bsv = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=5, samples=600)
    obs = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=5, samples=600, residual=True)
    lo, hi = bsv.band_percentile
    i = int(np.argmax(bsv.ce_quantiles[50]))
    assert (obs.ce_quantiles[hi][i] - obs.ce_quantiles[lo][i]) >= (bsv.ce_quantiles[hi][i] - bsv.ce_quantiles[lo][i])
    assert obs.band_includes_residual


# --------------------------------------------------------------------------- #
# V2 — separation index + variance decomposition
# --------------------------------------------------------------------------- #
def _fake_result(model_id, median, lo_val, hi_val, *, bsv_var=1.0, resid_var=0.1, tier="B"):
    n = 11
    t = np.linspace(0, 10, n)
    arr = lambda v: np.full(n, float(v))
    r = SimulationResult(model_id=model_id, t=t, cp=arr(median), ce=arr(median), tier=tier)
    r.ce_quantiles = {5: arr(lo_val), 50: arr(median), 95: arr(hi_val)}
    r.cp_quantiles = r.ce_quantiles
    r.band_percentile = (5, 95)
    r.band_tier = tier
    r.ce_bsv_var = arr(bsv_var); r.cp_bsv_var = arr(bsv_var)
    r.ce_resid_var = arr(resid_var); r.cp_resid_var = arr(resid_var)
    return r


def test_separation_disjoint_bands():
    a = _fake_result("m.a.high", median=5.0, lo_val=4.0, hi_val=6.0)
    b = _fake_result("m.b.low", median=1.0, lo_val=0.5, hi_val=1.5)
    d = _band_divergence([a, b], "ce")
    sep = d["separation"]
    assert sep["bands_disjoint_at_tstar"] is True
    assert sep["value"] > 0
    assert sep["fraction_trajectory_disjoint"] == pytest.approx(1.0)
    assert sep["driver_high"] == "m.a.high" and sep["driver_low"] == "m.b.low"


def test_separation_overlapping_bands():
    a = _fake_result("m.a.high", median=5.0, lo_val=0.5, hi_val=6.0)
    b = _fake_result("m.b.low", median=1.0, lo_val=0.5, hi_val=1.5)
    d = _band_divergence([a, b], "ce")
    assert d["separation"]["bands_disjoint_at_tstar"] is False
    assert d["separation"]["fraction_trajectory_disjoint"] == pytest.approx(0.0)


def test_variance_share_sums_to_one():
    a = _fake_result("m.a.high", median=5.0, lo_val=4.0, hi_val=6.0, bsv_var=2.0, resid_var=0.5)
    b = _fake_result("m.b.low", median=1.0, lo_val=0.5, hi_val=1.5, bsv_var=2.0, resid_var=0.5)
    vs = _band_divergence([a, b], "ce")["variance_share"]
    assert vs["structural"] + vs["bsv"] + vs["residual"] == pytest.approx(1.0, abs=1e-3)
    assert vs["structural"] > 0  # two separated medians => real structural variance


def test_compare_bands_names_excluded(ds):
    t = np.linspace(0, 30, 50)
    cmp = hypnos.compare(ds, drug="propofol", patient=PATIENT, schedule=SCHED, t=t,
                         bands=True, seed=7, samples=200)
    excluded = {e["model_id"] for e in cmp.excluded_from_bands}
    assert MARSH in excluded
    vs = cmp.divergence["ce"]["variance_share"]
    assert 0.0 <= vs["bsv"] <= 1.0
    # only Eleveld is band-eligible => no separation index (needs >= 2)
    assert "separation" not in cmp.divergence["ce"]


# --------------------------------------------------------------------------- #
# V3 — exports carry the random-effects layer
# --------------------------------------------------------------------------- #
def test_nonmem_emits_omega_and_sigma(ds):
    _, text = export_model("nonmem", ds[ELEVELD], ds)
    assert "$OMEGA\n" in text
    assert "0.61" in text  # V1 omega2
    assert "EXP(ETA(2))" in text
    assert "EXP(EPS(1))" in text  # log residual
    assert "$SIGMA  ; residual error" in text


def test_nonmem_no_bsv_keeps_fixed(ds):
    _, text = export_model("nonmem", ds[MARSH], ds)
    assert "$OMEGA 0 FIX" in text
    assert "$SIGMA 0 FIX" in text


def test_tci_json_carries_variability(ds):
    _, text = export_model("tci_json", ds[ELEVELD], ds)
    doc = json.loads(text)
    assert doc["variability"]["variability_status"] == "diagonal"
    assert doc["variability"]["bsv"]["V1"]["omega2"] == pytest.approx(0.610)
    assert doc["variability"]["residual_error"]["model"] == "log"


def test_tci_json_no_variability_for_marsh(ds):
    _, text = export_model("tci_json", ds[MARSH], ds)
    doc = json.loads(text)
    assert doc["variability"]["variability_status"] == "none"
