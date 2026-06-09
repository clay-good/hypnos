"""Tests for the v0.2 population-variability layer.

Covers V0 (schema + curation + validate consistency), V1 (sample_individual +
seeded Monte-Carlo bands + the never-synthesize rule), V2 (separation index +
variance decomposition), and V3 (NONMEM $OMEGA/$SIGMA + tci_json passthrough).
"""
import json
import math
import re

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
# V1 / §14 — PD effect band: PK BSV propagated through the (fixed) PD link
# --------------------------------------------------------------------------- #
ELEVELD_BIS = "pd_effect.propofol.eleveld_bis"


def test_effect_band_present_ordered_and_reproducible(ds):
    t = np.linspace(0, 30, 80)
    r1 = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t,
                         pd_model=ELEVELD_BIS, bands=True, seed=9, samples=400)
    r2 = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t,
                         pd_model=ELEVELD_BIS, bands=True, seed=9, samples=400)
    assert r1.effect_quantiles is not None and set(r1.effect_quantiles) == {5, 50, 95}
    lo, hi = r1.band_percentile
    assert np.array_equal(r1.effect_quantiles[lo], r2.effect_quantiles[lo])  # seeded
    assert np.all(r1.effect_quantiles[lo] <= r1.effect_quantiles[50] + 1e-9)
    assert np.all(r1.effect_quantiles[50] <= r1.effect_quantiles[hi] + 1e-9)
    # BIS stays within a sane scale and the band has real width somewhere
    assert float(r1.effect_quantiles[lo].min()) >= 0.0
    assert float(r1.effect_quantiles[hi].max()) <= 100.0
    assert float((r1.effect_quantiles[hi] - r1.effect_quantiles[lo]).max()) > 1.0


def test_effect_band_labeled_lower_bound(ds):
    t = np.linspace(0, 30, 60)
    r = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t,
                        pd_model=ELEVELD_BIS, bands=True, seed=3, samples=200)
    # PD-parameter BSV is uncurated → the effect band must declare itself a lower bound
    assert any("LOWER BOUND on true effect" in w for w in r.warnings)


def test_no_effect_band_without_pd_model(ds):
    t = np.linspace(0, 30, 40)
    r = hypnos.simulate(ds, ELEVELD, patient=PATIENT, schedule=SCHED, t=t, bands=True, seed=3, samples=100)
    assert r.effect_quantiles is None  # no PD link composed → no effect band


def test_no_effect_band_for_no_bsv_pk(ds):
    # Marsh publishes no BSV → never-synthesize: no PK band, hence no effect band either
    t = np.linspace(0, 30, 40)
    r = hypnos.simulate(ds, MARSH, patient=PATIENT, schedule=SCHED, t=t,
                        pd_model="pd_effect.propofol.bis_sigmoid", bands=True, seed=3, samples=50)
    assert r.effect_quantiles is None
    assert r.effect is not None  # the deterministic effect line is still drawn


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


# --------------------------------------------------------------------------- #
# V3 — random-effects projection across the pharmacometric formats
# --------------------------------------------------------------------------- #
def test_omega_diagonal_is_canonically_ordered(ds):
    from hypnos.export._variability import VC_ORDER, omega_diagonal

    syms = [s for s, _, _ in omega_diagonal(ds[ELEVELD])]
    assert syms == VC_ORDER  # all seven, clearances/volumes interleaved then ke0


def test_residual_spec_normalizes_log(ds):
    from hypnos.export._variability import residual_spec

    spec = residual_spec(ds[ELEVELD])
    assert spec.model == "log"
    assert spec.log_sd == pytest.approx(0.191)
    assert residual_spec(ds[MARSH]) is None  # Marsh curates no residual error


def test_rxode2_emits_population_block(ds):
    _, text = export_model("rxode2", ds[ELEVELD], ds)
    assert "_pop <- rxode2(" in text          # the population companion model
    assert "Cl1 <- 1.922638103 * exp(eta.Cl1)" in text
    assert "_omega <- lotri(" in text
    assert "eta.Cl1 ~ 0.265" in text
    assert "cp ~ lnorm(0.191)" in text        # log residual endpoint


def test_rxode2_population_block_round_trips_micro_constants(ds):
    # the V/Cl population block must derive the SAME typical micro-constants the
    # typical-value model emits (so a population sim collapses to the line at eta=0).
    _, text = export_model("rxode2", ds[ELEVELD], ds)
    # parse the "# hypnos.params:" k=v line (the typical micro-constants)
    kv = dict(tok.split("=") for tok in
              re.search(r"# hypnos\.params: (.+)", text).group(1).split())
    # recover Cl/V reference values from the population block and recompute k10
    cl1 = float(re.search(r"Cl1 <- ([\d.eE+-]+) \* exp", text).group(1))
    v1 = float(re.search(r"V1  <- ([\d.eE+-]+) \* exp", text).group(1))
    assert cl1 / v1 == pytest.approx(float(kv["k10"]), rel=1e-9)


def test_rxode2_no_population_block_for_marsh(ds):
    _, text = export_model("rxode2", ds[MARSH], ds)
    assert "_pop <- rxode2(" not in text
    assert "lotri(" not in text


def test_pumas_emits_random_block(ds):
    _, text = export_model("pumas", ds[ELEVELD], ds)
    assert "_pop = @model begin" in text
    assert "@random begin" in text
    assert "η_Cl1 ~ Normal(0.0, sqrt(ω²_Cl1))" in text
    assert "Cl1 = 1.922638103 * exp(η_Cl1)" in text
    assert "LogNormal(log(cp), 0.191)" in text


def test_pumas_no_random_block_for_marsh(ds):
    _, text = export_model("pumas", ds[MARSH], ds)
    assert "@random begin" not in text


def test_pharmml_emits_random_effects(ds):
    import xml.dom.minidom as minidom

    _, text = export_model("pharmml", ds[ELEVELD], ds)
    minidom.parseString(text)  # well-formed
    assert '<VariabilityModel bandTier="B" variabilityStatus="diagonal">' in text
    assert 'parameter="Cl1"' in text and 'variance="0.265"' in text
    assert '<ResidualError model="log"' in text and 'logSd="0.191"' in text


def test_pharmml_no_variability_for_marsh(ds):
    _, text = export_model("pharmml", ds[MARSH], ds)
    assert "<VariabilityModel" not in text


# --------------------------------------------------------------------------- #
# V3 — NONMEM $OMEGA BLOCK (off-diagonal Omega), exercised with a curated block
# --------------------------------------------------------------------------- #
def _eleveld_with_block(ds, between, correlation, complete=True):
    """Deep-copy the real Eleveld record and inject an omega_block (no dataset edit)."""
    import copy

    raw = copy.deepcopy(ds[ELEVELD].raw)
    raw["variability_status"] = "full"
    raw["omega_block"] = {
        "correlations": [{"between": list(between), "correlation": correlation,
                          "citation": "eleveld-2018-propofol"}],
        "complete": complete, "tier": "C",
        "extraction": {"review_status": "unverified"},
    }
    return Model(raw)


def test_contiguous_block_builds_front_anchored_covariance(ds):
    from hypnos.export._variability import contiguous_block

    m = _eleveld_with_block(ds, ("Cl1", "V1"), 0.6)
    syms, cov = contiguous_block(m)
    assert syms == ["Cl1", "V1"]
    assert cov[0][0] == pytest.approx(0.265)
    assert cov[1][1] == pytest.approx(0.610)
    assert cov[0][1] == pytest.approx(0.6 * math.sqrt(0.265 * 0.610))  # r*sqrt(om_i*om_j)


def test_contiguous_block_rejects_incomplete(ds):
    from hypnos.export._variability import contiguous_block

    assert contiguous_block(_eleveld_with_block(ds, ("Cl1", "V1"), 0.6, complete=False)) is None


def test_contiguous_block_rejects_non_front_anchored(ds):
    from hypnos.export._variability import contiguous_block

    # V1–Cl2 occupy eta positions 1–2, not anchored at eta 1 → no valid BLOCK
    assert contiguous_block(_eleveld_with_block(ds, ("V1", "Cl2"), 0.5)) is None


def test_nonmem_emits_omega_block_when_curated(ds):
    m = _eleveld_with_block(ds, ("Cl1", "V1"), 0.6)
    text = __import__("hypnos.export.nonmem", fromlist=["build"]).build(m, ds)
    assert "$OMEGA BLOCK(2)" in text
    cov = 0.6 * math.sqrt(0.265 * 0.610)
    assert f"{cov:.6g}" in text          # off-diagonal covariance emitted
    assert "$OMEGA  ; diagonal BSV" in text   # remaining five etas stay diagonal


def test_nonmem_diagonal_when_block_non_contiguous(ds):
    m = _eleveld_with_block(ds, ("V1", "Cl2"), 0.5)  # not front-anchored → fall back
    text = __import__("hypnos.export.nonmem", fromlist=["build"]).build(m, ds)
    assert "$OMEGA BLOCK(" not in text  # no BLOCK directive emitted
    assert "not emitted as a $OMEGA BLOCK" in text
