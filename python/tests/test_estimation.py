"""The estimation-uncertainty layer — v0.3 E0 (vocabulary + traps) and the
reducible/irreducible decomposition (v0.3 §7).

v0.3 makes *estimation* uncertainty (the SE/RSE on the typical value — reducible,
shrinks with more data) first-class and **distinct** from between-subject
variability (the BSV CV — irreducible). The separation is enforced structurally
(its own block beside `variability`) plus the numeric traps a machine can catch
(scale, RSE↔SE, CI↔SE). Eleveld 2018 now carries curated estimation values (sourced from its open-access
Table 2 99% profile-likelihood CIs); the rest await human PDF transcription — so these
tests exercise the machinery on constructed records, the same way the schema closes the conflation
by construction.
"""
import numpy as np
import pytest

import hypnos
from hypnos.models import EstimationUncertainty, Model
from hypnos.validate import _check_estimation


# --------------------------------------------------------------------------- #
# EstimationUncertainty.rse_from_se — the consistency reference (Trap 2)
# --------------------------------------------------------------------------- #
def test_rse_from_se_natural_scale():
    e = EstimationUncertainty(se=0.107, scale="natural")
    assert e.rse_from_se(1.79) == pytest.approx(5.98, abs=0.02)


def test_rse_from_se_log_scale():
    # on the log scale, RSE% ≈ 100·SE_log (a different relation — the whole point of Trap 2)
    e = EstimationUncertainty(se=0.06, scale="log")
    assert e.rse_from_se(1.79) == pytest.approx(6.0)


def test_rse_from_se_none_without_se():
    assert EstimationUncertainty(scale="natural").rse_from_se(1.79) is None


# --------------------------------------------------------------------------- #
# parsing: estimation lives BESIDE variability, never inside it
# --------------------------------------------------------------------------- #
def _model(est=None, cov=None, status=None, central=1.79):
    p = {"symbol": "Cl", "value": {"central": central, "units": "L/min"},
         "tier": "A", "extraction": {"review_status": "unverified"}}
    if est is not None:
        p["estimation_uncertainty"] = est
    raw = {"id": "hypnotics_iv.propofol.x", "subsystem": "hypnotics_iv", "purpose": "pk",
           "tier": "A", "primary_citation": "eleveld-2018-propofol",
           "structure": {"compartments": 3, "parameterization": "volumes_clearances"},
           "parameters": [p], "extraction": {"review_status": "unverified"}}
    if cov is not None:
        raw["estimate_covariance"] = cov
    if status is not None:
        raw["uncertainty_status"] = status
    return Model(raw=raw)


def test_parses_estimation_beside_variability():
    m = _model(est={"se": 0.107, "scale": "natural", "rse_percent": 6.0,
                    "method": "asymptotic_covariance", "tier": "B",
                    "extraction": {"review_status": "unverified"}})
    e = m.param("Cl").estimation
    assert e is not None and e.se == 0.107 and e.scale == "natural"
    assert e.method == "asymptotic_covariance"
    # the BSV slot is independent and empty here — the two are never conflated
    assert m.param("Cl").variability is None


def test_estimation_band_tier_and_status():
    m = _model(est={"se": 0.107, "scale": "natural", "method": "asymptotic_covariance",
                    "tier": "C", "extraction": {"review_status": "unverified"}},
               status="marginal")
    assert m.has_published_estimation is True
    assert m.uncertainty_status == "marginal"
    assert m.estimation_tier == "C"
    assert m.estimation_band_tier == "C"          # worst of structural A and estimation C


def test_no_estimation_default():
    m = _model()
    assert m.has_published_estimation is False
    assert m.uncertainty_status == "none"
    assert m.estimation_band_tier is None          # never-synthesize: no band


# --------------------------------------------------------------------------- #
# validate traps
# --------------------------------------------------------------------------- #
CITES = {"eleveld-2018-propofol"}


def test_clean_estimation_passes():
    m = _model(est={"se": 0.107, "scale": "natural", "rse_percent": 5.98,
                    "ci95": {"low": 1.58, "high": 2.00},
                    "method": "asymptotic_covariance", "tier": "B",
                    "primary_citation": "eleveld-2018-propofol",
                    "extraction": {"review_status": "unverified"}},
               status="marginal")
    assert _check_estimation(m, CITES) == []


def test_trap2_missing_scale():
    m = _model(est={"se": 0.107, "method": "asymptotic_covariance", "tier": "B",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    probs = _check_estimation(m, CITES)
    assert any("scale" in p for p in probs)


def test_trap2_rse_disagrees_with_se():
    m = _model(est={"se": 0.107, "scale": "natural", "rse_percent": 25.0,  # should be ~6
                    "method": "asymptotic_covariance", "tier": "B",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    probs = _check_estimation(m, CITES)
    assert any("rse_percent" in p and "disagrees" in p for p in probs)


def test_trap3_ci_inconsistent_with_se():
    m = _model(est={"se": 0.107, "scale": "natural", "ci95": {"low": 1.0, "high": 1.1},
                    "method": "asymptotic_covariance", "tier": "B",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    probs = _check_estimation(m, CITES)
    assert any("ci95 width" in p for p in probs)


def test_estimation_citation_must_resolve():
    m = _model(est={"se": 0.107, "scale": "natural", "method": "asymptotic_covariance",
                    "tier": "B", "primary_citation": "nope-not-real",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    probs = _check_estimation(m, CITES)
    assert any("estimation_uncertainty cites unknown" in p for p in probs)


def test_uncertainty_status_none_with_estimation_flagged():
    m = _model(est={"se": 0.107, "scale": "natural", "method": "asymptotic_covariance",
                    "tier": "B", "extraction": {"review_status": "unverified"}},
               status="none")
    probs = _check_estimation(m, CITES)
    assert any("uncertainty_status 'none'" in p for p in probs)


def test_uncertainty_status_correlated_requires_covariance():
    m = _model(est={"se": 0.107, "scale": "natural", "method": "asymptotic_covariance",
                    "tier": "B", "extraction": {"review_status": "unverified"}},
               status="correlated")
    probs = _check_estimation(m, CITES)
    assert any("requires an estimate_covariance" in p for p in probs)


def test_estimate_covariance_parses_and_validates():
    cov = {"correlations": [{"between": ["Cl", "V1"], "correlation": -0.42}],
           "complete": False, "method": "asymptotic_covariance",
           "covariance_step_succeeded": True, "tier": "C",
           "primary_citation": "eleveld-2018-propofol",
           "extraction": {"review_status": "unverified"}}
    m = _model(est={"se": 0.107, "scale": "natural", "method": "asymptotic_covariance",
                    "tier": "B", "extraction": {"review_status": "unverified"}},
               cov=cov, status="correlated")
    assert m.estimate_covariance.correlations[0]["correlation"] == -0.42
    assert m.estimate_covariance.covariance_step_succeeded is True
    assert _check_estimation(m, CITES) == []


def test_full_dataset_validates_clean():
    # Eleveld carries sourced estimation values; the layer is active and clean.
    assert hypnos.validate_dataset() == []


# --------------------------------------------------------------------------- #
# the reducible/irreducible decomposition (v0.3 §7) on real compare output
# --------------------------------------------------------------------------- #
def test_reducibility_rollup_present_and_consistent():
    ds = hypnos.load()
    cmp = hypnos.compare(
        ds, drug="propofol",
        patient={"age": 72, "weight": 60, "height": 162, "sex": "F"},
        schedule=[("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")],
        t=np.linspace(0, 30, 181), bands=True, seed=7)
    d = cmp.divergence["ce"]
    assert "reducibility" in d
    r = d["reducibility"]
    vs = d["variance_share"]
    # without a confidence band the estimation component contributes 0, so reducible ==
    # structural; irreducible == bsv + residual; they sum to 1
    assert r["reducible"] == pytest.approx(vs["structural"], abs=1e-3)
    assert r["irreducible"] == pytest.approx(vs["bsv"] + vs["residual"], abs=1e-3)
    assert r["reducible"] + r["irreducible"] == pytest.approx(1.0, abs=1e-3)
    # Eleveld now carries curated estimation uncertainty (E0/E1) — stated honestly
    assert r["estimation_curated"] is True
    assert "estimation" not in vs                       # not folded in without a confidence band


# --------------------------------------------------------------------------- #
# E1 — estimation confidence bands (v0.3)
# --------------------------------------------------------------------------- #

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"   # carries curated estimation SEs
MARSH = "hypnotics_iv.propofol.marsh_1991"       # no estimation uncertainty


def _sim(ds, mid, **kw):
    return hypnos.simulate(ds, mid, patient=dict(age=40, weight=75, height=175, sex="M"),
                           schedule=[("bolus", 0.0, "2 mg/kg")], t=np.linspace(0, 30, 200),
                           bands=["confidence"], seed=7, **kw)


def test_confidence_band_drawn_when_estimation_curated():
    ds = hypnos.load()
    res = _sim(ds, ELEVELD)
    assert res.cp_confidence_quantiles is not None
    assert res.confidence_band_tier == "A"
    q = res.cp_confidence_quantiles
    i = int(np.argmax(q[50]))
    assert q[5][i] < q[50][i] < q[95][i]              # a real band around the median


def test_confidence_band_is_narrower_than_prediction_band():
    # estimation uncertainty (reducible, how well the mean is pinned) is far tighter than
    # between-subject spread (irreducible) — the load-bearing v0.3 distinction.
    ds = hypnos.load()
    res = hypnos.simulate(ds, ELEVELD, patient=dict(age=40, weight=75, height=175, sex="M"),
                          schedule=[("bolus", 0.0, "2 mg/kg")], t=np.linspace(0, 30, 200),
                          bands=["confidence", "prediction"], seed=7)
    i = int(np.argmax(res.cp))
    conf_w = res.cp_confidence_quantiles[95][i] - res.cp_confidence_quantiles[5][i]
    pred_w = res.cp_quantiles[95][i] - res.cp_quantiles[5][i]
    assert conf_w < pred_w


def test_confidence_band_never_synthesizes():
    ds = hypnos.load()
    res = _sim(ds, MARSH)
    assert res.cp_confidence_quantiles is None        # no curated SE -> no band
    assert any("never-synthesize" in w for w in res.warnings)


def test_confidence_band_is_reproducible():
    ds = hypnos.load()
    a = _sim(ds, ELEVELD).cp_confidence_quantiles[95]
    b = _sim(ds, ELEVELD).cp_confidence_quantiles[95]
    np.testing.assert_array_equal(a, b)               # seeded -> byte-identical


def test_four_way_decomposition_with_confidence_band():
    # E2: requesting the confidence band folds the estimation component into the decomposition.
    ds = hypnos.load()
    cmp = hypnos.compare(
        ds, drug="propofol", patient={"age": 50, "weight": 80, "height": 178, "sex": "M"},
        schedule=[("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")],
        t=np.linspace(0, 30, 181), bands=["prediction", "confidence"], seed=7)
    d = cmp.divergence["ce"]
    vs = d["variance_share"]
    assert "estimation" in vs and vs["estimation"] >= 0.0
    # estimation is reducible (more data) -> folded into the reducible share
    r = d["reducibility"]
    assert r["reducible"] == pytest.approx(vs["structural"] + vs["estimation"], abs=1e-3)
    assert "estimation (more data per model)" in r["note"]
