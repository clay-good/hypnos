"""The estimation-uncertainty layer — v0.3 E0 (vocabulary + traps) and the
reducible/irreducible decomposition (v0.3 §7).

v0.3 makes *estimation* uncertainty (the SE/RSE on the typical value — reducible,
shrinks with more data) first-class and **distinct** from between-subject
variability (the BSV CV — irreducible). The separation is enforced structurally
(its own block beside `variability`) plus the numeric traps a machine can catch
(scale, RSE↔SE, CI↔SE). No model carries curated estimation values yet — they
await human PDF transcription of each RSE table — so these tests exercise the
machinery on constructed records, the same way the schema closes the conflation
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
    # no model carries estimation values yet; the layer is dormant and clean.
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
    # reducible == structural; irreducible == bsv + residual; they sum to 1
    assert r["reducible"] == pytest.approx(vs["structural"], abs=1e-3)
    assert r["irreducible"] == pytest.approx(vs["bsv"] + vs["residual"], abs=1e-3)
    assert r["reducible"] + r["irreducible"] == pytest.approx(1.0, abs=1e-3)
    # no estimation curated for any model yet — stated honestly, never silently
    assert r["estimation_curated"] is False
