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


# --------------------------------------------------------------------------- #
# E3 — every export carries the estimation object (v0.3 §8), round-trip exact (§9)
# --------------------------------------------------------------------------- #
import copy  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import xml.dom.minidom as minidom  # noqa: E402

from hypnos.export import export_model  # noqa: E402
from hypnos.export._estimation import (  # noqa: E402
    independence_caveat, theta_covariance, theta_estimates,
)

ESTIMATION_FORMATS = ["nonmem", "pharmml", "sbml", "tci_json", "rxode2", "pumas"]


def _curated(ds):
    """symbol -> (central, se) straight from the record — the round-trip reference."""
    return {p.symbol: (p.central, p.estimation.se) for p in ds[ELEVELD].parameters
            if p.estimation and p.estimation.se is not None}


def _correlated_eleveld(ds, succeeded=True):
    """Eleveld with a constructed estimate correlation (no model publishes one yet)."""
    raw = copy.deepcopy(ds[ELEVELD].raw)
    raw["uncertainty_status"] = "correlated"
    raw["estimate_covariance"] = {
        "correlations": [{"between": ["Cl1", "V1"], "correlation": -0.42}],
        "complete": False, "method": "asymptotic_covariance",
        "covariance_step_succeeded": succeeded, "tier": "C",
        "extraction": {"review_status": "unverified"}}
    return Model(raw=raw)


def test_theta_estimates_canonical_order_and_recomputed_rse():
    ds = hypnos.load()
    ests = theta_estimates(ds[ELEVELD])
    assert [e.symbol for e in ests] == ["Cl1", "V1", "Cl2", "V2", "Cl3", "V3"]   # ke0 has no SE
    for e in ests:
        assert e.rse_percent == pytest.approx(100.0 * e.se / e.central)
    assert theta_estimates(ds[MARSH]) == []          # never-synthesize


def test_theta_estimates_log_scale_converts_to_natural():
    m = _model(est={"se": 0.06, "scale": "log", "method": "bootstrap", "tier": "B",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    (e,) = theta_estimates(m)
    assert e.se == 0.06 and e.se_natural == pytest.approx(1.79 * 0.06)
    assert e.rse_percent == pytest.approx(6.0)


def test_theta_covariance_marginal_is_diagonal_with_recorded_caveat():
    ds = hypnos.load()
    syms, cov = theta_covariance(ds[ELEVELD])
    cur = _curated(ds)
    for i, s in enumerate(syms):
        assert cov[i][i] == pytest.approx(cur[s][1] ** 2)
        assert all(cov[i][j] == 0.0 for j in range(len(syms)) if j != i)
    assert "independence assumption" in independence_caveat(ds[ELEVELD])
    assert theta_covariance(ds[MARSH]) is None


def test_theta_covariance_carries_published_correlation_only():
    ds = hypnos.load()
    m = _correlated_eleveld(ds)
    syms, cov = theta_covariance(m)
    cur = _curated(ds)
    i, j = syms.index("Cl1"), syms.index("V1")
    assert cov[i][j] == cov[j][i] == pytest.approx(-0.42 * cur["Cl1"][1] * cur["V1"][1])
    off = [cov[a][b] for a in range(len(syms)) for b in range(a) if {a, b} != {i, j}]
    assert all(v == 0.0 for v in off)                # unpublished pairs are never invented
    assert independence_caveat(m) is None


@pytest.mark.parametrize("fmt", ESTIMATION_FORMATS)
def test_export_carries_estimation_status_and_tier(fmt):
    ds = hypnos.load()
    _, text = export_model(fmt, ds[ELEVELD], ds)
    if fmt == "tci_json":
        est = json.loads(text)["estimation"]
        assert (est["uncertainty_status"], est["estimation_band_tier"]) == ("marginal", "A")
    else:
        assert re.search(r"uncertaintyStatus\W+marginal", text)
        assert re.search(r"estimationBandTier\W+A\b", text)


@pytest.mark.parametrize("fmt", ESTIMATION_FORMATS)
def test_export_never_synthesizes_estimation(fmt):
    ds = hypnos.load()
    _, text = export_model(fmt, ds[MARSH], ds)
    if fmt == "tci_json":
        assert json.loads(text)["estimation"]["uncertainty_status"] == "none"
    else:
        for marker in ("uncertaintyStatus", "theta_cov", "StandardOutput",
                       "parameterStandardError"):
            assert marker not in text


def test_tci_json_estimation_round_trips_verbatim():
    ds = hypnos.load()
    theta = json.loads(export_model("tci_json", ds[ELEVELD], ds)[1])["estimation"]["theta"]
    raw = {p["symbol"]: p for p in ds[ELEVELD].raw["parameters"]}
    assert set(theta) == set(_curated(ds))
    for sym, block in theta.items():
        assert block["estimation_uncertainty"] == raw[sym]["estimation_uncertainty"]
        assert block["estimate"] == raw[sym]["value"]["central"]


def test_pharmml_standard_output_round_trips():
    ds = hypnos.load()
    dom = minidom.parseString(export_model("pharmml", ds[ELEVELD], ds)[1])
    cur = _curated(ds)
    mle = {e.getAttribute("symbol"): float(e.getAttribute("value"))
           for e in dom.getElementsByTagName("Estimate")}
    se = {e.getAttribute("symbol"): float(e.getAttribute("value"))
          for e in dom.getElementsByTagName("SE")}
    rse = {e.getAttribute("symbol"): float(e.getAttribute("percent"))
           for e in dom.getElementsByTagName("RSE")}
    assert set(se) == set(cur)
    for sym, (central, s) in cur.items():
        assert mle[sym] == pytest.approx(central) and se[sym] == pytest.approx(s)
        assert rse[sym] == pytest.approx(100.0 * s / central, rel=1e-5)
    # the estimation object is a sibling of — never nested in — the between-subject model
    so = dom.getElementsByTagName("StandardOutput")[0]
    assert so.parentNode is dom.getElementsByTagName("VariabilityModel")[0].parentNode
    rows = dom.getElementsByTagName("Row")
    assert float(rows[-1].firstChild.data.split()[-1]) == pytest.approx(cur["V3"][1] ** 2)


def test_sbml_estimation_rdf_round_trips_beside_variability():
    ds = hypnos.load()
    text = export_model("sbml", ds[ELEVELD], ds)[1]
    dom = minidom.parseString(text)
    bag = dom.getElementsByTagName("hypnos:parameterStandardError")[0]
    got = {li.getAttribute("hypnos:parameter"): (float(li.getAttribute("hypnos:estimate")),
                                                 float(li.getAttribute("hypnos:se")))
           for li in bag.getElementsByTagName("rdf:li")}
    cur = _curated(ds)
    assert set(got) == set(cur)
    for sym in cur:
        assert got[sym] == pytest.approx(cur[sym])
    # an RSE never lands in the Ω bag, and a CV never in the SE bag (Trap 1)
    omega = dom.getElementsByTagName("hypnos:betweenSubjectVariability")[0]
    assert not any(li.hasAttribute("hypnos:se") for li in omega.getElementsByTagName("rdf:li"))
    assert not any(li.hasAttribute("hypnos:omega2") for li in bag.getElementsByTagName("rdf:li"))
    assert "sees the point estimate only" in text


def test_nonmem_cov_companion_round_trips():
    ds = hypnos.load()
    text = export_model("nonmem", ds[ELEVELD], ds)[1]
    cur = _curated(ds)
    got = {m.group(1): (float(m.group(3)), float(m.group(2))) for m in re.finditer(
        r"^;\s+(\w+): SE (\S+) \(scale natural\) on curated estimate (\S+),", text, re.M)}
    assert set(got) == set(cur)
    for sym in cur:
        assert got[sym] == pytest.approx(cur[sym], rel=1e-5)
    assert ".cov-style lower triangle" in text
    # the block is comments only — it must not change the runnable control stream
    est_lines = [ln for ln in text.splitlines() if "estimat" in ln.lower() and "$EST" not in ln]
    assert est_lines and all(ln.startswith(";") for ln in est_lines)


def test_rxode2_theta_cov_round_trips_and_is_not_omega():
    ds = hypnos.load()
    text = export_model("rxode2", ds[ELEVELD], ds)[1]
    cur = _curated(ds)
    block = text[text.index("_theta_cov <- lotri({"):]
    got = {m.group(1): float(m.group(2))
           for m in re.finditer(r"^\s+theta\.(\w+) ~ (\S+)$", block, re.M)}
    assert set(got) == set(cur)
    for sym, (_, s) in cur.items():
        assert got[sym] == pytest.approx(s ** 2, rel=1e-9)
    assert "NOT Omega" in text
    # the BSV matrix keeps its own name and its own (η) rows
    omega = text[text.index("_omega <- lotri({"):text.index("# ── estimation uncertainty")]
    assert "eta.V1 ~ 0.61" in omega and "theta." not in omega


def test_pumas_theta_cov_round_trips():
    ds = hypnos.load()
    text = export_model("pumas", ds[ELEVELD], ds)[1]
    cur = _curated(ds)
    block = text[text.index("_theta_cov = ("):]
    names = re.search(r"names = \[(.*?)\]", block).group(1).replace(":", "").split(", ")
    rows = [[float(x) for x in ln.strip().rstrip(";").split()]
            for ln in block[block.index("cov = [") + 7:block.index("    ],")].strip().splitlines()]
    assert names == [e.symbol for e in theta_estimates(ds[ELEVELD])]
    for i, sym in enumerate(names):
        assert rows[i][i] == pytest.approx(cur[sym][1] ** 2, rel=1e-9)
    assert "NOT the @random" in text


@pytest.mark.parametrize("fmt", ["nonmem", "pharmml", "sbml", "rxode2", "pumas"])
def test_export_carries_published_estimate_correlation(fmt):
    ds = hypnos.load()
    text = export_model(fmt, _correlated_eleveld(ds), ds)[1]
    cur = _curated(ds)
    covar = -0.42 * cur["Cl1"][1] * cur["V1"][1]
    if fmt != "sbml":                                 # SBML carries r itself, not a matrix
        assert f"{covar:.6g}"[:9] in text             # the off-diagonal r·seᵢ·seⱼ
    assert "covariance step succeeded" in text
    assert "independence assumption" not in text
    if fmt in ("pharmml", "sbml"):
        minidom.parseString(text)
        assert re.search(r'between="Cl1 V1"\W+(hypnos:)?(correlation|coefficient)="-0.42"', text)
    if fmt == "rxode2":                               # one lower-triangular lotri block
        assert "theta.Cl1 + theta.V1 + theta.Cl2" in text


def test_export_states_a_failed_covariance_step():
    ds = hypnos.load()
    text = export_model("nonmem", _correlated_eleveld(ds, succeeded=False), ds)[1]
    assert "covariance step FAILED" in text and "MISSING, not zero" in text


# --------------------------------------------------------------------------- #
# The confidence band must scale with the PATIENT's parameters (relative SE)
# --------------------------------------------------------------------------- #
from hypnos.export.registry import instantiate  # noqa: E402
from hypnos.reference import sample_parameter_vector  # noqa: E402

_REFERENCE = dict(age=35, weight=70, height=170, sex="M")      # Eleveld's reference individual
_CHILD = dict(age=6, weight=20, height=115, sex="M")
_OBESE = dict(age=50, weight=160, height=175, sex="M")


def test_relative_ses_equal_the_curated_rse():
    ds = hypnos.load()
    m = ds[ELEVELD]
    rel = m.estimation_relative_ses()
    for p in m.parameters:
        if p.estimation and p.estimation.se is not None:
            assert rel[p.symbol] == pytest.approx(p.estimation.rse_percent / 100.0, rel=1e-3)
    assert "ke0" not in rel                         # no curated SE -> no entry (never-synthesize)
    assert ds[MARSH].estimation_relative_ses() == {}


def test_relative_ses_treat_a_log_scale_se_as_already_relative():
    m = _model(est={"se": 0.06, "scale": "log", "method": "bootstrap", "tier": "B",
                    "extraction": {"review_status": "unverified"}}, status="marginal")
    assert m.estimation_relative_ses() == {"Cl": 0.06}


@pytest.mark.parametrize("patient", [_REFERENCE, _CHILD, _OBESE])
def test_sampled_parameters_scatter_by_the_curated_rse_for_any_patient(patient):
    # a curated SE is on the reference-individual θ; a child's volumes are several times
    # smaller, so applying the ABSOLUTE SE would inflate the relative scatter (the bug)
    ds = hypnos.load()
    m = ds[ELEVELD]
    rel = m.estimation_relative_ses()
    typical = instantiate(m, patient)
    base = typical.as_volumes_clearances()
    rng = np.random.default_rng(11)
    draws = [sample_parameter_vector(typical, rel, rng, relative=True).as_volumes_clearances()
             for _ in range(4000)]
    for sym, rse in rel.items():
        got = np.std([d[sym] for d in draws]) / base[sym]
        assert got == pytest.approx(rse, rel=0.06), f"{sym} at {patient['weight']} kg"
    assert all(d["ke0"] == pytest.approx(base["ke0"]) for d in draws[:50])   # no SE -> held fixed


def test_confidence_band_relative_width_does_not_depend_on_patient_size():
    ds = hypnos.load()

    def width(patient):
        r = hypnos.simulate(ds, ELEVELD, patient=patient, t=np.linspace(0, 30, 121),
                            schedule=[("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")],
                            bands=["confidence"], seed=7)
        q = r.cp_confidence_quantiles
        return float((q[95][-1] - q[5][-1]) / q[50][-1])

    ref, child, obese = width(_REFERENCE), width(_CHILD), width(_OBESE)
    # the same few-percent estimation uncertainty for every patient: before the fix the
    # child's band was ~2.5x the reference's and the obese adult's ~0.6x
    assert child == pytest.approx(ref, rel=0.25) and obese == pytest.approx(ref, rel=0.25)
    assert 0.02 < ref < 0.08


def test_absolute_mode_is_unchanged_for_callers_passing_natural_scale_ses():
    ds = hypnos.load()
    typical = instantiate(ds[ELEVELD], _REFERENCE)
    a = sample_parameter_vector(typical, {"V1": 0.5}, np.random.default_rng(3))
    b = sample_parameter_vector(typical, {"V1": 0.5}, np.random.default_rng(3), relative=False)
    assert a.V1 == b.V1 and a.V1 != typical.V1


# --------------------------------------------------------------------------- #
# v0.3 §2.2 / Trap 5 — second-order SEs on Ω and Σ (variance vs SD in the key name)
# --------------------------------------------------------------------------- #
from hypnos.verification import _checklist_for  # noqa: E402


def _second_order_eleveld(ds, *, omega2_se=0.08, sigma=None):
    raw = copy.deepcopy(ds[ELEVELD].raw)
    v1 = next(p for p in raw["parameters"] if p["symbol"] == "V1")
    v1["variability"]["bsv"]["omega2_se"] = omega2_se
    raw["residual_error"]["log"].update(sigma if sigma is not None else {"sd_se": 0.004})
    return Model(raw=raw)


def test_schema_accepts_second_order_ses_and_rejects_unknown_scales():
    import jsonschema
    from pathlib import Path
    ds = hypnos.load()
    schema = json.loads((Path(__file__).resolve().parents[2] / "dataset" / "schema"
                         / "model.schema.json").read_text())
    ok = _second_order_eleveld(ds)
    jsonschema.validate(ok.raw, schema)
    bad = copy.deepcopy(ok.raw)
    bad["residual_error"]["log"]["sigma_se"] = 0.004      # scale must be in the key name
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_omega2_rse_is_on_the_variance():
    ds = hypnos.load()
    m = _second_order_eleveld(ds, omega2_se=0.061)
    v = m.param("V1").variability
    assert v.omega2_rse_percent == pytest.approx(10.0)            # 0.061 / 0.61
    assert m.param("V2").variability.omega2_rse_percent is None   # not curated -> no number
    assert m.residual_error.second_order_ses() == [("log", "sd", 0.004, 0.191)]


def test_clean_second_order_block_validates_and_real_dataset_is_untouched():
    ds = hypnos.load()
    assert _check_estimation(_second_order_eleveld(ds), CITES) == []
    assert all(p.variability is None or p.variability.omega2_se is None
               for m in ds for p in m.parameters)                 # none curated: an honest gap


def test_trap5_se_on_a_scale_the_value_is_not_curated_on():
    ds = hypnos.load()
    # Eleveld curates the log residual as an SD; an SE filed on the VARIANCE has no anchor
    probs = _check_estimation(_second_order_eleveld(ds, sigma={"variance_se": 0.002}), CITES)
    assert any("variance-vs-SD mismatch (Trap 5)" in p for p in probs)
    probs = _check_estimation(_second_order_eleveld(ds, omega2_se=-0.1), CITES)
    assert any("omega2_se is negative" in p for p in probs)


def test_verify_checklist_asks_variance_or_sd():
    ds = hypnos.load()
    def trap5(model):
        return [it for it in _checklist_for(model, ds)
                if it.group == "estimation" and "Trap 5" in it.label]
    items = trap5(_second_order_eleveld(ds))
    labels = " | ".join(it.label for it in items)
    assert "V1 SE on ω² — on the VARIANCE, not the SD" in labels
    assert "residual log SE — on the SD" in labels
    assert any("RSE 13.1% on the variance" in it.value for it in items)
    assert trap5(ds[ELEVELD]) == []


def test_exports_carry_omega2_se_only_where_curated():
    ds = hypnos.load()
    m = _second_order_eleveld(ds)
    assert json.loads(export_model("tci_json", m, ds)[1])["variability"]["bsv"]["V1"]["omega2_se"] == 0.08
    sbml_text = export_model("sbml", m, ds)[1]
    minidom.parseString(sbml_text)
    assert re.search(r'hypnos:parameter="V1" hypnos:omega2="0\.61" hypnos:cvPercent="[\d.]+" '
                     r'hypnos:omega2Se="0\.08"', sbml_text)
    assert sbml_text.count("hypnos:omega2Se") == 1                 # V2..ke0 curate none
    pharmml_text = export_model("pharmml", m, ds)[1]
    minidom.parseString(pharmml_text)
    assert pharmml_text.count('varianceSe="0.08"') == 1
    assert "SE on this variance 0.08" in export_model("nonmem", m, ds)[1]
    # the residual SE rides verbatim in the lossless format
    assert json.loads(export_model("tci_json", m, ds)[1])["variability"]["residual_error"]["log"]["sd_se"] == 0.004
    for fmt in ("sbml", "pharmml", "nonmem", "tci_json"):
        plain = export_model(fmt, ds[ELEVELD], ds)[1]            # the real record curates none
        assert "omega2_se" not in plain and "omega2Se" not in plain and "varianceSe" not in plain
