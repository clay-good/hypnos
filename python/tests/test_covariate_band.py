"""Tests for the v0.7 C2 covariate-value band + the fifth variance component.

The covariate band propagates a caller-supplied covariate-VALUE distribution
(e.g. an estimated weight) through forward simulation, seeded and reproducible.
The never-invent rule is load-bearing: absent a distribution, no value band is
drawn. The variance decomposition gains a fifth (covariate) component, split into
equation-choice and value-uncertainty — and the legacy bands=True path is unchanged.
"""
import numpy as np
import pytest

import hypnos
from hypnos import covariates as cov

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
ELDERLY = {"age": 72, "weight": 60, "height": 162, "sex": "F"}


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def _t():
    return np.linspace(0, 30, 121)


# --------------------------------------------------------------------------- #
# point-patient collapse — backward compatibility
# --------------------------------------------------------------------------- #
def test_point_patient_collapse_matches_scalar(ds):
    """A distribution-valued covariate collapsed to its mean must reproduce the
    scalar curve exactly (the deterministic path is untouched)."""
    t = _t()
    scalar = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": 80, "height": 175, "sex": "M"},
                             schedule=SCHED, t=t)
    dist = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": {"mean": 80}, "height": 175, "sex": "M"},
                           schedule=SCHED, t=t)
    np.testing.assert_allclose(dist.ce, scalar.ce, rtol=0, atol=0)


def test_scalar_patient_unaffected_by_band_machinery(ds):
    t = _t()
    a = hypnos.simulate(ds, ELEVELD, patient=ELDERLY, schedule=SCHED, t=t)
    b = hypnos.simulate(ds, ELEVELD, patient=ELDERLY, schedule=SCHED, t=t, bands=True, seed=1, samples=100)
    np.testing.assert_allclose(a.ce, b.ce, rtol=0, atol=0)   # the median line is the point estimate


# --------------------------------------------------------------------------- #
# the covariate-value band
# --------------------------------------------------------------------------- #
def test_value_band_drawn_only_with_a_distribution(ds):
    t = _t()
    # scalar weight -> no value band (never-invent), but equation-choice variance present
    scalar = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": 130, "height": 170, "sex": "M"},
                             schedule=SCHED, t=t, bands=["covariate"], seed=7, samples=200)
    assert scalar.ce_covariate_band is None
    assert scalar.ce_cov_equation_var is not None
    assert any("no covariate-value distribution" in w for w in scalar.warnings)
    # distribution -> value band drawn
    dist = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": {"mean": 130, "sd": 12}, "height": 170, "sex": "M"},
                           schedule=SCHED, t=t, bands=["covariate"], seed=7, samples=200)
    assert dist.ce_covariate_band is not None
    assert dist.covariate_band_tier is not None


def test_value_band_is_reproducible(ds):
    t = _t()
    p = {"age": 50, "weight": {"mean": 130, "sd": 12}, "height": 170, "sex": "M"}
    a = hypnos.simulate(ds, ELEVELD, patient=p, schedule=SCHED, t=t, bands=["covariate"], seed=7, samples=300)
    b = hypnos.simulate(ds, ELEVELD, patient=p, schedule=SCHED, t=t, bands=["covariate"], seed=7, samples=300)
    for q in (5, 50, 95):
        np.testing.assert_allclose(a.ce_covariate_band[q], b.ce_covariate_band[q], rtol=0, atol=0)


def test_value_band_widens_with_uncertainty(ds):
    t = _t()
    def width(sd):
        r = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": {"mean": 130, "sd": sd}, "height": 170, "sex": "M"},
                            schedule=SCHED, t=t, bands=["covariate"], seed=7, samples=600)
        i = int(np.argmax(r.ce))
        return r.ce_covariate_band[95][i] - r.ce_covariate_band[5][i]
    assert width(20) > width(5)             # more weight uncertainty -> wider band


def test_band_requires_seed(ds):
    with pytest.raises(ValueError, match="seed"):
        hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": {"mean": 130, "sd": 12}, "height": 170, "sex": "M"},
                        schedule=SCHED, t=_t(), bands=["covariate"])


def test_effect_covariate_band_when_pd_attached(ds):
    t = _t()
    r = hypnos.simulate(ds, ELEVELD, patient={"age": 50, "weight": {"mean": 130, "sd": 12}, "height": 170, "sex": "M"},
                        schedule=SCHED, t=t, pd_model="pd_effect.propofol.eleveld_bis",
                        bands=["covariate"], seed=7, samples=200)
    assert r.effect_covariate_band is not None


# --------------------------------------------------------------------------- #
# the fifth variance component in compare()
# --------------------------------------------------------------------------- #
def test_compare_adds_covariate_component_and_split(ds):
    t = _t()
    p = {"age": 50, "weight": {"mean": 130, "sd": 13}, "height": 170, "sex": "M"}
    cmp = hypnos.compare(ds, drug="propofol", patient=p, schedule=SCHED, t=t,
                         bands=["prediction", "covariate"], seed=7, samples=400)
    d = cmp.divergence["ce"]
    vs = d["variance_share"]
    assert "covariate" in vs
    # the five shares sum to 1
    s = vs["structural"] + vs["bsv"] + vs["residual"] + vs["covariate"]
    assert s == pytest.approx(1.0, abs=1e-3)
    split = d["covariate_split"]
    assert split["equation_choice"] + split["value_uncertainty"] == pytest.approx(vs["covariate"], abs=1e-3)
    # reducible includes the covariate component
    r = d["reducibility"]
    assert r["reducible"] == pytest.approx(vs["structural"] + vs["covariate"], abs=1e-3)
    assert r["reducible"] + r["irreducible"] == pytest.approx(1.0, abs=1e-3)
    assert d.get("covariate_band_tier") is not None


def test_compare_legacy_bands_true_unchanged(ds):
    """The v0.2 path (bands=True) keeps its exact 3-way decomposition — no covariate key."""
    t = _t()
    cmp = hypnos.compare(ds, drug="propofol", patient=ELDERLY, schedule=SCHED, t=t,
                         bands=True, seed=7, samples=400)
    vs = cmp.divergence["ce"]["variance_share"]
    assert "covariate" not in vs
    assert vs["structural"] + vs["bsv"] + vs["residual"] == pytest.approx(1.0, abs=1e-3)
    assert "covariate_split" not in cmp.divergence["ce"]


# --------------------------------------------------------------------------- #
# sample_covariate_vector — the seeded, never-invent primitive
# --------------------------------------------------------------------------- #
def test_sample_covariate_vector_recovers_marginals():
    rng = np.random.default_rng(3)
    p = {"age": 50, "weight": {"mean": 70, "sd": 6}, "height": 170, "sex": "M"}
    draws = np.array([cov.sample_covariate_vector(p, rng)["weight"] for _ in range(8000)])
    assert draws.mean() == pytest.approx(70, abs=0.5)
    assert draws.std() == pytest.approx(6, abs=0.5)


def test_sample_holds_scalars_fixed():
    rng = np.random.default_rng(1)
    p = {"age": 50, "weight": {"mean": 70, "sd": 6}, "height": 170, "sex": "M"}
    pv = cov.sample_covariate_vector(p, rng)
    assert pv["age"] == 50 and pv["height"] == 170 and pv["sex"] == "M"  # untouched


def test_cli_simulate_covariate_band(ds, capsys):
    from hypnos.cli import main
    rc = main(["simulate", ELEVELD, "--age", "50", "--weight", "130", "--height", "170",
               "--sex", "M", "--covariate-band", "--weight-sd", "12", "--samples", "200", "--seed", "7"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "covariate band (value uncertainty)" in out
