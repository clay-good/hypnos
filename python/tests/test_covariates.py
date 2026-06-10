"""Tests for the v0.7 covariate-model-uncertainty layer (C0).

The covariate equations are the part of a model most often re-implemented from
memory, so they get the same discipline as everything else: the library values
must equal the reference kernels they generalize, the James LBM inversion must be
a *tested* property (peak-then-decline), and the validator must enforce every
binding resolves.
"""
import pytest

import hypnos
from hypnos import covariates as cov
from hypnos.models import Model
from hypnos.reference import ffm_al_sallami, ffm_janmahasatian, lbm_james
from hypnos.validate import _check_covariate_model

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --------------------------------------------------------------------------- #
# Library loads + values equal the reference kernels they generalize
# --------------------------------------------------------------------------- #
def test_library_loads(ds):
    assert set(ds.covariate_equations) == {"james_1976", "janmahasatian_2005", "al_sallami_2015"}


def test_evaluate_matches_reference_kernels(ds):
    p = dict(age=50, weight=80, height=175, sex="M")
    assert cov.evaluate("james_1976", p, ds=ds).value == pytest.approx(lbm_james(80, 175, "M"))
    assert cov.evaluate("janmahasatian_2005", p, ds=ds).value == pytest.approx(
        ffm_janmahasatian(80, 175, "M"))
    assert cov.evaluate("al_sallami_2015", p, ds=ds).value == pytest.approx(
        ffm_al_sallami(80, 175, 50, "M"))


def test_evaluate_sex_branches_differ(ds):
    pm = dict(age=40, weight=70, height=170, sex="M")
    pf = dict(age=40, weight=70, height=170, sex="F")
    assert cov.evaluate("james_1976", pm, ds=ds).value != cov.evaluate("james_1976", pf, ds=ds).value


def test_evaluate_accepts_a_covariate_distribution_mean(ds):
    """A caller-supplied {mean, sd} weight collapses to its mean for a point eval
    (forward-compatible with the C2 covariate band; never invents a value)."""
    point = cov.evaluate("james_1976", dict(age=50, weight=80, height=175, sex="M"), ds=ds)
    dist = cov.evaluate("james_1976", dict(age=50, weight={"mean": 80, "sd": 6}, height=175, sex="M"), ds=ds)
    assert dist.value == pytest.approx(point.value)


def test_unknown_equation_raises(ds):
    with pytest.raises(KeyError):
        cov.evaluate("devine_1974", dict(weight=70, height=170, sex="M"), ds=ds)


# --------------------------------------------------------------------------- #
# The James LBM inversion is a TESTED property of the library record
# --------------------------------------------------------------------------- #
def test_james_lbm_peaks_then_declines(ds):
    """At fixed height, James LBM rises, peaks, and then DECREASES with weight —
    the non-physical inversion (v0.7 §9: a tested property, not an assertion)."""
    height = 170.0
    weights = [w for w in range(40, 200, 2)]
    lbms = [lbm_james(w, height, "M") for w in weights]
    peak_i = max(range(len(lbms)), key=lambda i: lbms[i])
    # the peak is interior (not at the lightest or heaviest weight)
    assert 0 < peak_i < len(lbms) - 1
    # strictly rising before the peak, strictly falling after it
    assert all(lbms[i] < lbms[i + 1] for i in range(peak_i))
    assert lbms[peak_i] > lbms[-1]


def test_james_inverted_flag_and_tier_down_in_obesity(ds):
    obese = dict(age=50, weight=130, height=170, sex="M")     # BMI ~45
    r = cov.evaluate("james_1976", obese, ds=ds)
    assert r.inverted is True
    assert r.out_of_envelope is True
    assert r.tier == "D"
    assert any("INVERTED" in w for w in r.warnings)


def test_james_well_behaved_in_normal_range(ds):
    normal = dict(age=50, weight=75, height=178, sex="M")     # BMI ~24
    r = cov.evaluate("james_1976", normal, ds=ds)
    assert r.inverted is False
    assert r.out_of_envelope is False
    assert r.tier == "B"
    assert r.warnings == []


def test_ffm_equations_are_monotone_not_inverted(ds):
    """The whole point of the divergence: the FFM equations do NOT invert where
    James does, so swapping equations changes the answer (and the trustworthiness)."""
    obese = dict(age=50, weight=130, height=170, sex="M")
    for eid in ("janmahasatian_2005", "al_sallami_2015"):
        r = cov.evaluate(eid, obese, ds=ds)
        assert r.inverted is False, eid


# --------------------------------------------------------------------------- #
# Model bindings + covariate-layer tier
# --------------------------------------------------------------------------- #
def test_schnider_binds_james_for_clearance(ds):
    m = ds[SCHNIDER]
    cm = m.covariate_model
    assert cm is not None
    assert m.covariate_sensitivity_status == "declared"
    di = cm.derived_inputs[0]
    assert di.equation == "james_1976"
    assert di.quantity == "lbm"
    assert "Cl1" in di.used_for
    assert di.verbatim is True


def test_covariate_layer_tier_tracks_the_equation(ds):
    m = ds[SCHNIDER]
    assert cov.covariate_layer_tier(m, dict(age=50, weight=75, height=178, sex="M"), ds) == "B"
    # obese: James inverted -> covariate layer forced to D
    assert cov.covariate_layer_tier(m, dict(age=50, weight=130, height=170, sex="M"), ds) == "D"


def test_model_without_covariate_model_has_none(ds):
    # Marsh scales on raw weight only — no derived equation, an honest gap
    marsh = ds["hypnotics_iv.propofol.marsh_1991"]
    assert marsh.covariate_model is None
    assert marsh.covariate_sensitivity_status == "none"
    assert cov.covariate_layer_tier(marsh, dict(age=50, weight=70, height=170, sex="M"), ds) is None


# --------------------------------------------------------------------------- #
# Validator enforces the bindings (a silent substitution cannot slip in)
# --------------------------------------------------------------------------- #
def test_validator_flags_unknown_equation(ds):
    raw = dict(ds[SCHNIDER].raw)
    raw = {**raw, "covariate_model": {"derived_inputs": [
        {"quantity": "lbm", "equation": "no_such_equation", "used_for": ["Cl1"],
         "verbatim": True, "tier": "B",
         "extraction": {"review_status": "unverified"}}]}}
    problems = _check_covariate_model(Model(raw), ds, set(ds.citations))
    assert any("not in covariate_equations" in p for p in problems)


def test_validator_flags_unknown_used_for_symbol(ds):
    raw = {**dict(ds[SCHNIDER].raw), "covariate_model": {"derived_inputs": [
        {"quantity": "lbm", "equation": "james_1976", "used_for": ["NOT_A_PARAM"],
         "verbatim": True, "tier": "B", "extraction": {"review_status": "unverified"}}]}}
    problems = _check_covariate_model(Model(raw), ds, set(ds.citations))
    assert any("not a parameter" in p for p in problems)


def test_validator_flags_status_mismatch(ds):
    raw = {**dict(ds[SCHNIDER].raw), "covariate_sensitivity_status": "none"}
    problems = _check_covariate_model(Model(raw), ds, set(ds.citations))
    assert any("expected 'declared'" in p for p in problems)


def test_dataset_validates_clean(ds):
    # the whole dataset (with every v0.7 binding) is valid end to end
    assert hypnos.validate_dataset(ds) == []


# --------------------------------------------------------------------------- #
# Verification checklist gains the covariate_equation group (the 5 traps)
# --------------------------------------------------------------------------- #
def test_verification_has_covariate_equation_group(ds):
    mv = hypnos.model_verification(ds, SCHNIDER)
    groups = {it.group for it in mv.checklist}
    assert "covariate_equation" in groups
    labels = " ".join(it.label for it in mv.checklist if it.group == "covariate_equation")
    for trap in ("Trap 1", "Trap 2", "Trap 3", "Trap 4", "Trap 5"):
        assert trap in labels
    md = hypnos.checklist_markdown(mv)
    assert "Covariate-model bindings" in md
