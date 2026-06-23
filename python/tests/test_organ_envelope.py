"""The organ-function (physiological) envelope — v0.5 §B / Phase S0.

Hypnos's silence on organ-failure patients is made to *speak*: a patient who
declares hepatic / renal / cardiac / albumin impairment greys every model with no
cited standing (Tier D + a named extrapolation), while a model that does have
standing (remifentanil's esterase clearance) carries an explaining note instead.
A normal simulation — no organ covariates — is unaffected.
"""
import numpy as np
import pytest

import hypnos
from hypnos.models import Envelope, Model, Range
from hypnos.simulate import compare, evaluate_safety

PROPOFOL = "hypnotics_iv.propofol.eleveld_2018"
REMI = "opioids.remifentanil.minto_1997"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def _organ_warnings(warnings):
    return [w for w in warnings if w.startswith("ORGAN")]


# --------------------------------------------------------------------------- #
# backward compatibility: no organ covariates => no organ behavior
# --------------------------------------------------------------------------- #
def test_normal_patient_unaffected(ds):
    patient = {"age": 40, "weight": 75, "height": 175, "sex": "M"}
    tier, warnings, excluded = evaluate_safety(ds[PROPOFOL], patient)
    assert tier == "A" and not excluded
    assert _organ_warnings(warnings) == []


def test_envelope_organ_check_empty_for_normal_patient(ds):
    assert ds[PROPOFOL].applicability_envelope.organ_check(
        {"age": 40, "weight": 75}) == []


# --------------------------------------------------------------------------- #
# hepatic: a model with no standing greys; remifentanil retains standing
# --------------------------------------------------------------------------- #
def test_cirrhotic_greys_propofol(ds):
    patient = {"age": 55, "weight": 75, "height": 175, "sex": "M", "child_pugh": "C"}
    tier, warnings, excluded = evaluate_safety(ds[PROPOFOL], patient)
    assert tier == "D" and excluded
    assert any("HEPATIC EXTRAPOLATION" in w for w in warnings)


def test_cirrhotic_remifentanil_retains_standing(ds):
    patient = {"age": 55, "weight": 75, "height": 175, "sex": "M", "child_pugh": "C"}
    tier, warnings, excluded = evaluate_safety(ds[REMI], patient)
    assert tier == "B" and not excluded     # standing preserved, not greyed
    assert any("ORGAN NOTE" in w and "esterases" in w for w in warnings)


def test_child_pugh_a_also_triggers(ds):
    # any Child-Pugh class signals chronic liver disease the model was not fit in.
    _, warnings, excluded = evaluate_safety(
        ds[PROPOFOL], {"age": 55, "weight": 75, "child_pugh": "A"})
    assert excluded and any("HEPATIC EXTRAPOLATION" in w for w in warnings)


# --------------------------------------------------------------------------- #
# renal: remifentanil tolerant (with metabolite caveat), propofol greyed
# --------------------------------------------------------------------------- #
def test_renal_failure_greys_propofol(ds):
    _, warnings, excluded = evaluate_safety(
        ds[PROPOFOL], {"age": 60, "weight": 75, "crcl_ml_min": 15})
    assert excluded and any("RENAL EXTRAPOLATION" in w for w in warnings)


def test_renal_failure_remifentanil_caveat(ds):
    tier, warnings, excluded = evaluate_safety(
        ds[REMI], {"age": 60, "weight": 75, "crcl_ml_min": 15})
    assert tier == "B" and not excluded
    assert any("CAVEAT" in w and "GR90291" in w for w in warnings)


# --------------------------------------------------------------------------- #
# cardiac + albumin: no model has standing (incl. remifentanil) -> all greyed
# --------------------------------------------------------------------------- #
def test_low_ejection_fraction_greys_all(ds):
    patient = {"age": 70, "weight": 75, "ejection_fraction_pct": 30}
    for mid in (PROPOFOL, REMI):
        tier, warnings, excluded = evaluate_safety(ds[mid], patient)
        assert excluded and any("CARDIAC EXTRAPOLATION" in w for w in warnings)


def test_hypoalbuminemia_greys_all(ds):
    patient = {"age": 70, "weight": 75, "albumin_g_dl": 2.4}
    for mid in (PROPOFOL, REMI):
        tier, warnings, excluded = evaluate_safety(ds[mid], patient)
        assert excluded and any("ALBUMIN EXTRAPOLATION" in w for w in warnings)


# --------------------------------------------------------------------------- #
# staging thresholds are exact cut-points
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("crcl, impaired", [(60, False), (59.9, True), (90, False)])
def test_crcl_threshold(ds, crcl, impaired):
    findings = ds[PROPOFOL].applicability_envelope.organ_check(
        {"crcl_ml_min": crcl})
    assert bool(findings) is impaired


@pytest.mark.parametrize("ef, impaired", [(40, False), (39, True)])
def test_ef_threshold(ds, ef, impaired):
    findings = ds[PROPOFOL].applicability_envelope.organ_check(
        {"ejection_fraction_pct": ef})
    assert bool(findings) is impaired


@pytest.mark.parametrize("alb, impaired", [(3.5, False), (3.4, True)])
def test_albumin_threshold(ds, alb, impaired):
    findings = ds[PROPOFOL].applicability_envelope.organ_check(
        {"albumin_g_dl": alb})
    assert bool(findings) is impaired


# --------------------------------------------------------------------------- #
# standing via a fitted numeric range (the forward-compatible mechanism)
# --------------------------------------------------------------------------- #
def test_standing_via_fitted_range():
    env = Envelope(crcl_ml_min=Range(min=10, max=120))
    findings = env.organ_check({"crcl_ml_min": 15})
    assert len(findings) == 1
    assert findings[0].extrapolation is False
    assert "within the model's fitted range" in findings[0].message


def test_below_fitted_range_is_extrapolation():
    env = Envelope(crcl_ml_min=Range(min=30, max=120))
    findings = env.organ_check({"crcl_ml_min": 15})  # below the fitted floor
    assert findings[0].extrapolation is True


# --------------------------------------------------------------------------- #
# compare(): the eligible set shrinks honestly under organ failure (§B6)
# --------------------------------------------------------------------------- #
def test_compare_remifentanil_survives_cirrhosis(ds):
    patient = {"age": 55, "weight": 75, "height": 175, "sex": "M", "child_pugh": "C"}
    t = np.linspace(0, 30, 121)
    schedule = [("infusion", 0.0, "0.2 ug/kg/min")]
    cmp = compare(ds, drug="remifentanil", patient=patient, schedule=schedule, t=t)
    assert cmp.included, "remifentanil models retain standing in hepatic failure"
    assert all("remifentanil" in r.model_id for r in cmp.included)


def test_compare_propofol_all_greyed_in_cirrhosis(ds):
    patient = {"age": 55, "weight": 75, "height": 175, "sex": "M", "child_pugh": "C"}
    t = np.linspace(0, 30, 121)
    schedule = [("infusion", 0.0, "6 mg/kg/h")]
    cmp = compare(ds, drug="propofol", patient=patient, schedule=schedule, t=t)
    assert cmp.included == []                 # nothing has hepatic standing
    assert cmp.excluded                        # everything named + greyed
    assert cmp.divergence["cp"] == {}          # degrades gracefully, no crash
    assert any("HEPATIC EXTRAPOLATION" in r for e in cmp.excluded for r in e["reasons"])


# --------------------------------------------------------------------------- #
# validate: an organ_tolerance citation must resolve
# --------------------------------------------------------------------------- #
def test_dataset_validates_with_organ_tolerance(ds):
    # the curated remifentanil tolerance (Dershwitz/Hoke) resolves cleanly.
    assert hypnos.validate_dataset(ds) == []


class _FakeDataset:
    """Minimal Dataset stand-in: iterable of models + a citations map."""
    def __init__(self, models):
        self._models = models
        self.citations = {"minto-1997-remifentanil": {}}
    def __iter__(self):
        return iter(self._models)


def test_validate_flags_unknown_organ_tolerance_citation():
    from hypnos.validate import validate_dataset
    bad = Model(raw={
        "id": "opioids.x.y", "subsystem": "opioids",
        "primary_citation": "minto-1997-remifentanil",
        "tier": "B", "purpose": "pk",
        "structure": {"compartments": 3, "parameterization": "volumes_clearances"},
        "parameters": [{"symbol": "V1", "value": {"central": 5.0, "units": "L"},
                        "tier": "B", "extraction": {"review_status": "unverified"}}],
        "extraction": {"review_status": "unverified"},
        "applicability_envelope": {
            "organ_tolerance": [{"axis": "hepatic", "basis": "x", "citation": "nope-not-real"}]},
    })
    probs = validate_dataset(_FakeDataset([bad]))
    assert any("organ_tolerance" in p and "nope-not-real" in p for p in probs)


# --------------------------------------------------------------------------- #
# v0.5 §B5 — quantitative disease_state_modifier (opt-in, Tier-D, cited)
# --------------------------------------------------------------------------- #

DEX = "alpha2_agonists.dexmedetomidine.hannivoort_2015"


def test_disease_modifier_curated(ds):
    mods = ds[DEX].disease_state_modifiers
    assert mods and mods[0].dimension == "child_pugh"
    assert mods[0].evidence_tier == "D" and mods[0].applied_by_default is False
    assert mods[0].scale_factor == 0.66 and mods[0].primary_citation


def test_disease_modifier_triggers_only_on_declaration(ds):
    mod = ds[DEX].disease_state_modifiers[0]
    assert mod.triggered_by({"child_pugh": "C"}) is True
    assert mod.triggered_by({"child_pugh": "A"}) is True      # "present" => any class
    assert mod.triggered_by({}) is False                      # undeclared -> never fires


def test_disease_modifier_opt_in_reduces_clearance(ds):
    t = np.linspace(0, 60, 200)
    sched = [("infusion", 0.0, "0.5 mcg/kg/h")]
    normal = hypnos.simulate(ds, DEX, patient=dict(age=50, weight=70), schedule=sched, t=t)
    hep = hypnos.simulate(ds, DEX, patient=dict(age=50, weight=70, child_pugh="C"),
                          schedule=sched, t=t, disease_state=True)
    assert hep.cp_peak > normal.cp_peak                       # reduced clearance -> higher conc
    assert hep.tier == "D"                                    # applying a modifier forces Tier D
    assert any("DISEASE MODIFIER" in w for w in hep.warnings)


def test_disease_modifier_off_by_default(ds):
    # without the opt-in flag, no modifier is applied even if the state is declared
    t = np.linspace(0, 60, 100)
    sched = [("infusion", 0.0, "0.5 mcg/kg/h")]
    res = hypnos.simulate(ds, DEX, patient=dict(age=50, weight=70, child_pugh="C"), schedule=sched, t=t)
    assert not any("DISEASE MODIFIER" in w for w in res.warnings)


def _dz_model(citation, applied=False):
    from hypnos.models import Model
    return Model({"id": "x.y.z", "drug": {"name": "propofol"}, "purpose": "pk", "tier": "B",
                  "primary_citation": "eleveld-2018-propofol",
                  "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                  "parameters": [], "extraction": {"review_status": "unverified"},
                  "disease_state_modifiers": [{
                      "condition": {"dimension": "child_pugh", "value": "C"},
                      "affected_parameter": "Cl", "adjustment": {"type": "scale_factor", "value": 0.6},
                      "evidence_tier": "D", "applied_by_default": applied,
                      "primary_citation": citation}]})


def test_validate_rejects_uncited_disease_modifier():
    from hypnos.validate import _check_disease_modifiers
    # the never-invent rule: a disease adjustment with no citation is rejected (§B5)
    assert any("needs a" in p for p in _check_disease_modifiers(_dz_model(""), set()))


def test_validate_rejects_default_applied_disease_modifier():
    from hypnos.validate import _check_disease_modifiers
    probs = _check_disease_modifiers(_dz_model("eleveld-2018-propofol", applied=True),
                                     {"eleveld-2018-propofol"})
    assert any("applied_by_default" in p for p in probs)
