"""External-validation metric engine — Varvel's framework (v0.4 §6).

The metric math is checked against hand-computed answers on small synthetic
fixtures (including the textbook edge cases the spec names: a single sample,
all-zero error, monotone drift for divergence, a non-positive prediction). The
``validate_against_cohort`` harness is exercised end-to-end against the *real*
reference kernel via a self-consistency fixture: a model's own prediction, fed
back as the observations, must score ~0 error — proving the alignment + solver
wiring with no invented clinical numbers.
"""
import math

import numpy as np
import pytest

import hypnos
from hypnos.analysis import (
    SubjectRecord,
    performance_error,
    pooled_performance,
    validate_against_cohort,
    varvel_metrics,
)


# --------------------------------------------------------------------------- #
# performance_error
# --------------------------------------------------------------------------- #
def test_performance_error_basic():
    # obs/pred = [110/100, 90/100, 105/100] -> PE = [+10, -10, +5] %
    pe = performance_error([110.0, 90.0, 105.0], [100.0, 100.0, 100.0])
    np.testing.assert_allclose(pe, [10.0, -10.0, 5.0])


def test_performance_error_zero_prediction_is_nan():
    # the prediction is the denominator; a non-positive prediction is undefined.
    pe = performance_error([5.0, 5.0, 5.0], [0.0, -1.0, 2.0])
    assert math.isnan(pe[0]) and math.isnan(pe[1])
    assert pe[2] == pytest.approx(150.0)


def test_performance_error_shape_mismatch_raises():
    with pytest.raises(ValueError):
        performance_error([1.0, 2.0], [1.0])


# --------------------------------------------------------------------------- #
# varvel_metrics
# --------------------------------------------------------------------------- #
def test_varvel_metrics_hand_computed():
    pe = [10.0, -10.0, 5.0]
    m = varvel_metrics(pe)
    # MDPE = median(10, -10, 5) = 5; MDAPE = median(10, 10, 5) = 10
    assert m.mdpe == pytest.approx(5.0)
    assert m.mdape == pytest.approx(10.0)
    # wobble = median(|PE - MDPE|) = median(|10-5|, |-10-5|, |5-5|) = median(5,15,0) = 5
    assert m.wobble == pytest.approx(5.0)
    assert m.n == 3
    assert math.isnan(m.divergence)  # no times supplied


def test_varvel_divergence_monotone_drift():
    # |PE| = [0, 10, 20] at t = [0, 60, 120] min -> slope 10%/60min = 10 %/h exactly.
    pe = [0.0, 10.0, 20.0]
    times = [0.0, 60.0, 120.0]
    m = varvel_metrics(pe, times)
    assert m.divergence == pytest.approx(10.0)


def test_varvel_all_zero_error():
    m = varvel_metrics([0.0, 0.0, 0.0], [0.0, 30.0, 60.0])
    assert m.mdpe == 0.0 and m.mdape == 0.0 and m.wobble == 0.0
    assert m.divergence == pytest.approx(0.0)


def test_varvel_single_sample():
    m = varvel_metrics([12.0], [0.0])
    assert m.mdpe == pytest.approx(12.0)
    assert m.mdape == pytest.approx(12.0)
    assert m.wobble == pytest.approx(0.0)
    assert math.isnan(m.divergence)  # < 2 samples -> undefined
    assert m.n == 1


def test_varvel_ignores_nan_samples():
    # a non-positive prediction produced a nan PE; it must drop out of the metrics.
    pe = performance_error([110.0, 5.0, 105.0], [100.0, 0.0, 100.0])
    m = varvel_metrics(pe, [0.0, 30.0, 60.0])
    assert m.n == 2
    assert m.mdpe == pytest.approx(7.5)  # median(10, 5)


def test_varvel_empty_is_nan():
    m = varvel_metrics([np.nan, np.nan])
    assert m.n == 0
    assert math.isnan(m.mdpe) and math.isnan(m.mdape)


# --------------------------------------------------------------------------- #
# pooled_performance
# --------------------------------------------------------------------------- #
def test_pooled_population_median_and_determinism():
    subjects = [
        varvel_metrics([10.0, 10.0]),   # mdpe 10, mdape 10
        varvel_metrics([20.0, 20.0]),   # mdpe 20, mdape 20
        varvel_metrics([30.0, 30.0]),   # mdpe 30, mdape 30
    ]
    a = pooled_performance(subjects, seed=7)
    b = pooled_performance(subjects, seed=7)
    assert a.mdpe == pytest.approx(20.0)   # median across subjects
    assert a.mdape == pytest.approx(20.0)
    assert a.n_subjects == 3
    # seeded bootstrap is byte-reproducible
    assert a.ci95["mdpe"] == b.ci95["mdpe"]
    lo, hi = a.ci95["mdpe"]
    assert lo <= 20.0 <= hi


def test_pooled_single_subject_degenerate_ci():
    pop = pooled_performance([varvel_metrics([15.0, 15.0])], seed=1)
    assert pop.n_subjects == 1
    assert pop.mdpe == pytest.approx(15.0)
    assert pop.ci95["mdpe"] == (pytest.approx(15.0), pytest.approx(15.0))


def test_pooled_empty():
    pop = pooled_performance([], seed=1)
    assert pop.n_subjects == 0
    assert math.isnan(pop.mdpe)


# --------------------------------------------------------------------------- #
# validate_against_cohort — end-to-end against the real kernel
# --------------------------------------------------------------------------- #
def _self_consistency_subjects(ds, model_id):
    """Build a cohort whose 'observations' ARE the model's own predictions.

    With observations equal to predictions the performance error is identically
    zero, so this fixture validates the harness wiring (dose reconstruction,
    forward simulation, time alignment, metric math) against the real solver —
    without curating or inventing any clinical concentration.
    """
    import numpy as np

    from hypnos.simulate import simulate

    cov = {"age": 40, "weight": 70, "height": 175, "sex": "M"}
    schedule = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    obs_t = np.array([1.0, 3.0, 7.0, 15.0, 30.0])
    res = simulate(ds, model_id, patient=cov, schedule=schedule,
                   t=np.union1d(np.linspace(0, 30, 200), obs_t))
    pred = np.interp(obs_t, res.t, res.cp)
    observations = [(float(t), float(v), "cp") for t, v in zip(obs_t, pred)]
    return [SubjectRecord(covariates=cov, schedule=schedule, observations=observations,
                          subject_id="s1")]


def test_validate_against_cohort_self_consistency():
    ds = hypnos.load()
    mid = "hypnotics_iv.propofol.eleveld_2018"
    subjects = _self_consistency_subjects(ds, mid)
    cv = validate_against_cohort(ds, mid, subjects, target="cp", seed=7)
    assert cv.n_subjects == 1
    assert cv.mode == "pk_concentration"
    # observations == predictions -> ~0 error to numerical/interpolation tolerance
    assert abs(cv.population.mdpe) < 1e-6
    assert cv.population.mdape < 1e-6


def test_validate_against_cohort_known_offset():
    """A uniform +20% on observations yields MDPE = MDAPE = 20% exactly."""
    ds = hypnos.load()
    mid = "hypnotics_iv.propofol.eleveld_2018"
    subjects = _self_consistency_subjects(ds, mid)
    inflated = [
        SubjectRecord(
            covariates=s.covariates, schedule=s.schedule,
            observations=[(t, 1.2 * v, k) for (t, v, k) in s.observations],
        )
        for s in subjects
    ]
    cv = validate_against_cohort(ds, mid, inflated, target="cp", seed=7)
    assert cv.population.mdpe == pytest.approx(20.0, abs=1e-3)
    assert cv.population.mdape == pytest.approx(20.0, abs=1e-3)


def test_validate_against_cohort_to_record_schema_shape():
    ds = hypnos.load()
    mid = "hypnotics_iv.propofol.eleveld_2018"
    cv = validate_against_cohort(ds, mid, _self_consistency_subjects(ds, mid),
                                 target="cp", seed=7, dataset="self_consistency")
    rec = cv.to_record()
    # the record validates against the schema's external_validation_entry $def
    assert rec["mode"] == "pk_concentration"
    assert rec["target"] == "cp"
    assert {m["name"] for m in rec["metrics"]} == {"MDPE", "MDAPE", "wobble", "divergence"}
    assert rec["provenance"]["computed_by"] == "hypnos"
    assert rec["reproducible"] is True


def test_validate_bis_target_requires_pd_model():
    ds = hypnos.load()
    with pytest.raises(ValueError):
        validate_against_cohort(ds, "hypnotics_iv.propofol.eleveld_2018",
                                [], target="bis")


def test_validate_bad_target_raises():
    ds = hypnos.load()
    with pytest.raises(ValueError):
        validate_against_cohort(ds, "hypnotics_iv.propofol.eleveld_2018",
                                [], target="auc")


# --------------------------------------------------------------------------- #
# schema + integrity layer for the external_validation block (v0.4 §4)
# --------------------------------------------------------------------------- #
import json  # noqa: E402
from pathlib import Path  # noqa: E402

from hypnos.load import find_dataset_dir  # noqa: E402
from hypnos.models import Model  # noqa: E402
from hypnos.validate import _check_external_validation  # noqa: E402


def _schema():
    root = find_dataset_dir()
    return json.loads((Path(root) / "schema" / "model.schema.json").read_text())


def _entry(**over):
    e = {
        "dataset": "self_consistency",
        "mode": "pk_concentration",
        "target": "cp",
        "cohort": {"n_subjects": 1, "filter": "all", "in_envelope": True},
        "metrics": [{"name": "MDAPE", "value": 18.9, "units": "%",
                     "ci95": {"low": 17.0, "high": 20.8}}],
        "provenance": {"computed_by": "hypnos", "seed": 7},
        "reproducible": True,
    }
    e.update(over)
    return e


def test_to_record_validates_against_schema():
    import jsonschema

    ds = hypnos.load()
    mid = "hypnotics_iv.propofol.eleveld_2018"
    rec = validate_against_cohort(ds, mid, _self_consistency_subjects(ds, mid),
                                  target="cp", seed=7).to_record()
    sub = {"$ref": "#/$defs/external_validation_entry", "$defs": _schema()["$defs"]}
    jsonschema.validate(rec, sub)  # raises on any schema violation


def test_check_external_validation_accepts_consistent():
    raw = {"id": "hypnotics_iv.propofol.eleveld_2018",
           "validation_status": "external_pk",
           "external_validation": [_entry()]}
    assert _check_external_validation(Model(raw=raw)) == []


def test_check_external_validation_flags_mode_target_mismatch():
    raw = {"id": "x.y.z", "validation_status": "external_pk",
           "external_validation": [_entry(target="bis")]}
    probs = _check_external_validation(Model(raw=raw))
    assert any("inconsistent with mode" in p for p in probs)


def test_check_external_validation_flags_status_none_with_entries():
    raw = {"id": "x.y.z", "validation_status": "none",
           "external_validation": [_entry()]}
    probs = _check_external_validation(Model(raw=raw))
    assert any("validation_status 'none'" in p for p in probs)


def test_check_external_validation_flags_ci_inversion():
    bad = _entry(metrics=[{"name": "MDPE", "value": 1.0, "units": "%",
                           "ci95": {"low": 5.0, "high": 1.0}}])
    raw = {"id": "x.y.z", "validation_status": "external_pk",
           "external_validation": [bad]}
    probs = _check_external_validation(Model(raw=raw))
    assert any("ci95 low" in p for p in probs)


def test_full_dataset_still_validates():
    # no curated model carries the new block yet; the dataset stays clean.
    assert hypnos.validate_dataset() == []


# --------------------------------------------------------------------------- #
# the generic CSV adapter + the self-consistency cohort (the CLI-exposed surface)
# --------------------------------------------------------------------------- #
from hypnos.analysis import (  # noqa: E402
    subjects_from_cohort_self_consistency,
    subjects_from_csv,
    subjects_from_vitaldb,
)


def _synthetic_vitaldb_case():
    """A VitalDB-shaped case: a stepped propofol rate track + a measured-BIS track."""
    # rate (mL/h) sampled every 10 s: a 40 mL/h induction infusion that steps to 20,
    # then 0 — three distinct rates, so three reconstructed infusion events.
    rate = [(0.0, 40.0), (10.0, 40.0), (20.0, 20.0), (30.0, 20.0), (40.0, 0.0)]
    bis = [(0.0, 95.0), (10.0, 70.0), (20.0, 45.0), (30.0, 48.0), (40.0, 60.0),
           (50.0, 300.0)]   # the 300 is an artifact -> must be filtered out
    return {"id": "42", "age": 55, "sex": "F", "height": 165, "weight": 62,
            "infusion_ml_h": rate, "bis": bis}


def test_subjects_from_vitaldb_reconstructs_schedule_and_observations():
    subs = subjects_from_vitaldb([_synthetic_vitaldb_case()], propofol_mg_per_ml=20.0)
    assert len(subs) == 1
    s = subs[0]
    assert s.subject_id == "42"
    assert s.covariates == {"age": 55.0, "height": 165.0, "weight": 62.0, "sex": "F"}
    # one infusion event per RATE CHANGE: 40 -> 20 -> 0  (3 events, in mg/h = mL/h*20)
    assert s.schedule == [("infusion", 0.0, "800 mg/h"),
                          ("infusion", round(20.0 / 60.0, 4), "400 mg/h"),
                          ("infusion", round(40.0 / 60.0, 4), "0 mg/h")]
    # observations are measured BIS in MINUTES, with the >100 artifact dropped
    assert all(kind == "bis" and 0 < v <= 100 for (_, v, kind) in s.observations)
    assert len(s.observations) == 5                      # the 300 artifact removed
    assert s.observations[0] == (0.0, 95.0, "bis")


def test_subjects_from_vitaldb_drops_cases_without_bis_or_infusion():
    no_bis = {"id": "1", "age": 50, "weight": 70, "infusion_ml_h": [(0.0, 10.0)], "bis": []}
    no_inf = {"id": "2", "age": 50, "weight": 70, "infusion_ml_h": [], "bis": [(0.0, 50.0)]}
    assert subjects_from_vitaldb([no_bis, no_inf]) == []


def test_subjects_from_vitaldb_validates_end_to_end():
    # the reconstructed cohort runs through the real PK->BIS stack + Varvel engine
    ds = hypnos.load()
    subs = subjects_from_vitaldb([_synthetic_vitaldb_case()])
    cv = validate_against_cohort(ds, "hypnotics_iv.propofol.eleveld_2018", subs,
                                 target="bis", pd_model="pd_effect.propofol.eleveld_bis", seed=0)
    assert cv.mode == "pd_bis" and cv.target == "bis" and cv.n_subjects == 1
    assert np.isfinite(cv.population.mdape)


def test_subjects_from_csv_groups_and_parses():
    rows = [
        {"subject": "s1", "time_min": "5", "observed": "3.1", "kind": "cp",
         "age": "50", "weight": "70", "height": "170", "sex": "M",
         "bolus": "2 mg/kg", "infusion": "6 mg/kg/h"},
        {"subject": "s1", "time_min": "15", "observed": "2.4", "kind": "cp",
         "age": "50", "weight": "70", "height": "170", "sex": "M"},
        {"subject": "s2", "time_min": "5", "observed": "2.8", "kind": "cp",
         "age": "68", "weight": "60", "height": "165", "sex": "F", "bolus": "2 mg/kg"},
    ]
    subjects = subjects_from_csv(rows)
    assert [s.subject_id for s in subjects] == ["s1", "s2"]          # order preserved
    assert len(subjects[0].observations) == 2                        # grouped by subject
    assert subjects[0].covariates == {"age": 50.0, "weight": 70.0, "height": 170.0, "sex": "M"}
    assert ("bolus", 0.0, "2 mg/kg") in subjects[0].schedule
    assert ("infusion", 0.0, "6 mg/kg/h") in subjects[0].schedule
    assert subjects[1].schedule == [("bolus", 0.0, "2 mg/kg")]       # infusion blank -> omitted


def test_subjects_from_csv_skips_malformed_rows():
    rows = [
        {"subject": "s1", "time_min": "5", "observed": "3.1", "kind": "cp"},
        {"subject": "s1", "time_min": "", "observed": "x", "kind": "cp"},   # malformed -> skipped
        {"subject": "s1", "observed": "2.0", "kind": "cp"},                 # no time -> skipped
    ]
    subjects = subjects_from_csv(rows)
    assert len(subjects) == 1 and len(subjects[0].observations) == 1


def test_self_consistency_recovers_the_known_offset():
    # a +X% offset of the model's own prediction MUST recover MDPE ≈ +X% and
    # MDAPE ≈ |X|% — the engine's known-answer check, no external data (v0.4 §8).
    ds = hypnos.load()
    for offset in (20.0, -12.0):
        subs = subjects_from_cohort_self_consistency(
            ds, "hypnotics_iv.propofol.eleveld_2018", target="cp", offset_pct=offset)
        cv = validate_against_cohort(ds, "hypnotics_iv.propofol.eleveld_2018", subs,
                                     target="cp", seed=0)
        assert cv.population.mdpe == pytest.approx(offset, abs=1e-6)
        assert cv.population.mdape == pytest.approx(abs(offset), abs=1e-6)
        assert cv.population.wobble == pytest.approx(0.0, abs=1e-6)   # no intra-subject scatter


# --------------------------------------------------------------------------- #
# Envelope stratification + the cross-model leaderboard (v0.4 VE2/VE3)
# --------------------------------------------------------------------------- #
from hypnos.analysis import (  # noqa: E402
    cross_model_leaderboard,
    partition_by_envelope,
)

ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
MARSH = "hypnotics_iv.propofol.marsh_1991"
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
BIS_PD = "pd_effect.propofol.bis_sigmoid"


def _subj(cov):
    return SubjectRecord(covariates=cov, schedule=[("bolus", 0.0, "2 mg/kg")],
                         observations=[(5.0, 40.0, "bis")], subject_id="x")


def test_partition_by_envelope_is_per_model():
    ds = hypnos.load()
    healthy = _subj(dict(age=40, weight=75, height=178, sex="M"))      # in-envelope everywhere
    obese = _subj(dict(age=40, weight=140, height=172, sex="M"))       # BMI 47 -> outside Schnider [20,42]
    in_s, out_s = partition_by_envelope(ds, SCHNIDER, [healthy, obese])
    assert healthy in in_s and obese in out_s                          # Schnider greys the obese subject
    in_e, out_e = partition_by_envelope(ds, ELEVELD, [healthy, obese])
    assert out_e == []                                                 # Eleveld's envelope covers both


def test_leaderboard_self_consistency_ranks_source_first():
    # cohort = Eleveld's own BIS predictions (0% offset) -> Eleveld must score ~0 and rank #1.
    ds = hypnos.load()
    subs = subjects_from_cohort_self_consistency(ds, ELEVELD, target="bis",
                                                 pd_model=BIS_PD, offset_pct=0.0)
    cands = [(ELEVELD, BIS_PD), (MARSH, BIS_PD), (SCHNIDER, BIS_PD)]
    lb = cross_model_leaderboard(ds, subs, cands, target="bis", seed=7, dataset="sc")
    assert [e.model_id for e in lb.entries][0] == ELEVELD             # ranked best (lowest MDAPE)
    assert lb.entries[0].overall.population.mdape == pytest.approx(0.0, abs=1e-6)
    # ordering is monotone in MDAPE (nan last)
    mdapes = [e.mdape for e in lb.entries]
    assert mdapes == sorted(mdapes)


def test_leaderboard_is_deterministic():
    ds = hypnos.load()
    subs = subjects_from_cohort_self_consistency(ds, ELEVELD, target="bis", pd_model=BIS_PD, offset_pct=10.0)
    cands = [(ELEVELD, BIS_PD), (MARSH, BIS_PD)]
    a = cross_model_leaderboard(ds, subs, cands, target="bis", seed=3).to_record()
    b = cross_model_leaderboard(ds, subs, cands, target="bis", seed=3).to_record()
    assert a == b                                                     # same (subjects, candidates, seed)


def test_leaderboard_skips_kernel_pending_candidate():
    # a kernel-pending model raises NotImplementedError in the engine -> skipped, never a
    # fabricated number (the never-invent rule carried into the leaderboard).
    ds = hypnos.load()
    subs = subjects_from_cohort_self_consistency(ds, ELEVELD, target="bis", pd_model=BIS_PD, offset_pct=0.0)
    cands = [(ELEVELD, BIS_PD), ("nmb_agents.succinylcholine.roy_2004", BIS_PD)]
    lb = cross_model_leaderboard(ds, subs, cands, target="bis", seed=1)
    ids = [e.model_id for e in lb.entries]
    assert ELEVELD in ids and "nmb_agents.succinylcholine.roy_2004" not in ids


def test_leaderboard_record_shape():
    ds = hypnos.load()
    subs = subjects_from_cohort_self_consistency(ds, ELEVELD, target="bis", pd_model=BIS_PD, offset_pct=5.0)
    lb = cross_model_leaderboard(ds, subs, [(ELEVELD, BIS_PD)], target="bis", seed=0, dataset="sc")
    rec = lb.to_record()
    assert rec["target"] == "bis" and rec["leaderboard"][0]["rank"] == 1
    assert rec["leaderboard"][0]["overall"]["metrics"][0]["name"] == "MDPE"


# --------------------------------------------------------------------------- #
# Open-TCI cp adapter (v0.4 VE2) + tier falsification (v0.4 §5)
# --------------------------------------------------------------------------- #
from hypnos.analysis import subjects_from_opentci  # noqa: E402
from hypnos.validate import _check_tier_falsification  # noqa: E402


def test_opentci_adapter_maps_to_cp_subjects():
    cases = [{
        "id": "ext1", "age": 40, "sex": "M", "height": 175, "weight": 75,
        "doses": [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "10 mg/kg/h")],
        "cp": [(5.0, 3.1), (10.0, 2.6), (None, 9.9), (20.0, -1.0)],   # bad rows dropped
    }, {"id": "empty", "age": 50, "weight": 70, "cp": [], "doses": []}]  # dropped (no cp/dose)
    subs = subjects_from_opentci(cases)
    assert len(subs) == 1                                     # the empty case is dropped
    s = subs[0]
    assert s.covariates["weight"] == 75.0 and s.covariates["sex"] == "M"
    assert [o[2] for o in s.observations] == ["cp", "cp"]     # only the 2 valid cp rows
    assert s.schedule[0] == ("bolus", 0.0, "2 mg/kg")


def test_opentci_cp_round_trips_through_the_engine():
    # build a known-answer cp cohort from Eleveld's own prediction (+0%) -> MDPE ~ 0.
    ds = hypnos.load()
    import numpy as np
    grid = np.linspace(0.0, 30.0, 200)
    sim = hypnos.simulate(ds, ELEVELD, patient=dict(age=40, weight=75, height=175, sex="M"),
                          schedule=[("bolus", 0.0, "2 mg/kg")], t=grid)
    obs_t = [5.0, 10.0, 20.0, 30.0]
    cp = [(t, float(np.interp(t, grid, sim.cp))) for t in obs_t]
    cases = [{"id": "k1", "age": 40, "sex": "M", "height": 175, "weight": 75,
              "doses": [("bolus", 0.0, "2 mg/kg")], "cp": cp}]
    subs = subjects_from_opentci(cases)
    cv = validate_against_cohort(ds, ELEVELD, subs, target="cp", dataset="opentci_propofol", seed=0)
    assert cv.target == "cp" and cv.mode == "pk_concentration"
    # observations ARE the model's own cp -> ~0 error (to interpolation-grid tolerance)
    assert cv.population.mdpe == pytest.approx(0.0, abs=0.05)
    assert cv.population.mdape == pytest.approx(0.0, abs=0.05)


def _val_model(tier, mdape, in_env=True):
    return Model({"id": "x.y.z", "drug": {"name": "propofol"}, "purpose": "pk", "tier": tier,
                  "primary_citation": "eleveld-2018-propofol",
                  "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                  "parameters": [], "extraction": {"review_status": "unverified"},
                  "validation_status": "external_pk",
                  "external_validation": [{
                      "dataset": "opentci_propofol", "mode": "pk_concentration", "target": "cp",
                      "cohort": {"in_envelope": in_env},
                      "metrics": [{"name": "MDAPE", "value": mdape, "units": "%"}],
                      "provenance": {"computed_by": "hypnos"}, "reproducible": True}]})


def test_tier_falsification_flags_in_envelope_overshoot():
    # a Tier-A model with in-envelope MDAPE 45% > the ~30% Tier-A band -> advisory flag
    probs = _check_tier_falsification(_val_model("A", 45.0))
    assert any("tier-mismatch" in p and "Tier A" in p for p in probs)


def test_tier_falsification_silent_when_within_band():
    assert _check_tier_falsification(_val_model("A", 22.0)) == []        # comfortably inside


def test_tier_falsification_ignores_out_of_envelope():
    # out-of-envelope poor accuracy is EXPECTED (the failure mode), not a tier mismatch
    assert _check_tier_falsification(_val_model("A", 60.0, in_env=False)) == []


def test_visual_predictive_check_bins_observed_and_predicted():
    from hypnos.analysis import subjects_from_cohort_self_consistency, visual_predictive_check
    ds = hypnos.load()
    subs = subjects_from_cohort_self_consistency(ds, ELEVELD, target="cp", offset_pct=0.0, n_subjects=3)
    v = visual_predictive_check(ds, ELEVELD, subs, target="cp", n_bins=5, seed=7)
    assert len(v.bin_centers) == 5
    assert v.n_per_bin.sum() > 0                              # observations were binned
    # Eleveld carries BSV -> a real predicted band (10/50/90), not a collapsed line
    assert 10 in v.predicted and 90 in v.predicted
    # observations were generated from the model itself (0% offset) -> observed median sits
    # inside the predicted 10-90 band wherever both are defined
    for b in range(5):
        if v.n_per_bin[b] > 0:
            assert v.predicted[10][b] <= v.observed[50][b] <= v.predicted[90][b] + 1e-6
