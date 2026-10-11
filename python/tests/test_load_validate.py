import pytest

import hypnos
from hypnos.validate import assert_valid, validate_dataset


def test_load_models():
    ds = hypnos.load()
    assert len(ds) >= 4
    assert "hypnotics_iv.propofol.schnider_1998" in ds
    m = ds["hypnotics_iv.propofol.schnider_1998"]
    assert m.drug_name == "propofol"
    assert m.n_compartments == 3
    assert m.has_effect_compartment


def test_drug_and_citation_lookup():
    ds = hypnos.load()
    assert ds.drug("propofol")["atc"] == "N01AX10"
    cit = ds.citation("schnider-1998-propofol-pk")
    assert cit["doi"].startswith("10.1097")


def test_dataset_is_valid():
    assert validate_dataset() == []
    assert_valid()  # raises on problems


def test_tier_invariant_detects_violation():
    # Construct a model whose record tier is better than its worst parameter tier.
    import copy
    from hypnos.models import worst_tier

    ds = hypnos.load()
    raw = copy.deepcopy(ds["hypnotics_iv.propofol.marsh_1991"].raw)
    raw["parameters"][0]["tier"] = "D"  # worst now D, record still B
    assert worst_tier([p["tier"] for p in raw["parameters"]]) == "D"


def test_validate_catches_dangling_predictive_performance_citation():
    # A performance number must be traceable to a real citation (spec §5/§9):
    # validate() flags a predictive_performance entry citing an unknown record.
    # Deep-copy so we never mutate the process-wide cached dataset.
    import copy

    ds = copy.deepcopy(hypnos.load())
    m = ds["opioids.remifentanil.minto_1997"]
    m.raw.setdefault("predictive_performance", []).append(
        {"metric": "MDAPE", "value": 99.0, "units": "%", "population": "x", "citation": "no-such-citation"}
    )
    problems = validate_dataset(ds)
    assert any("no-such-citation" in p and "predictive_performance" in p for p in problems)


def test_summary_counts():
    s = hypnos.summary(hypnos.load())
    assert s["n_models"] >= 4
    assert s["kernels_implemented"] >= 2
    assert set(s["by_tier"]).issubset({"A", "B", "C", "D"})


# --------------------------------------------------------------------------- #
# Library-level input hardening: unknown ids and physically impossible patients
# --------------------------------------------------------------------------- #
def test_unknown_model_error_is_a_keyerror_that_suggests_the_nearest_ids():
    import hypnos
    from hypnos.load import UnknownModelError
    ds = hypnos.load()
    with pytest.raises(KeyError) as exc:                 # still a KeyError for dict-style callers
        ds["hypnotics_iv.propofol.marsh"]
    assert isinstance(exc.value, UnknownModelError)
    assert "did you mean hypnotics_iv.propofol.marsh_1991" in str(exc.value)
    with pytest.raises(UnknownModelError, match="hypnos models"):
        ds["zzz"]
    assert "hypnotics_iv.propofol.marsh_1991" in ds and "zzz" not in ds


@pytest.mark.parametrize("patient,needle", [
    (dict(age=40, weight=0, height=175, sex="M"), "weight must be a finite number > 0"),
    (dict(age=40, weight=-70, height=175, sex="M"), "weight must be"),
    (dict(age=40, weight=float("nan"), height=175, sex="M"), "weight must be"),
    (dict(age=40, weight=float("inf"), height=175, sex="M"), "weight must be"),
    (dict(age=40, weight="70", height=175, sex="M"), "weight must be"),
    (dict(age=40, weight=70, height=0, sex="M"), "height must be"),
    (dict(age=-1, weight=70, height=175, sex="M"), "age must be a finite number >= 0"),
    (dict(age=40, weight=70, height=175, sex="X"), "sex must be 'M' or 'F'"),
    (dict(age=40, weight={"mean": 70, "sd": -2}, height=175, sex="M"), "SD must be >= 0"),
])
def test_simulate_refuses_a_physically_impossible_patient(patient, needle):
    import hypnos
    ds = hypnos.load()
    with pytest.raises(ValueError, match=needle):
        hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018", patient=patient,
                        schedule=[("bolus", 0.0, "2 mg/kg")], t=[0.0, 1.0, 2.0])


def test_unusual_but_possible_patients_are_the_envelopes_job_not_an_error():
    # a neonate (age 0) and a 250 kg adult are simulated — tiered down with a reason, not refused
    import hypnos
    ds = hypnos.load()
    kw = dict(schedule=[("bolus", 0.0, "2 mg/kg")], t=[0.0, 1.0, 2.0])
    neonate = hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018",
                              patient=dict(age=0, weight=3.2, height=50, sex="F"), **kw)
    heavy = hypnos.simulate(ds, "hypnotics_iv.propofol.marsh_1991",
                            patient=dict(age=40, weight=250, height=180, sex="male"), **kw)
    assert neonate.tier == "D" and heavy.tier == "D"
    assert any("ENVELOPE" in w for w in heavy.warnings)


def test_band_arguments_are_validated():
    import hypnos
    ds = hypnos.load()
    kw = dict(patient=dict(age=40, weight=70, height=175, sex="M"),
              schedule=[("bolus", 0.0, "2 mg/kg")], t=[0.0, 1.0, 2.0], bands=True, seed=1)
    for bad in (dict(samples=0), dict(percentile=(95, 5)), dict(percentile=(5, 101))):
        with pytest.raises(ValueError, match="a band needs"):
            hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018", **kw, **bad)
