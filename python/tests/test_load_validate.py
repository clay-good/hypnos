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


def test_summary_counts():
    s = hypnos.summary(hypnos.load())
    assert s["n_models"] >= 4
    assert s["kernels_implemented"] >= 2
    assert set(s["by_tier"]).issubset({"A", "B", "C", "D"})
