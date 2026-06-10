"""Protein-binding / free-fraction failure mode — v0.5 §B3 / Phase S1.

For a binding-sensitive drug (propofol ~98%, fentanyl ~84%, dexmedetomidine
~94% bound), hypoalbuminemia raises the free (active) fraction, so a model fit in
normal-albumin patients UNDER-estimates effect from a given total concentration.
Hypnos surfaces this as a *cited* failure mode — named, never silently modeled
(the never-invent rule: the free-fraction shift is flagged, not fabricated into a
clearance number).
"""
import pytest

import hypnos
from hypnos.simulate import evaluate_safety, simulate
import numpy as np

PROPOFOL = "hypnotics_iv.propofol.eleveld_2018"
REMI = "opioids.remifentanil.minto_1997"
HYPO = {"age": 55, "weight": 75, "height": 175, "sex": "M", "albumin_g_dl": 2.4}


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --------------------------------------------------------------------------- #
# the curated drug-level chemistry
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("drug, frac", [("propofol", 0.98), ("fentanyl", 0.84),
                                        ("dexmedetomidine", 0.94)])
def test_binding_sensitive_drugs_curated(ds, drug, frac):
    pb = ds.drug(drug)["protein_binding"]
    assert pb["binding_sensitive"] is True
    assert pb["fraction_bound"] == frac
    assert ds.citation(pb["citation"]) is not None   # citation resolves


def test_remifentanil_not_binding_sensitive(ds):
    # ~70% bound: not flagged — the never-invent rule, no claim without standing.
    assert "protein_binding" not in ds.drug("remifentanil")


# --------------------------------------------------------------------------- #
# the caveat surfaces only for a binding-sensitive drug + hypoalbuminemia
# --------------------------------------------------------------------------- #
def _binding_warnings(warnings):
    return [w for w in warnings if w.startswith("BINDING-SENSITIVE")]


def test_caveat_fires_for_binding_sensitive_drug(ds):
    m = ds[PROPOFOL]
    _, warnings, _ = evaluate_safety(m, HYPO, ds.drug(m.drug_name))
    bw = _binding_warnings(warnings)
    assert len(bw) == 1
    assert "98%" in bw[0] and "zamacona-1997-propofol-binding" in bw[0]


def test_no_caveat_for_non_binding_sensitive_drug(ds):
    m = ds[REMI]
    _, warnings, _ = evaluate_safety(m, HYPO, ds.drug(m.drug_name))
    assert _binding_warnings(warnings) == []
    # but the generic albumin extrapolation still fires
    assert any("ALBUMIN EXTRAPOLATION" in w for w in warnings)


def test_no_caveat_without_drug_meta(ds):
    # backward compatible: callers that don't pass drug_meta get no binding caveat.
    m = ds[PROPOFOL]
    _, warnings, _ = evaluate_safety(m, HYPO)
    assert _binding_warnings(warnings) == []


def test_no_caveat_at_normal_albumin(ds):
    m = ds[PROPOFOL]
    normal = {"age": 55, "weight": 75, "height": 175, "sex": "M", "albumin_g_dl": 4.2}
    _, warnings, _ = evaluate_safety(m, normal, ds.drug(m.drug_name))
    assert _binding_warnings(warnings) == []


def test_no_caveat_when_albumin_absent(ds):
    m = ds[PROPOFOL]
    _, warnings, _ = evaluate_safety(m, {"age": 55, "weight": 75}, ds.drug(m.drug_name))
    assert _binding_warnings(warnings) == []


# --------------------------------------------------------------------------- #
# end-to-end through simulate() (which threads drug_meta automatically)
# --------------------------------------------------------------------------- #
def test_simulate_surfaces_binding_caveat(ds):
    res = simulate(ds, PROPOFOL, patient=HYPO,
                   schedule=[("infusion", 0.0, "6 mg/kg/h")], t=np.linspace(0, 20, 81))
    assert any(w.startswith("BINDING-SENSITIVE") for w in res.warnings)


# --------------------------------------------------------------------------- #
# validate: a drug-level protein_binding citation must resolve
# --------------------------------------------------------------------------- #
def test_dataset_validates_with_binding_citations(ds):
    assert hypnos.validate_dataset(ds) == []


def test_per_model_bibtex_includes_binding_and_organ_citations(ds):
    # a per-model bibtex export must carry the drug's protein-binding source and the
    # model's organ-tolerance sources, not just the structural-parameter citations.
    from hypnos.export import bibtex
    prop = bibtex.build(ds, [ds["hypnotics_iv.propofol.eleveld_2018"]])
    assert "zamacona-1997-propofol-binding" in prop
    remi = bibtex.build(ds, [ds["opioids.remifentanil.minto_1997"]])
    assert "dershwitz-1996-remifentanil-hepatic" in remi
    assert "hoke-1997-remifentanil-renal" in remi


def test_validate_flags_unknown_protein_binding_citation():
    from hypnos.validate import validate_dataset

    class _FakeDS:
        citations = {}
        drugs = {"x": {"name": "x", "protein_binding": {"binding_sensitive": True,
                                                        "citation": "nope-not-real"}}}
        def __iter__(self):
            return iter([])   # no models
    probs = validate_dataset(_FakeDS())
    assert any("protein_binding" in p and "nope-not-real" in p for p in probs)
