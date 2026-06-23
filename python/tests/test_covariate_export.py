"""v0.7 C3 — covariate-aware exports.

The derived-covariate equation is part of the model definition, so every export must
carry it: NONMEM computes LBM/FFM the model's way (named, enveloped), the emitted
computation round-trips against the library, SBML carries it as RDF, and TCI-JSON passes
the structured block through verbatim. A model with no covariate_model emits no block.
"""
import json
import xml.dom.minidom as minidom

import pytest

import hypnos
from hypnos.export import nonmem, sbml, tci_json
from hypnos.export._covariate import evaluate_nonmem_expr, to_nonmem_expr

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"   # LBM via james_1976 on Cl1
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"     # FFM via al_sallami_2015 on V3
MARSH = "hypnotics_iv.propofol.marsh_1991"         # no covariate_model


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_nonmem_computes_the_named_covariate(ds):
    ctl = nonmem.build(ds[SCHNIDER], ds)
    assert "LBM via james_1976" in ctl                  # named by equation id
    assert "IF (SEX.EQ.1) THEN" in ctl                  # sex-conditional, not a baked number
    assert "validity envelope: bmi_kg_m2" in ctl        # the envelope rides along
    assert "WT" in ctl.split("$INPUT")[1].split("\n")[0]  # input columns added


def test_no_covariate_model_emits_no_block(ds):
    ctl = nonmem.build(ds[MARSH], ds)
    assert "derived covariate" not in ctl               # an explicit gap, never a fabricated equation


@pytest.mark.parametrize("mid", [SCHNIDER, ELEVELD,
                                 "opioids.remifentanil.minto_1997",
                                 "opioids.remifentanil.kim_2017",
                                 "opioids.remifentanil.eleveld_2017"])
def test_emitted_expression_round_trips_against_library(ds, mid):
    # the exported NONMEM computation must equal hypnos.covariates.evaluate (v0.7 §8).
    m = ds[mid]
    patient = {"weight": 80.0, "height": 175.0, "age": 45.0, "sex": "M"}
    for di in m.covariate_model.derived_inputs:
        rec = ds.covariate_equations[di.equation]
        male = to_nonmem_expr(rec["sex_specific"]["male"])
        emitted = evaluate_nonmem_expr(male, patient)
        lib = hypnos.evaluate_covariate(di.equation, patient).value
        assert emitted == pytest.approx(lib, rel=1e-9), f"{mid}/{di.equation}"


def test_sbml_carries_covariate_rdf_and_stays_well_formed(ds):
    xml = sbml.build(ds[SCHNIDER], ds)
    assert 'hypnos:covariateEquation hypnos:quantity="lbm"' in xml
    assert "hypnos:covariateValidityEnvelope" in xml
    assert "hypnos:covariateSensitivityStatus" in xml
    minidom.parseString(xml)                            # must remain valid XML


def test_tci_json_passes_covariate_model_through(ds):
    doc = json.loads(tci_json.build(ds[SCHNIDER], ds)[1]) if isinstance(
        tci_json.build(ds[SCHNIDER], ds), tuple) else tci_json.build_dict(ds[SCHNIDER], ds)
    assert doc["covariate_model"] is not None           # lossless structured block
    assert doc["covariate_sensitivity_status"] == "declared"
    di = doc["covariate_model"]["derived_inputs"][0]
    assert di["equation"] == "james_1976" and di["quantity"] == "lbm"
