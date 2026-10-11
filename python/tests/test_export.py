import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET

import json

import numpy as np
import pytest

import hypnos
from hypnos.export import FORMATS, export_model
from hypnos.export.registry import instantiate
from hypnos.reference import Dosing, MicroParams, simulate

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
MARSH = "hypnotics_iv.propofol.marsh_1991"
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"

# Every PK model is exported (CI regenerates the matrix and bundles it in the
# .omex), including the kernel-pending ones (fentanyl, rocuronium) whose exporters
# take a distinct "no instantiated parameters" branch — so the well-formedness
# guarantee must span all of them, not just the kernel-backed adult models.
_PK_MODEL_IDS = sorted(m.id for m in hypnos.load() if m.purpose == "pk")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# Model exports carry the banner; bibtex (citations) and omex (binary archive,
# tested in test_combine_bibtex) are excluded.
_BANNER_FORMATS = [f for f in FORMATS if f not in ("bibtex", "csv", "omex")]


@pytest.mark.parametrize("fmt", _BANNER_FORMATS)
@pytest.mark.parametrize("mid", [SCHNIDER, MARSH])
def test_every_export_carries_clinical_use(ds, fmt, mid):
    _, text = export_model(fmt, ds[mid], ds)
    assert "PROHIBITED" in text  # universal clinicalUse flag


@pytest.mark.parametrize("fmt", ["sbml", "pharmml"])
@pytest.mark.parametrize("mid", _PK_MODEL_IDS)
def test_xml_exports_are_well_formed(ds, fmt, mid):
    _, text = export_model(fmt, ds[mid], ds)
    minidom.parseString(text)  # raises if malformed — covers kernel-pending models too


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def test_sbml_round_trips_against_kernel(ds):
    m = ds[SCHNIDER]
    patient = {"age": 50, "weight": 77, "height": 177, "sex": "M"}
    _, text = export_model("sbml", m, ds, patient)

    # Recover the micro-rate constants from the SBML <parameter> block.
    root = ET.fromstring(text)
    params = {}
    for el in root.iter():
        if _local(el.tag) == "parameter":
            params[el.attrib["id"]] = float(el.attrib["value"])
    recovered = MicroParams(
        V1=params["V1"], k10=params["k10"], k12=params["k12"], k21=params["k21"],
        k13=params["k13"], k31=params["k31"], ke0=params["ke0"],
    )

    direct = instantiate(m, patient)
    dosing = Dosing(boluses=((0.0, 150.0),), infusions=((0.0, 6.0),))
    t = np.linspace(0, 60, 121)
    a = simulate(direct, dosing, t)
    b = simulate(recovered, dosing, t)
    # algebraic round-trip tolerance (spec §6: ~1e-6)
    assert np.allclose(a.cp, b.cp, rtol=1e-6, atol=1e-9)
    assert np.allclose(a.ce, b.ce, rtol=1e-6, atol=1e-9)


def test_tci_json_round_trips_against_kernel(ds):
    m = ds[MARSH]
    patient = {"age": 50, "weight": 80, "height": 180, "sex": "M"}
    _, text = export_model("tci_json", m, ds, patient)
    doc = json.loads(text)
    mc = doc["instantiated_parameters"]["micro_rate_constants"]
    recovered = MicroParams(**{k: mc[k] for k in
                               ("V1", "k10", "k12", "k21", "k13", "k31", "ke0")})
    direct = instantiate(m, patient)
    dosing = Dosing(boluses=((0.0, 160.0),))
    t = np.linspace(0, 30, 100)
    assert np.allclose(simulate(direct, dosing, t).cp,
                       simulate(recovered, dosing, t).cp, rtol=1e-9, atol=1e-12)


def test_nonmem_theta_matches_instantiation(ds):
    m = ds[SCHNIDER]
    patient = {"age": 50, "weight": 77, "height": 177, "sex": "M"}
    _, text = export_model("nonmem", m, ds, patient)
    assert "$THETA" in text and "ADVAN11" in text
    vc = instantiate(m, patient).as_volumes_clearances()
    # CL theta line: "<value>   ; 1 CL  (L/min)"
    cl_line = next(ln for ln in text.splitlines() if "; 1 CL" in ln)
    cl_val = float(cl_line.split(";")[0].strip())
    assert abs(cl_val - vc["Cl1"]) < 1e-3


def test_tci_json_marks_pending_kernel(ds):
    # fentanyl (Shafer) is still kernel-pending; its TCI-JSON omits instantiated params
    _, text = export_model("tci_json", ds["opioids.fentanyl.shafer_1990"], ds)
    doc = json.loads(text)
    assert doc["instantiated_parameters"] is None
    assert "pending" in doc["kernel_status"]


# --------------------------------------------------------------------------- #
# v0.1 §7 — an export carries the tier PROPAGATED to the patient it is instantiated at
# --------------------------------------------------------------------------- #
import re  # noqa: E402
import zipfile  # noqa: E402
from io import BytesIO  # noqa: E402

from hypnos.export import combine  # noqa: E402
from hypnos.simulate import simulate as _simulate_api  # noqa: E402

KATARIA = "hypnotics_iv.propofol.kataria_1994"      # pediatric (3–11 y), declared tier C
_MODEL_FORMATS = ["nonmem", "pharmml", "sbml", "tci_json", "rxode2", "pumas"]
_CHILD = {"age": 6, "weight": 22, "height": 116, "sex": "M"}
_OBESE = {"age": 50, "weight": 140, "height": 165, "sex": "M"}     # BMI ~51


def _exported_tier(fmt, text):
    if fmt == "tci_json":
        return json.loads(text)["provenance"]["hypnos:confidenceTier"]
    m = re.search(r"onfidenceTier\W+(\w)\b", text) or re.search(r"\(tier (\w),", text)
    return m.group(1)


@pytest.mark.parametrize("fmt", _MODEL_FORMATS)
def test_pediatric_model_at_the_adult_reference_patient_is_exported_tier_d(ds, fmt):
    # the default export instantiates Kataria for a 50-year-old: that artifact is an
    # extrapolation and must say so, not carry the model's in-envelope tier
    text = export_model(fmt, ds[KATARIA], ds)[1]
    assert _exported_tier(fmt, text) == "D"
    assert "nvelopeWarning" in text and "outside derivation envelope" in text
    assert re.search(r"eclaredTier\W+C\b|declared C\)", text)       # the declared tier stays visible


@pytest.mark.parametrize("fmt", _MODEL_FORMATS)
def test_same_model_in_envelope_keeps_its_declared_tier_and_no_warning(ds, fmt):
    text = export_model(fmt, ds[KATARIA], ds, _CHILD)[1]
    assert _exported_tier(fmt, text) == ds[KATARIA].tier == "C"
    assert "nvelopeWarning" not in text and "eclaredTier" not in text
    assert "OUTSIDE THIS MODEL'S ENVELOPE" not in text


@pytest.mark.parametrize("fmt", _MODEL_FORMATS)
@pytest.mark.parametrize("mid,patient", [(SCHNIDER, _OBESE), (KATARIA, None), (SCHNIDER, None),
                                         (KATARIA, _CHILD)])
def test_exported_tier_is_never_more_confident_than_the_simulation(ds, fmt, mid, patient):
    # one shared safety evaluation: the label on the artifact equals simulate()'s tier
    from hypnos.export._common import REFERENCE_PATIENT
    pat = patient or dict(REFERENCE_PATIENT)
    sim = _simulate_api(ds, mid, patient=pat, schedule=[("bolus", 0.0, "2 mg/kg")],
                        t=np.linspace(0, 5, 11))
    assert _exported_tier(fmt, export_model(fmt, ds[mid], ds, patient)[1]) == sim.tier


def test_schnider_in_obesity_export_names_the_failure_mode(ds):
    text = export_model("sbml", ds[SCHNIDER], ds, _OBESE)[1]
    ET.fromstring(text)
    assert "<hypnos:confidenceTier>D</hypnos:confidenceTier>" in text
    assert "<hypnos:declaredTier>B</hypnos:declaredTier>" in text
    assert "FAILURE MODE" in text and "bmi" in text.lower()


def test_envelope_warnings_never_change_the_exported_numbers(ds):
    # the label changes; the instantiated parameters are still exactly the kernel's
    doc = json.loads(export_model("tci_json", ds[KATARIA], ds)[1])
    from hypnos.export._common import REFERENCE_PATIENT
    direct = instantiate(ds[KATARIA], dict(REFERENCE_PATIENT))
    assert doc["instantiated_parameters"]["micro_rate_constants"]["k10"] == pytest.approx(direct.k10)
    assert doc["provenance"]["hypnos:declaredTier"] == "C"
    assert len(doc["provenance"]["hypnos:envelopeWarning"]) >= 2


def test_omex_metadata_carries_the_propagated_tier(ds):
    z = zipfile.ZipFile(BytesIO(combine.build_model_archive(ds[KATARIA], ds)))
    meta = z.read("metadata.rdf").decode()
    assert "<hypnos:confidenceTier>D</hypnos:confidenceTier>" in meta
    assert "<hypnos:declaredTier>C</hypnos:declaredTier>" in meta
    z2 = zipfile.ZipFile(BytesIO(combine.build_model_archive(ds[KATARIA], ds, _CHILD)))
    assert "<hypnos:confidenceTier>C</hypnos:confidenceTier>" in z2.read("metadata.rdf").decode()


def test_cli_export_accepts_a_patient(ds, tmp_path):
    from hypnos.cli import main
    base = ["export", "--format", "tci_json", "--model", KATARIA]
    assert main(base + ["--output", str(tmp_path / "adult")]) == 0
    assert main(base + ["--output", str(tmp_path / "child"), "--age", "6", "--weight", "22",
                        "--height", "116"]) == 0
    name = "hypnotics_iv_propofol_kataria_1994.tci.json"
    adult = json.loads((tmp_path / "adult" / name).read_text())
    child = json.loads((tmp_path / "child" / name).read_text())
    assert adult["provenance"]["hypnos:confidenceTier"] == "D"
    assert child["provenance"]["hypnos:confidenceTier"] == "C"
    assert child["instantiated_for"]["age"] == 6.0
    # no flags -> the library default patient, byte-identical to export_model()
    assert (tmp_path / "adult" / name).read_text() == export_model("tci_json", ds[KATARIA], ds)[1]


# --------------------------------------------------------------------------- #
# SBML structural consistency (the rules libSBML would reject, checked without it)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mid", _PK_MODEL_IDS)
@pytest.mark.parametrize("patient", [None, _CHILD, _OBESE], ids=["reference", "child", "obese"])
def test_sbml_rules_reference_only_declared_non_constant_parameters(ds, mid, patient):
    dom = minidom.parseString(export_model("sbml", ds[mid], ds, patient)[1])
    params = {p.getAttribute("id"): p for p in dom.getElementsByTagName("parameter")}
    ids = [p.getAttribute("id") for p in dom.getElementsByTagName("parameter")]
    assert len(ids) == len(set(ids)), "duplicate parameter ids"
    for p in params.values():
        float(p.getAttribute("value"))                      # every parameter has a numeric value
    ruled = []
    for tag in ("rateRule", "assignmentRule"):
        for rule in dom.getElementsByTagName(tag):
            var = rule.getAttribute("variable")
            ruled.append(var)
            assert var in params, f"{tag} targets undeclared {var}"
            assert params[var].getAttribute("constant") == "false", f"{var} is ruled but constant"
            for ci in rule.getElementsByTagName("ci"):
                assert ci.firstChild.data.strip() in params, f"{tag} {var} reads an undeclared symbol"
    assert len(ruled) == len(set(ruled)), "a variable is the target of more than one rule"
    # every units attribute names a unitDefinition in this document or an SBML base unit
    units = {u.getAttribute("id") for u in dom.getElementsByTagName("unitDefinition")}
    builtin = {"dimensionless", "litre", "kilogram", "second", "mole", "gram", "metre"}
    for p in params.values():
        if p.hasAttribute("units"):
            assert p.getAttribute("units") in units | builtin, f"undefined unit on {p.getAttribute('id')}"


@pytest.mark.parametrize("fmt", _MODEL_FORMATS)
def test_export_refuses_a_physically_impossible_patient(ds, fmt):
    with pytest.raises(ValueError, match="weight must be a finite number > 0"):
        export_model(fmt, ds[SCHNIDER], ds, {"age": 50, "weight": 0, "height": 177, "sex": "M"})
