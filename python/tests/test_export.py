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
