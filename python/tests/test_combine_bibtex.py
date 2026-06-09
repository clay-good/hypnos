import io
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
import zipfile

import pytest

import hypnos
from hypnos.export import bibtex, combine

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def _zip(data):
    return zipfile.ZipFile(io.BytesIO(data))


def test_model_omex_contents(ds):
    z = _zip(combine.build_model_archive(ds[SCHNIDER], ds))
    names = z.namelist()
    assert "manifest.xml" in names
    assert any(n.endswith(".sbml.xml") for n in names)
    assert any(n.endswith(".pharmml.xml") for n in names)
    assert any(n.endswith(".tci.json") for n in names)
    assert "metadata.rdf" in names and "citations.bib" in names


def test_omex_manifest_is_valid_and_marks_master(ds):
    z = _zip(combine.build_model_archive(ds[SCHNIDER], ds))
    manifest = z.read("manifest.xml").decode()
    minidom.parseString(manifest)  # well-formed
    root = ET.fromstring(manifest)
    contents = [el.attrib for el in root if el.tag.endswith("content")]
    # the manifest must describe itself and an omex archive
    locs = {c["location"] for c in contents}
    assert "./manifest.xml" in locs and "." in locs
    # exactly one master, and it is the SBML
    masters = [c for c in contents if c.get("master") == "true"]
    assert len(masters) == 1 and masters[0]["location"].endswith(".sbml.xml")


def test_omex_carries_clinical_use_and_sbml_is_wellformed(ds):
    z = _zip(combine.build_model_archive(ds[SCHNIDER], ds))
    assert "PROHIBITED" in z.read("metadata.rdf").decode()
    sbml_name = next(n for n in z.namelist() if n.endswith(".sbml.xml"))
    minidom.parseString(z.read(sbml_name))  # round-trip: parses


def test_omex_is_deterministic(ds):
    a = combine.build_model_archive(ds[SCHNIDER], ds)
    b = combine.build_model_archive(ds[SCHNIDER], ds)
    assert a == b  # byte-identical: reproducibility guarantee
    da = combine.build_dataset_archive(ds)
    db = combine.build_dataset_archive(ds)
    assert da == db


def test_dataset_omex_bundles_models_and_readme(ds):
    z = _zip(combine.build_dataset_archive(ds))
    names = z.namelist()
    assert "README.txt" in names and "citations.bib" in names and "metadata.rdf" in names
    assert sum(n.endswith(".sbml.xml") for n in names) >= 2
    assert "NOT FOR CLINICAL USE" in z.read("README.txt").decode()


def test_bibtex_dataset_export(ds):
    bib = bibtex.build(ds)
    assert bib.count("@article{") == len(ds.citations)
    assert "doi = {10.1097/00000542-199805000-00006}" in bib  # Schnider DOI


def test_bibtex_for_model_includes_primary(ds):
    bib = bibtex.build_for_model(ds[SCHNIDER], ds)
    assert "@article{schnider-1998-propofol-pk" in bib


def test_eleveld_predictive_performance_backfilled(ds):
    perf = ds[ELEVELD].predictive_performance
    metrics = {p["metric"] for p in perf}
    assert {"MDPE", "MDAPE"} <= metrics
    for p in perf:
        assert p["citation"] == "eleveld-2018-propofol"
