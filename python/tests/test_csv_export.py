import csv
import io

import pytest

import hypnos
from hypnos.export import FORMATS, csv_flat, export_model

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_csv_in_formats():
    assert "csv" in FORMATS


def test_dataset_csv_one_row_per_parameter(ds):
    rows = list(csv.DictReader(io.StringIO(csv_flat.build(ds))))
    expected = sum(len(m.parameters) for m in ds.models)
    assert len(rows) == expected
    assert set(rows[0].keys()) >= {"model_id", "symbol", "central", "units",
                                   "param_tier", "covariate_model", "doi"}


def test_csv_resolves_doi_and_quotes_commas(ds):
    rows = list(csv.DictReader(io.StringIO(csv_flat.build(ds))))
    cl1 = next(r for r in rows if r["model_id"] == SCHNIDER and r["symbol"] == "Cl1")
    # DOI joined from the citation record
    assert cl1["doi"] == "10.1097/00000542-199805000-00006"
    assert cl1["param_tier"] == "B" and cl1["units"] == "L/min"
    # the covariate equation contains commas; CSV quoting must round-trip it intact
    assert cl1["covariate_model"] == "Cl1 = 1.89 + 0.0456*(WGT-77) - 0.0681*(LBM-59) + 0.0264*(HGT-177)"


def test_per_model_csv_is_subset(ds):
    rows = list(csv.DictReader(io.StringIO(csv_flat.build_for_model(ds[SCHNIDER], ds))))
    assert rows and all(r["model_id"] == SCHNIDER for r in rows)
    assert len(rows) == len(ds[SCHNIDER].parameters)


def test_export_model_csv(ds):
    fname, text = export_model("csv", ds[SCHNIDER], ds)
    assert fname == "parameters.csv"
    assert text.startswith("model_id,")


def test_kernel_pending_param_has_blank_central(ds):
    # rocuronium V1 has central=null (kernel pending) -> blank cell, not "None"
    # (fentanyl's params were source-filled, so rocuronium is now the null-valued example)
    rows = list(csv.DictReader(io.StringIO(csv_flat.build(ds))))
    roc = [r for r in rows if r["model_id"] == "nmb_agents.rocuronium.wierda_1991" and r["symbol"] == "V1"]
    assert roc and roc[0]["central"] == ""
