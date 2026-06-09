"""Published predictive-performance surfacing (`performance_table` + `hypnos performance`).

The metrics are the numeric counterpart to the editorial tier (spec §5). These
tests assert the data is well-formed, traceable to a real citation, and surfaced —
not that any particular number is "good" (a Tier-D-regime MDAPE can be huge by design).
"""
import pytest

import hypnos
from hypnos.cli import main
from hypnos.filter import performance_table, summary

_ALLOWED_METRICS = {"MDPE", "MDAPE", "wobble", "divergence", "RMSE", "other"}


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_performance_table_rows_are_wellformed_and_cited(ds):
    rows = performance_table(ds)
    assert rows, "expected at least some backfilled predictive-performance metrics"
    for r in rows:
        assert r["metric"] in _ALLOWED_METRICS
        assert isinstance(r["value"], (int, float))
        assert r["model_id"] in ds
        # every metric must resolve to a real citation record with a DOI — the
        # honesty stance: a performance number is never asserted bare.
        assert r["citation"] is not None, f"{r['model_id']} {r['metric']} has no citation"
        assert ds.citation(r["citation"]) is not None, f"dangling citation {r['citation']!r}"
        assert r["doi"], f"{r['citation']} has no DOI"


def test_backfilled_models_are_present(ds):
    by_model = {}
    for r in performance_table(ds):
        by_model.setdefault(r["model_id"], set()).add(r["metric"])
    # Marsh gained its first metric (external head-to-head MDAPE)
    assert "MDAPE" in by_model.get("hypnotics_iv.propofol.marsh_1991", set())
    # Minto carries the quantified James-LBM-in-obesity failure mode (all four Varvel metrics)
    minto = by_model.get("opioids.remifentanil.minto_1997", set())
    assert {"MDPE", "MDAPE", "divergence", "wobble"} <= minto


def test_minto_obesity_failure_mode_is_quantified(ds):
    """La Colla 2010: the original James-LBM Minto set is clinically unacceptable
    in morbid obesity — the published number behind a documented failure mode."""
    rows = [r for r in performance_table(ds, drug="remifentanil")
            if r["metric"] == "MDAPE" and r["citation"] == "lacolla-2010-minto-obese"]
    assert len(rows) == 1
    assert rows[0]["value"] == pytest.approx(53.4)
    assert "obese" in (rows[0]["population"] or "")


def test_performance_drug_filter(ds):
    assert all(r["model_id"].split(".")[1] == "propofol"
               for r in performance_table(ds, drug="propofol"))


def test_summary_counts_models_with_performance(ds):
    s = summary(ds)
    n = s["models_with_predictive_performance"]
    assert n == sum(1 for m in ds if m.predictive_performance)
    assert n >= 5


def test_cli_performance_all(capsys):
    assert main(["performance"]) == 0
    out = capsys.readouterr().out
    assert "MDAPE" in out and "marsh_1991" in out
    assert "lacolla-2010-minto-obese" in out


def test_cli_performance_by_drug(capsys):
    assert main(["performance", "--drug", "propofol"]) == 0
    out = capsys.readouterr().out
    assert "propofol" in out and "remifentanil" not in out


def test_cli_performance_no_metrics_for_drug():
    # rocuronium has no published performance metrics in the dataset
    assert main(["performance", "--drug", "rocuronium"]) == 2


def test_predictive_mdape_excludes_out_of_envelope(ds):
    # Minto carries an in-envelope MDAPE (24.6%) AND an out-of-envelope one (53.4%,
    # morbid obesity). For an *included* (in-envelope) model only the former applies.
    vals = [e["value"] for e in ds["opioids.remifentanil.minto_1997"].predictive_mdape]
    assert 24.6 in vals and 53.4 not in vals
    # a model with no MDAPE returns an empty list
    assert ds["opioids.remifentanil.kim_2017"].predictive_mdape == []
    # Marsh has exactly one (external head-to-head)
    assert [e["value"] for e in ds["hypnotics_iv.propofol.marsh_1991"].predictive_mdape] == [25.0]


def test_compare_cli_shows_accuracy(capsys):
    # the divergence view now reports each included model's published inaccuracy
    assert main(["compare", "--drug", "propofol", "--age", "50", "--weight", "80",
                 "--height", "178", "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "MDAPE" in out
    # Minto's in-envelope number is shown for an in-envelope remifentanil patient,
    # never the morbid-obesity failure-mode value
    assert main(["compare", "--drug", "remifentanil", "--age", "50", "--weight", "80",
                 "--height", "178", "--sex", "M"]) == 0
    out = capsys.readouterr().out
    assert "24.6" in out and "53.4" not in out
