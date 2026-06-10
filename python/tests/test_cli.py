import json
import re

import pytest

from hypnos.cli import main


def _band_after(text, marker):
    """Pull the first ``[lo, hi]`` numeric pair that follows ``marker`` in ``text``."""
    seg = text.split(marker, 1)[1]
    m = re.search(r"\[\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\]", seg)
    assert m, f"no band found after {marker!r}"
    return float(m.group(1)), float(m.group(2))


def test_version(capsys):
    assert main(["version"]) == 0
    out = capsys.readouterr().out
    assert "hypnos" in out and "PROHIBITED" in out


def test_validate(capsys):
    assert main(["validate"]) == 0
    assert "OK" in capsys.readouterr().out


def test_info_is_json(capsys):
    assert main(["info"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["n_models"] >= 4


def test_simulate_cli(capsys):
    rc = main(["simulate", "hypnotics_iv.propofol.schnider_1998",
               "--age", "72", "--weight", "60", "--height", "162", "--sex", "F",
               "--pd", "pd_effect.propofol.bis_sigmoid"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "tier (propagated): C" in out


def test_compare_cli(capsys):
    rc = main(["compare", "--drug", "propofol",
               "--age", "40", "--weight", "140", "--height", "172", "--sex", "M"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "excluded for envelope" in out


def test_compare_bands_cli(capsys):
    rc = main(["compare", "--drug", "propofol", "--age", "72", "--weight", "60",
               "--height", "162", "--sex", "F", "--bands", "--samples", "200", "--seed", "7"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "band-tier B" in out          # Eleveld earns a band
    assert "variance share" in out       # variance decomposition surfaced
    assert "excluded from band metrics" in out  # Marsh/Schnider named, not fabricated


def test_compare_band_matches_simulate_at_peak(capsys):
    """compare and simulate report the SAME prediction band for the same model/seed.

    Regression guard: compare must report the band at the median-peak *instant*
    (q[lo][i], q[hi][i]), not the independent temporal maxima of each percentile
    curve — those peak at different times and yield an incoherent interval.
    """
    args = ["--age", "72", "--weight", "60", "--height", "162", "--sex", "F",
            "--bands", "--samples", "500", "--seed", "7"]
    assert main(["simulate", "hypnotics_iv.propofol.eleveld_2018", *args]) == 0
    sim_lo, sim_hi = _band_after(capsys.readouterr().out, "band-tier B  Ce peak")

    assert main(["compare", "--drug", "propofol", *args]) == 0
    cmp_lo, cmp_hi = _band_after(capsys.readouterr().out, "eleveld_2018")

    # compare prints 2dp, simulate 3dp; agree to compare's precision
    assert cmp_lo == pytest.approx(sim_lo, abs=0.01)
    assert cmp_hi == pytest.approx(sim_hi, abs=0.01)


def test_export_cli(tmp_path, capsys):
    rc = main(["export", "--format", "sbml", "--output", str(tmp_path),
               "--model", "hypnotics_iv.propofol.schnider_1998"])
    assert rc == 0
    files = list(tmp_path.glob("*.sbml.xml"))
    assert len(files) == 1
    assert "PROHIBITED" in files[0].read_text()
