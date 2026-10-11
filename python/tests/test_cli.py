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


def test_validate_cohort_self_consistency(capsys):
    # the known-answer fixture: a +20% offset must recover MDPE ≈ 20% from the CLI
    rc = main(["validate-cohort", "--model", "hypnotics_iv.propofol.eleveld_2018",
               "--self-consistency", "--offset", "20"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "external validation" in out and "MDPE" in out
    m = re.search(r"MDPE\s+([0-9.\-]+)%", out)
    assert m and float(m.group(1)) == pytest.approx(20.0, abs=0.01)


def test_validate_cohort_json_record(capsys):
    rc = main(["validate-cohort", "--model", "hypnotics_iv.propofol.eleveld_2018",
               "--self-consistency", "--offset", "10", "--json"])
    assert rc == 0
    rec = json.loads(capsys.readouterr().out)
    assert rec["reproducible"] is True and rec["mode"] == "pk_concentration"
    mdpe = next(m for m in rec["metrics"] if m["name"] == "MDPE")
    assert mdpe["value"] == pytest.approx(10.0, abs=0.01)


def test_validate_cohort_from_csv(tmp_path, capsys):
    csv = tmp_path / "cohort.csv"
    csv.write_text(
        "subject,time_min,observed,kind,age,weight,height,sex,bolus,infusion\n"
        "s1,5,3.1,cp,50,70,170,M,2 mg/kg,6 mg/kg/h\n"
        "s1,15,2.4,cp,50,70,170,M,2 mg/kg,6 mg/kg/h\n"
        "s1,30,2.0,cp,50,70,170,M,2 mg/kg,6 mg/kg/h\n", encoding="utf-8")
    rc = main(["validate-cohort", "--model", "hypnotics_iv.propofol.schnider_1998",
               "--observations", str(csv)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "1 subject(s)" in out and "MDAPE" in out


def test_validate_cohort_requires_a_source(capsys):
    rc = main(["validate-cohort", "--model", "hypnotics_iv.propofol.eleveld_2018"])
    assert rc == 2
    assert "self-consistency" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# Bad input is one line on stderr and exit 2 — never a traceback
# --------------------------------------------------------------------------- #
_MARSH = "hypnotics_iv.propofol.marsh_1991"
_ELEVELD = "hypnotics_iv.propofol.eleveld_2018"

_BAD_INPUT = [
    (["simulate", "nope.x.y"], "unknown model id 'nope.x.y'"),
    (["simulate", "hypnotics_iv.propofol.marsh"], "did you mean hypnotics_iv.propofol.marsh_1991"),
    (["simulate", _MARSH, "--weight", "0"], "patient weight must be a finite number > 0"),
    (["simulate", _MARSH, "--weight", "-5"], "patient weight must be a finite number > 0"),
    (["simulate", "hypnotics_iv.propofol.schnider_1998", "--height", "0"], "patient height"),
    (["simulate", _ELEVELD, "--age", "-1"], "patient age must be a finite number >= 0"),
    (["simulate", "hypnotics_iv.propofol.schnider_1998", "--sex", "X"], "sex must be 'M' or 'F'"),
    (["simulate", _MARSH, "--n", "0"], "--n must be at least 2"),
    (["simulate", _MARSH, "--tmax", "0"], "--tmax must be > 0"),
    (["simulate", "opioids.fentanyl.shafer_1990"], "kernel is not implemented"),
    (["simulate", "volatiles.sevoflurane.mac"], "expects a PK model"),
    (["simulate", _ELEVELD, "--bands", "--samples", "0"], "samples >= 1"),
    (["simulate", _ELEVELD, "--bands", "--percentile", "150"], "--percentile expects 'lo,hi'"),
    (["simulate", _ELEVELD, "--bands", "--percentile", "95,5"], "0 <= lo < hi <= 100"),
    (["simulate", _ELEVELD, "--pd", "nope"], "unknown model id 'nope'"),
    (["simulate", _ELEVELD, "--covariate-band", "--weight-sd", "-3"], "SD must be >= 0"),
    (["simulate", _MARSH, "--developmental"], "pma_weeks"),
    (["compare", "--drug", "nope"], "no PK models for drug 'nope'"),
    (["compare", "--drug", "propofol", "--weight", "0"], "patient weight"),
    (["tpeak", "nope"], "unknown model id"),
    (["tpeak", _MARSH, "--weight", "0"], "patient weight"),
    (["decrement", "nope"], "unknown model id"),
    (["mac", "--agent", "nope", "--age", "40"], "unknown model id 'volatiles.nope.mac'"),
    (["mac", "--agent", "sevoflurane", "--age", "-3"], "--age must be >= 0"),
    (["la", "--drug", "bupivacaine", "--dose", "0"], "--dose must be > 0"),
    (["la", "--drug", "bupivacaine", "--tmax", "0"], "--tmax must be > 0"),
    (["reverse", "--agent", "nope"], "unknown reversal agent 'nope' (known: naloxone, neostigmine, sugammadex)"),
    (["interact", "--propofol", "banana"], "unknown model id 'banana'"),
    (["interact", "--surface", "nope"], "unknown model id 'nope'"),
    (["covariates", "--equation", "james_1976", "--height", "0"], "patient height"),
    (["developmental", "--drug", "nope"], "no executable PK models for drug 'nope'"),
    (["export", "--format", "nonmem", "--output", "/dev/null/x"], "Not a directory"),
    (["export", "--format", "nonmem", "--output", "unused", "--model", "nope"], "unknown model id"),
]


@pytest.mark.parametrize("argv,needle", _BAD_INPUT, ids=[" ".join(a) for a, _ in _BAD_INPUT])
def test_bad_input_is_a_one_line_error_not_a_traceback(argv, needle, capsys, monkeypatch):
    monkeypatch.delenv("HYPNOS_DEBUG", raising=False)
    assert main(argv) == 2                      # returns; does not raise
    err = capsys.readouterr().err
    assert needle in err
    assert "Traceback" not in err


def test_hypnos_debug_reraises_for_a_real_traceback(monkeypatch):
    monkeypatch.setenv("HYPNOS_DEBUG", "1")
    with pytest.raises(KeyError):
        main(["simulate", "nope.x.y"])


@pytest.mark.parametrize("agent", ["sugammadex", "naloxone", "neostigmine"])
def test_reverse_accepts_the_short_agent_name_its_help_advertises(agent, capsys):
    assert main(["reverse", "--agent", agent]) == 0
    out = capsys.readouterr().out
    assert f"reversal: reversal.{agent}." in out


def test_reverse_still_accepts_a_full_model_id(capsys):
    assert main(["reverse", "--agent", "reversal.sugammadex.encapsulation"]) == 0
    assert "encapsulation" in capsys.readouterr().out
