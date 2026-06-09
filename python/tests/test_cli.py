import json

import pytest

from hypnos.cli import main


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


def test_export_cli(tmp_path, capsys):
    rc = main(["export", "--format", "sbml", "--output", str(tmp_path),
               "--model", "hypnotics_iv.propofol.schnider_1998"])
    assert rc == 0
    files = list(tmp_path.glob("*.sbml.xml"))
    assert len(files) == 1
    assert "PROHIBITED" in files[0].read_text()
