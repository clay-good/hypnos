"""Shared dosing presets, simulatable-drug discovery, and the `hypnos models` CLI."""
import pytest

import hypnos
from hypnos.cli import main
from hypnos.filter import pk_drugs
from hypnos.presets import DEFAULT_SCHEDULES, default_schedule_for


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_presets_are_drug_appropriate():
    assert default_schedule_for("propofol")[0][2] == "2 mg/kg"
    assert "mcg/kg" in default_schedule_for("remifentanil")[0][2]
    assert default_schedule_for("rocuronium")[0][2] == "0.6 mg/kg"
    # unknown drug falls back to the propofol default (never an empty schedule)
    assert default_schedule_for("nope") == DEFAULT_SCHEDULES["propofol"]


def test_cli_and_presets_share_one_source():
    # the CLI alias must be the shared preset function, not a private copy
    from hypnos import cli
    assert cli._default_schedule_for is default_schedule_for


def test_pk_drugs_lists_simulatable_only(ds):
    drugs = pk_drugs(ds)
    # drugs with >=1 executable PK kernel
    assert "propofol" in drugs and "remifentanil" in drugs and "dexmedetomidine" in drugs
    # fentanyl PK kernel is pending -> not simulatable; volatiles aren't PK
    assert "fentanyl" not in drugs and "sevoflurane" not in drugs
    assert drugs == sorted(drugs)


def test_cli_models_all(capsys):
    assert main(["models"]) == 0
    out = capsys.readouterr().out
    assert "simulatable drugs" in out
    assert "opioids.remifentanil.eleveld_2017" in out


def test_cli_models_by_drug(capsys):
    assert main(["models", "--drug", "remifentanil"]) == 0
    out = capsys.readouterr().out
    assert "remifentanil" in out and "propofol" not in out


def test_cli_models_unknown_drug():
    assert main(["models", "--drug", "nope"]) == 2
