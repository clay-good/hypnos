"""Smoke test for the figure regeneration in scripts/regenerate.py.

The figures are deterministic projections of the dataset/kernels — the same
guarantee as the exports — so they must regenerate from the live data without
hand-editing. This exercises every figure generator (catching API drift that
would otherwise silently rot a README figure). matplotlib lives in the
``notebooks`` extra, so the test skips where it is not installed.
"""
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

import hypnos

ROOT = Path(__file__).resolve().parents[2]
_REGEN = ROOT / "scripts" / "regenerate.py"


def _load_regen():
    spec = importlib.util.spec_from_file_location("hypnos_regenerate", _REGEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_all_figures_regenerate(tmp_path, monkeypatch):
    import matplotlib
    matplotlib.use("Agg")

    regen = _load_regen()
    monkeypatch.setattr(regen, "IMAGES", tmp_path)   # never clobber the committed figures
    ds = hypnos.load()
    assert regen.regenerate_figures(ds) is True

    expected = {"divergence.png", "synergy.png", "pediatric.png", "mac_age.png",
                "washin.png", "washout.png", "variability.png", "effect_band.png",
                "la_double_uncertainty.png", "la_cardiotoxicity.png"}
    produced = {p.name for p in tmp_path.glob("*.png")}
    assert expected <= produced
    # every figure is a non-trivial PNG (renders something, not an empty file)
    for name in expected:
        assert (tmp_path / name).stat().st_size > 5000
