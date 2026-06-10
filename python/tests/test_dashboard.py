"""Smoke test for the Streamlit dashboard via Streamlit's AppTest harness.

The dashboard is a thin presentation layer over tested package logic, so the
risk it adds is purely in the rendering plumbing (Altair ribbons, the v0.2
band/separation/variance readout). This test runs the *whole* script end to end
— with and without the prediction-band toggle — and asserts it renders without
raising.

It is skipped wherever the dashboard extra (streamlit/altair/pandas) is not
installed — notably the CI ``test`` job, which installs only ``.[dev]``. Run it
locally with ``pip install -e '.[dashboard,dev]'``.
"""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
pytest.importorskip("altair")
pytest.importorskip("pandas")

from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parents[2] / "dashboard" / "app.py")


def test_dashboard_renders_without_bands():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception
    # the divergence view is the headline; it must produce at least one chart
    assert at.title  # page title rendered


def test_dashboard_renders_with_prediction_bands():
    at = AppTest.from_file(APP, default_timeout=120).run()
    assert not at.exception
    bands_toggle = next(c for c in at.checkbox if c.label.startswith("Seeded"))
    bands_toggle.set_value(True).run()
    assert not at.exception
    # the uncertainty-aware divergence section must render once bands are on
    assert any("distinguishable" in s.value for s in at.subheader)
    # default drug is propofol → Eleveld is band-eligible, Marsh/Schnider are not,
    # so the never-synthesize "excluded from band math" warning must appear.
    warnings = " ".join(w.value for w in at.warning)
    assert "excluded from band math" in warnings
    # the PD effect (BIS) band — PK BSV through the Hill link — must surface too
    # (Eleveld PK + Eleveld BIS), as an honest lower bound.
    subheaders = " ".join(s.value for s in at.subheader)
    assert "PD effect band" in subheaders
    captions = " ".join(c.value for c in at.caption)
    assert "lower bound" in captions
