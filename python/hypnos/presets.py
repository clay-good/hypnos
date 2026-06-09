"""Drug-appropriate default dose schedules (single source of truth).

A propofol regimen (2 mg/kg + 6 mg/kg/h) applied to remifentanil would be a
~1000x overdose, so each drug gets a clinically sensible default. Both the CLI
and the dashboard use these, so the two never drift. Doses are illustrative
defaults for *simulation*, not recommendations (Hypnos emits no dose advice).
"""
from __future__ import annotations

from typing import List, Tuple

Schedule = List[Tuple[str, float, str]]

_PROPOFOL_DEFAULT: Schedule = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]

DEFAULT_SCHEDULES = {
    "propofol": _PROPOFOL_DEFAULT,
    "remifentanil": [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")],
    "dexmedetomidine": [("infusion", 0.0, "6 mcg/kg/h"), ("infusion", 10.0, "0.5 mcg/kg/h")],
    "fentanyl": [("bolus", 0.0, "2 mcg/kg")],
    "rocuronium": [("bolus", 0.0, "0.6 mg/kg")],
}


def default_schedule_for(drug: str) -> Schedule:
    """A clinically sensible default dosing schedule for the given drug."""
    return DEFAULT_SCHEDULES.get(drug, _PROPOFOL_DEFAULT)
