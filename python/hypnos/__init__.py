"""Hypnos — a curated, citation-backed, tier-annotated dataset of anesthetic and
perioperative PK/PD models, with envelope-aware forward simulation and exports
into the standard pharmacometric formats.

NOT a clinical decision-support tool. NOT a TCI pump driver. NOT a dosing
calculator. For research, method development, education, and simulation only.
See spec §10.
"""
from __future__ import annotations

__version__ = "0.1.0"

CLINICAL_USE = "PROHIBITED — research/education/simulation only"

from .analysis import DecrementTime, PeakEffect, decrement_time, time_to_peak_effect  # noqa: E402
from .filter import performance_table, select, summary  # noqa: E402
from .inhalational import (  # noqa: E402
    MacResult,
    WashinResult,
    WashoutResult,
    mac,
    washin,
    washin_comparison,
    washout,
    washout_comparison,
)
from .load import Dataset, load  # noqa: E402
from .models import Model, worst_tier  # noqa: E402
from .simulate import (  # noqa: E402
    Comparison,
    InteractionResult,
    SimulationResult,
    compare,
    simulate,
    simulate_interaction,
)
from .validate import assert_valid, validate_dataset  # noqa: E402
from .verification import (  # noqa: E402
    ModelVerification,
    checklist_markdown,
    model_verification,
    next_to_verify,
    verification_summary,
)

__all__ = [
    "__version__",
    "CLINICAL_USE",
    "load",
    "Dataset",
    "Model",
    "select",
    "summary",
    "performance_table",
    "simulate",
    "compare",
    "simulate_interaction",
    "Comparison",
    "SimulationResult",
    "InteractionResult",
    "mac",
    "MacResult",
    "washin",
    "washin_comparison",
    "WashinResult",
    "washout",
    "washout_comparison",
    "WashoutResult",
    "time_to_peak_effect",
    "PeakEffect",
    "decrement_time",
    "DecrementTime",
    "validate_dataset",
    "assert_valid",
    "worst_tier",
    "verification_summary",
    "model_verification",
    "next_to_verify",
    "checklist_markdown",
    "ModelVerification",
]
