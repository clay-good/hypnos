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

from .filter import select, summary  # noqa: E402
from .load import Dataset, load  # noqa: E402
from .models import Model, worst_tier  # noqa: E402
from .simulate import Comparison, SimulationResult, compare, simulate  # noqa: E402
from .validate import assert_valid, validate_dataset  # noqa: E402

__all__ = [
    "__version__",
    "CLINICAL_USE",
    "load",
    "Dataset",
    "Model",
    "select",
    "summary",
    "simulate",
    "compare",
    "Comparison",
    "SimulationResult",
    "validate_dataset",
    "assert_valid",
    "worst_tier",
]
