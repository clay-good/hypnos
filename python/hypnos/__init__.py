"""Hypnos — a curated, citation-backed, tier-annotated dataset of anesthetic and
perioperative PK/PD models, with envelope-aware forward simulation and exports
into the standard pharmacometric formats.

NOT a clinical decision-support tool. NOT a TCI pump driver. NOT a dosing
calculator. For research, method development, education, and simulation only.
See spec §10.
"""
from __future__ import annotations

__version__ = "0.6.0"

CLINICAL_USE = "PROHIBITED — research/education/simulation only"

from .analysis import (  # noqa: E402
    CohortValidation,
    DecrementTime,
    PeakEffect,
    PopulationPerformance,
    SubjectRecord,
    VarvelResult,
    decrement_time,
    performance_error,
    pooled_performance,
    subjects_from_cohort_self_consistency,
    subjects_from_csv,
    time_to_peak_effect,
    validate_against_cohort,
    varvel_metrics,
)
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
from .models import (  # noqa: E402
    EstimateCovariance,
    EstimationUncertainty,
    Model,
    OmegaBlock,
    Parameter,
    ParameterVariability,
    ResidualError,
    ToxicityThreshold,
    worst_tier,
)
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
    "Parameter",
    "ParameterVariability",
    "ResidualError",
    "OmegaBlock",
    "EstimationUncertainty",
    "EstimateCovariance",
    "ToxicityThreshold",
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
    "performance_error",
    "varvel_metrics",
    "pooled_performance",
    "validate_against_cohort",
    "subjects_from_csv",
    "subjects_from_cohort_self_consistency",
    "VarvelResult",
    "PopulationPerformance",
    "SubjectRecord",
    "CohortValidation",
    "validate_dataset",
    "assert_valid",
    "worst_tier",
    "verification_summary",
    "model_verification",
    "next_to_verify",
    "checklist_markdown",
    "ModelVerification",
]
