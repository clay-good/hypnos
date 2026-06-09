"""Shared helpers for exporters."""
from __future__ import annotations

from typing import Any, Dict, Optional

from ..models import Model

# A canonical adult reference patient. Exports instantiate the (population)
# covariate model at this individual so the emitted artifact is a concrete,
# simulatable model. Every export states the patient it was instantiated for.
REFERENCE_PATIENT: Dict[str, Any] = {"age": 50, "weight": 77, "height": 177, "sex": "M"}


def resolve_patient(model: Model, patient: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if patient is not None:
        return dict(patient)
    return dict(REFERENCE_PATIENT)


def safe_name(model: Model) -> str:
    return model.id.replace(".", "_")
