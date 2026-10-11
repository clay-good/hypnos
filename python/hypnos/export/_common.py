"""Shared helpers for exporters."""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from ..models import Model, worst_tier

# A canonical adult reference patient. Exports instantiate the (population)
# covariate model at this individual so the emitted artifact is a concrete,
# simulatable model. Every export states the patient it was instantiated for.
REFERENCE_PATIENT: Dict[str, Any] = {"age": 50, "weight": 77, "height": 177, "sex": "M"}


def resolve_patient(model: Model, patient: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if patient is not None:
        from ..covariates import point_patient
        point_patient(patient)           # refuses a physically impossible patient (ValueError)
        return dict(patient)
    return dict(REFERENCE_PATIENT)


def instantiated_standing(model: Model, ds, patient: Dict[str, Any]) -> Tuple[str, List[str]]:
    """``(propagated tier, warnings)`` for the individual an export is instantiated at.

    An export's numbers are a model evaluated for ONE stated patient. If that patient is
    outside the model's applicability envelope (a pediatric model at the adult reference
    patient, Schnider at BMI 50, an organ-failure patient with no cited standing), the
    artifact must say so and carry the tiered-down label — the same worst-tier-wins rule
    ``simulate`` applies (v0.1 §5/§7). Uses the one shared safety evaluation, so an export
    can never be labeled more confidently than the simulation it describes."""
    from ..simulate import evaluate_safety

    drug = (ds.drug(model.drug_name) if ds is not None else None) or {}
    floor, warnings, _ = evaluate_safety(model, patient, drug)
    return worst_tier([model.tier, floor]), list(warnings)


def safe_name(model: Model) -> str:
    return model.id.replace(".", "_")
