"""TCI-sim JSON exporter.

A clean, documented JSON the open-TCI / simulator community can ingest:
population parameters + covariate equations + applicability envelope + the
mandatory NOT-FOR-CLINICAL-USE flag. Self-describing and round-trippable
against the reference kernel.
"""
from __future__ import annotations

import json
from typing import Any, Dict, Optional

from . import annotate
from ._common import resolve_patient, safe_name
from .registry import KERNELS


def build_dict(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    pat = resolve_patient(model, patient)
    doc: Dict[str, Any] = {
        "format": "hypnos.tci-sim/v1",
        "model_id": model.id,
        "label": model.label,
        "drug": model.drug,
        "purpose": model.purpose,
        "structure": model.structure,
        "covariate_equations": {p.symbol: p.covariate_model for p in model.parameters},
        "applicability_envelope": model.raw.get("applicability_envelope", {}),
        "known_failure_modes": model.raw.get("known_failure_modes", []),
        "instantiated_for": pat,
        "provenance": annotate.provenance(model, ds),
    }
    # v0.2 population-variability layer — lossless passthrough (it is JSON). A
    # consumer can pin reproducibility AND know how much of the NLME object is
    # curated; the never-synthesize rule means absence is a true "not curated".
    if model.has_published_variability:
        doc["variability"] = {
            "variability_status": model.variability_status,
            "band_tier": model.band_tier,
            "bsv": {p.symbol: {"omega2": p.variability.omega2, "cv_percent": p.variability.cv_percent}
                    for p in model.parameters if p.variability and p.variability.omega2 is not None},
            "residual_error": model.raw.get("residual_error"),
            "omega_block": model.raw.get("omega_block"),
        }
    else:
        doc["variability"] = {"variability_status": "none",
                              "note": "no published between-subject variability (no band)"}
    # v0.6 LA: the site-absorption block and the toxicity-threshold RANGES pass
    # through verbatim (TCI-JSON is lossless; v0.6 §9). Threshold ranges export AS
    # ranges — there is no projection that collapses them to a single value. The
    # drug-level protein binding rides along so the total->free story is consumable.
    if model.absorption is not None:
        doc["absorption"] = model.raw.get("absorption")
    if model.has_toxicity_thresholds:
        doc["toxicity_thresholds"] = model.raw.get("toxicity_thresholds")
        doc["safety_critical"] = True
        if ds is not None:
            drug = ds.drug(model.drug_name) or {}
            if drug.get("protein_binding") is not None:
                doc["protein_binding"] = drug["protein_binding"]
            if drug.get("cardiotoxicity_class") is not None:
                doc["cardiotoxicity_class"] = drug["cardiotoxicity_class"]   # v0.6 LA2
    if model.kernel_implemented and model.kernel_function in KERNELS:
        params = KERNELS[model.kernel_function](pat)
        doc["instantiated_parameters"] = {
            "micro_rate_constants": {
                "V1": params.V1, "k10": params.k10, "k12": params.k12,
                "k21": params.k21, "k13": params.k13, "k31": params.k31, "ke0": params.ke0,
            },
            "volumes_clearances": params.as_volumes_clearances(),
            "units": {"volume": "L", "clearance": "L/min", "rate_constant": "1/min"},
        }
    else:
        doc["instantiated_parameters"] = None
        doc["kernel_status"] = "pending — reference kernel not yet implemented/verified"
    return doc


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    return json.dumps(build_dict(model, ds, patient), indent=2) + "\n"


def filename(model) -> str:
    return f"{safe_name(model)}.tci.json"
