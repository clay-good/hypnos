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
        # v0.7 C3: the structured covariate_model block passes through VERBATIM (lossless) so a
        # downstream simulator knows WHICH named LBM/FFM equation to use — directly addressing
        # the documented pump-to-pump body-size substitution problem (§1/§8).
        "covariate_model": model.raw.get("covariate_model"),
        "covariate_sensitivity_status": model.covariate_sensitivity_status,
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
    # v0.3 estimation-uncertainty layer — verbatim passthrough, in its OWN block beside
    # `variability` (an estimation RSE is not a between-subject CV; v0.3 Trap 1). The SEs
    # are on the curated typical values, so those ride along with them.
    if model.has_published_estimation:
        doc["estimation"] = {
            "uncertainty_status": model.uncertainty_status,
            "estimation_band_tier": model.estimation_band_tier,
            "theta": {p["symbol"]: {"estimate": p["value"].get("central"),
                                    "estimation_uncertainty": p["estimation_uncertainty"]}
                      for p in model.raw.get("parameters", [])
                      if (p.get("estimation_uncertainty") or {}).get("se") is not None},
            "estimate_covariance": model.raw.get("estimate_covariance"),
        }
    else:
        doc["estimation"] = {"uncertainty_status": "none",
                             "note": "no published estimation uncertainty (no confidence band)"}
    # v0.4 §8: the Hypnos-COMPUTED external_validation[] blocks pass through verbatim, with
    # their manifest/provenance, so a consumer sees reproducible metrics — kept apart from
    # the publisher-reported predictive_performance.
    doc["validation_status"] = model.validation_status
    if model.has_external_validation:
        doc["external_validation"] = model.external_validation
    # v0.5 §D: the reversal coupling and the opt-in Tier-D disease-state modifiers pass
    # through verbatim. The instantiated parameters below are NEVER pre-adjusted by a
    # modifier — a consumer sees the adjustment, its citation, and its opt-in flag.
    if model.reversal_mechanism:
        doc["reversal_mechanism"] = model.raw["reversal_mechanism"]
    if model.raw.get("disease_state_modifiers"):
        doc["disease_state_modifiers"] = model.raw["disease_state_modifiers"]
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
