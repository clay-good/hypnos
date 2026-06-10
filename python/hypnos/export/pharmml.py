"""PharmML exporter (projection).

PharmML is the standardized pharmacometric markup language — the "SBML of PK/PD"
and the durable interop anchor. A full PharmML document is large; Hypnos emits a
well-formed, namespaced PharmML projection that carries the structural model
(compartment count + parameterization), the instantiated population parameters,
the verbatim covariate equations, and the Hypnos provenance/clinical-use
annotation. It is round-trippable against the reference kernel via its parameter
block.
"""
from __future__ import annotations

from typing import Any, Dict, Optional
from xml.sax.saxutils import escape

from . import annotate
from ._common import resolve_patient, safe_name
from ._variability import omega_correlations, omega_diagonal, residual_spec
from .registry import KERNELS

_NS = 'xmlns:pharmml="http://www.pharmml.org/pharmml/0.8/PharmML"'


def _variability_xml(model) -> str:
    """First-class NLME random effects (spec §8): η → RandomEffect / VariabilityLevel,
    ε → the residual-error model. Empty when the model publishes no variability."""
    if not model.has_published_variability:
        return ""
    parts = [
        f'    <VariabilityModel bandTier="{model.band_tier}" '
        f'variabilityStatus="{escape(model.variability_status)}">',
        '      <VariabilityLevel symbol="indiv" type="betweenSubject"/>',
    ]
    for sym, om2, cv in omega_diagonal(model):
        parts.append(
            f'      <RandomEffect symbol="eta_{escape(sym)}" parameter="{escape(sym)}" '
            f'level="indiv" distribution="Normal" transformation="log" '
            f'variance="{om2:.10g}" cvPercent="{cv:.4g}"/>'
        )
    for a, b, r in omega_correlations(model):
        parts.append(
            f'      <Correlation between="{escape(a)} {escape(b)}" coefficient="{r:.6g}"/>'
        )
    spec = residual_spec(model)
    if spec is not None:
        terms = ""
        if spec.log_sd is not None:
            terms += f' logSd="{spec.log_sd:.10g}"'
        if spec.prop_var is not None:
            terms += f' proportionalVariance="{spec.prop_var:.10g}"'
        if spec.add_sd is not None:
            terms += f' additiveSd="{spec.add_sd:.10g}"'
        parts.append(
            f'      <ResidualError model="{escape(spec.model)}" '
            f'description="{escape(spec.label)}"{terms}/>'
        )
    parts.append('    </VariabilityModel>')
    return "\n" + "\n".join(parts)


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    pat = resolve_patient(model, patient)
    name = safe_name(model)
    prov = annotate.provenance(model, ds)

    cov_lines = "\n".join(
        f'      <CovariateEquation symbol="{escape(p.symbol)}">{escape(p.covariate_model or "")}</CovariateEquation>'
        for p in model.parameters
    )

    if model.kernel_implemented and model.kernel_function in KERNELS:
        vc = KERNELS[model.kernel_function](pat).as_volumes_clearances()
        param_lines = "\n".join(
            f'      <PopulationParameter symbol="{k}" value="{v:.10g}"/>'
            for k, v in vc.items()
        )
        kernel_note = ""
    else:
        param_lines = "      <!-- kernel pending: no instantiated parameters -->"
        kernel_note = '\n    <KernelStatus>pending</KernelStatus>'

    safety_xml = (
        '      <hypnos:safetyCritical>true</hypnos:safetyCritical>\n'
        if "hypnos:safetyCritical" in prov else ""
    )
    # v0.6 §9: the LA toxicity-threshold RANGES travel as annotation, AS ranges.
    la_xml = "".join(
        f'      <hypnos:toxicityThresholdRange endpoint="{escape(th.endpoint)}" '
        f'basis="{escape(th.basis)}" low="{th.low:.10g}" high="{th.high:.10g}" '
        f'units="{escape(th.units)}" tier="{th.tier}"/>\n'
        for th in model.toxicity_thresholds
    )
    prov_xml = (
        "    <Annotation>\n"
        f'      <hypnos:clinicalUse>{escape(prov["hypnos:clinicalUse"])}</hypnos:clinicalUse>\n'
        f'      <hypnos:confidenceTier>{prov["hypnos:confidenceTier"]}</hypnos:confidenceTier>\n'
        f'      <hypnos:datasetVersion>{prov["hypnos:datasetVersion"]}</hypnos:datasetVersion>\n'
        + safety_xml
        + la_xml
        + "".join(
            f'      <bqmodel:isDerivedFrom>{escape(u)}</bqmodel:isDerivedFrom>\n'
            for u in prov["bqmodel:isDerivedFrom"]
        )
        + "    </Annotation>"
    )

    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<pharmml:PharmML {_NS} '
        'xmlns:hypnos="https://w3id.org/hypnos/terms#" '
        'xmlns:bqmodel="http://biomodels.net/model-qualifiers/" '
        'writtenVersion="0.8" implementedBy="hypnos">\n'
        f'  <pharmml:ModelDefinition id="{name}" name="{escape(model.label)}">\n'
        f'    <StructuralModel compartments="{model.n_compartments}" '
        f'parameterization="{model.structure["parameterization"]}" '
        f'effectCompartment="{str(model.has_effect_compartment).lower()}"/>\n'
        f'    <InstantiatedFor age="{pat.get("age")}" weight="{pat.get("weight")}" '
        f'height="{pat.get("height")}" sex="{pat.get("sex")}"/>\n'
        '    <PopulationParameters>\n' + param_lines + '\n    </PopulationParameters>\n'
        '    <CovariateModel>\n' + cov_lines + '\n    </CovariateModel>'
        + _variability_xml(model)
        + kernel_note + '\n'
        + prov_xml + '\n'
        '  </pharmml:ModelDefinition>\n'
        '</pharmml:PharmML>\n'
    )


def filename(model) -> str:
    return f"{safe_name(model)}.pharmml.xml"
