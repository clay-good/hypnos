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
from .registry import KERNELS

_NS = 'xmlns:pharmml="http://www.pharmml.org/pharmml/0.8/PharmML"'


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

    prov_xml = (
        "    <Annotation>\n"
        f'      <hypnos:clinicalUse>{escape(prov["hypnos:clinicalUse"])}</hypnos:clinicalUse>\n'
        f'      <hypnos:confidenceTier>{prov["hypnos:confidenceTier"]}</hypnos:confidenceTier>\n'
        f'      <hypnos:datasetVersion>{prov["hypnos:datasetVersion"]}</hypnos:datasetVersion>\n'
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
        + kernel_note + '\n'
        + prov_xml + '\n'
        '  </pharmml:ModelDefinition>\n'
        '</pharmml:PharmML>\n'
    )


def filename(model) -> str:
    return f"{safe_name(model)}.pharmml.xml"
