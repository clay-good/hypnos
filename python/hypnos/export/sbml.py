"""SBML L3v2 exporter.

A PK model *is* a compartmental ODE system, so it maps cleanly onto SBML and
keeps continuity with the systems-biology toolchain (COPASI, Tellurium) and
with Nidus. The exported model carries the instantiated micro-rate constants as
global parameters (recoverable for round-trip validation) and expresses the
dynamics as SBML rate rules. A MIRIAM/``hypnos:`` RDF annotation rides in the
model ``<annotation>``.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from . import annotate
from ._common import resolve_patient, safe_name
from .registry import KERNELS

_HEADER = '<?xml version="1.0" encoding="UTF-8"?>\n'


def build(model, ds=None, patient: Optional[Dict[str, Any]] = None) -> str:
    pat = resolve_patient(model, patient)
    if not (model.kernel_implemented and model.kernel_function in KERNELS):
        # Emit a stub SBML carrying provenance only; no dynamics without a kernel.
        return (
            _HEADER
            + '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">\n'
            + f'  <model id="{safe_name(model)}" name="{model.label}">\n'
            + annotate.rdf_annotation_xml(model, ds)
            + "\n  </model>\n</sbml>\n"
        )

    p = KERNELS[model.kernel_function](pat)
    vc = p.as_volumes_clearances()
    name = safe_name(model)

    def param(pid, val, units="dimensionless"):
        return f'      <parameter id="{pid}" value="{val:.10g}" constant="true" units="{units}"/>'

    params = "\n".join([
        param("k10", p.k10, "per_min"),
        param("k12", p.k12, "per_min"),
        param("k21", p.k21, "per_min"),
        param("k13", p.k13, "per_min"),
        param("k31", p.k31, "per_min"),
        param("ke0", p.ke0, "per_min"),
        param("V1", vc["V1"], "litre"),
        param("V2", vc["V2"], "litre"),
        param("V3", vc["V3"], "litre"),
    ])

    # state variables as parameters with rate rules (amounts A1..A3, effect Ce)
    states = "\n".join([
        param("A1", 0.0, "mg").replace('constant="true"', 'constant="false"'),
        param("A2", 0.0, "mg").replace('constant="true"', 'constant="false"'),
        param("A3", 0.0, "mg").replace('constant="true"', 'constant="false"'),
        param("Ce", 0.0, "mg_per_l").replace('constant="true"', 'constant="false"'),
    ])

    def rate_rule(var, mathml):
        return (
            f'      <rateRule variable="{var}">\n'
            f'        <math xmlns="http://www.w3.org/1998/Math/MathML">\n{mathml}\n'
            f'        </math>\n      </rateRule>'
        )

    # dA1/dt = -(k10+k12+k13)A1 + k21 A2 + k31 A3
    a1 = (
        "          <apply><plus/>\n"
        "            <apply><times/><apply><minus/><apply><plus/><ci>k10</ci><ci>k12</ci><ci>k13</ci></apply></apply><ci>A1</ci></apply>\n"
        "            <apply><times/><ci>k21</ci><ci>A2</ci></apply>\n"
        "            <apply><times/><ci>k31</ci><ci>A3</ci></apply>\n"
        "          </apply>"
    )
    a2 = (
        "          <apply><minus/>\n"
        "            <apply><times/><ci>k12</ci><ci>A1</ci></apply>\n"
        "            <apply><times/><ci>k21</ci><ci>A2</ci></apply>\n"
        "          </apply>"
    )
    a3 = (
        "          <apply><minus/>\n"
        "            <apply><times/><ci>k13</ci><ci>A1</ci></apply>\n"
        "            <apply><times/><ci>k31</ci><ci>A3</ci></apply>\n"
        "          </apply>"
    )
    ce = (
        "          <apply><times/><ci>ke0</ci>\n"
        "            <apply><minus/><apply><divide/><ci>A1</ci><ci>V1</ci></apply><ci>Ce</ci></apply>\n"
        "          </apply>"
    )
    rules = "\n".join([rate_rule("A1", a1), rate_rule("A2", a2), rate_rule("A3", a3), rate_rule("Ce", ce)])

    units = (
        '    <listOfUnitDefinitions>\n'
        '      <unitDefinition id="per_min"><listOfUnits>'
        '<unit kind="second" exponent="-1" scale="0" multiplier="60"/></listOfUnits></unitDefinition>\n'
        '      <unitDefinition id="mg_per_l"><listOfUnits>'
        '<unit kind="gram" exponent="1" scale="-3" multiplier="1"/>'
        '<unit kind="litre" exponent="-1" scale="0" multiplier="1"/></listOfUnits></unitDefinition>\n'
        '    </listOfUnitDefinitions>'
    )

    return (
        _HEADER
        + '<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">\n'
        + f'  <model id="{name}" name="{model.label}">\n'
        + annotate.rdf_annotation_xml(model, ds) + "\n"
        + units + "\n"
        + "    <listOfParameters>\n" + params + "\n" + states + "\n    </listOfParameters>\n"
        + "    <listOfRules>\n" + rules + "\n    </listOfRules>\n"
        + "  </model>\n</sbml>\n"
    )


def filename(model) -> str:
    return f"{safe_name(model)}.sbml.xml"
