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
from xml.sax.saxutils import escape, quoteattr

from . import annotate
from ._common import resolve_patient, safe_name
from ._estimation import (
    covariance_step_note, estimate_correlations, independence_caveat, theta_covariance,
    theta_estimates,
)
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


def _estimation_xml(model) -> str:
    """The estimation object in the PharmML Standard Output's shape (v0.3 §8).

    ``Estimation/PopulationEstimates/MLE`` carries the CURATED typical values the SEs are
    on (not the patient-instantiated ``PopulationParameters`` above);
    ``PrecisionPopulationEstimates`` carries ``StandardError`` / ``RelativeStandardError``
    and the estimate ``CovarianceMatrix`` + any published ``CorrelationMatrix`` entries.
    Empty when the model publishes no estimation uncertainty (never-synthesize)."""
    ests = theta_estimates(model)
    if not ests:
        return ""
    method = escape(ests[0].method)
    parts = [
        f'    <StandardOutput uncertaintyStatus="{escape(model.uncertainty_status)}" '
        f'estimationBandTier="{model.estimation_band_tier}">',
        '      <!-- ESTIMATION uncertainty (SE on the typical value; reducible with more data) '
        '— NOT the between-subject VariabilityModel above -->',
        '      <Estimation>',
        '        <PopulationEstimates>',
        '          <MLE>',
    ]
    parts += [f'            <Estimate symbol="{escape(e.symbol)}" value="{e.central:.10g}"/>'
              for e in ests]
    parts += ['          </MLE>', '        </PopulationEstimates>',
              f'        <PrecisionPopulationEstimates method="{method}">',
              '          <StandardError>']
    for e in ests:
        ci = f' ci95Low="{e.ci95[0]:.10g}" ci95High="{e.ci95[1]:.10g}"' if e.ci95 else ""
        parts.append(f'            <SE symbol="{escape(e.symbol)}" value="{e.se:.10g}" '
                     f'scale="{escape(e.scale)}" tier="{e.tier}"{ci}/>')
    parts += ['          </StandardError>', '          <RelativeStandardError>']
    parts += [f'            <RSE symbol="{escape(e.symbol)}" percent="{e.rse_percent:.6g}"/>'
              for e in ests if e.rse_percent is not None]
    parts.append('          </RelativeStandardError>')
    syms, cov = theta_covariance(model)
    parts.append(f'          <CovarianceMatrix scale="natural" rowNames="{escape(" ".join(syms))}">')
    for i, row in enumerate(cov):
        cells = " ".join(f"{row[j]:.10g}" for j in range(i + 1))
        parts.append(f'            <Row symbol="{escape(syms[i])}">{cells}</Row>')
    parts.append('          </CovarianceMatrix>')
    corrs = estimate_correlations(model)
    if corrs:
        parts.append('          <CorrelationMatrix>')
        parts += [f'            <Correlation between="{escape(a)} {escape(b)}" '
                  f'coefficient="{r:.6g}"/>' for a, b, r in corrs]
        parts.append('          </CorrelationMatrix>')
    parts.append('        </PrecisionPopulationEstimates>')
    for note in (covariance_step_note(model), independence_caveat(model)):
        if note:
            parts.append(f'        <Note>{escape(note)}</Note>')
    parts += ['      </Estimation>', '    </StandardOutput>']
    return "\n" + "\n".join(parts)


def _validation_xml(model) -> str:
    """Hypnos-computed external-validation metrics as a goodness-of-fit ``Diagnostics``
    block (v0.4 §8), a sibling of the v0.3 estimation object. Empty when none are curated."""
    entries = annotate._validation_entries(model)
    if not entries:
        return ""
    parts = [f'    <Diagnostics validationStatus="{escape(model.validation_status)}">',
             f'      <!-- {escape(annotate._NOT_CLEARANCE)} -->']
    for e in entries:
        attrs = (f'dataset={quoteattr(str(e["dataset"]))} mode="{escape(str(e["mode"]))}" '
                 f'target="{escape(str(e["target"]))}"')
        if e["n_subjects"] is not None:
            attrs += f' nSubjects="{e["n_subjects"]}"'
        if e["in_envelope"] is not None:
            attrs += f' inEnvelope="{str(e["in_envelope"]).lower()}"'
        parts.append(f'      <ExternalValidation {attrs}>')
        parts += [f'        <Metric name="{escape(str(n))}" value="{v:.10g}" units={quoteattr(u)}/>'
                  for n, v, u in e["metrics"]]
        parts.append('      </ExternalValidation>')
    parts.append('    </Diagnostics>')
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
    # v0.6 LA2: the drug's cardiotoxicity class (stereochemistry / CNS-to-CVS margin).
    cc = (ds.drug(model.drug_name) or {}).get("cardiotoxicity_class") if ds is not None else None
    if cc:
        la_xml += (
            f'      <hypnos:cardiotoxicityClass rank="{escape(cc.get("rank", ""))}" '
            f'stereochemistry="{escape(cc.get("stereochemistry", ""))}" '
            f'cnsToCvsMargin="{escape(cc.get("cns_to_cvs_margin", ""))}"/>\n'
        )
    ffm = ((ds.drug(model.drug_name) or {}).get("protein_binding") or {}).get(
        "free_fraction_model") if ds is not None else None
    if ffm:
        la_xml += (
            f'      <hypnos:freeFractionModel type="{escape(ffm.get("type", ""))}" '
            f'bindingCapacityUgMl="{ffm.get("binding_capacity_ug_ml")}" '
            f'tier="{escape(ffm.get("tier", ""))}"/>\n'
        )
    # v0.5 §D: reversal coupling + opt-in Tier-D disease modifiers (annotation, not parameters).
    rev = model.reversal_mechanism
    if rev:
        la_xml += (
            f'      <hypnos:reversalMechanism type="{escape(str(rev.get("type")))}" '
            f'tier="{escape(str(rev.get("tier", "D")))}" '
            f'targets="{escape(" ".join(model.reversal_targets))}"'
            + "".join(f" {annotate._camel(k)}={quoteattr(str(v))}"
                      for k, v in annotate._reversal_facts(rev))
            + "/>\n")
    for dz in model.disease_state_modifiers:
        sf = dz.scale_factor
        la_xml += (
            f'      <hypnos:diseaseStateModifier condition={quoteattr(annotate._condition_text(dz))} '
            f'affectedParameter="{escape(dz.affected_parameter)}"'
            + (f' scaleFactor="{sf:.10g}"' if sf is not None else "")
            + f' evidenceTier="{dz.evidence_tier}" '
            f'appliedByDefault="{str(dz.applied_by_default).lower()}"/>\n')
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
        + _estimation_xml(model)
        + _validation_xml(model)
        + kernel_note + '\n'
        + prov_xml + '\n'
        '  </pharmml:ModelDefinition>\n'
        '</pharmml:PharmML>\n'
    )


def filename(model) -> str:
    return f"{safe_name(model)}.pharmml.xml"
