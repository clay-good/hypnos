"""Provenance & safety annotations attached to every export.

Two durability tiers, mirroring the spec (§7):

* **MIRIAM-style RDF** (``bqmodel:isDerivedFrom`` -> DOI/PMID) survives even if a
  downstream tool strips the custom ``hypnos:`` predicates;
* **custom ``hypnos:`` predicates** carry the propagated confidence tier, the
  dataset version, and the universal, machine-readable
  ``hypnos:clinicalUse = "PROHIBITED"`` flag.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional
from xml.sax.saxutils import quoteattr

from .. import CLINICAL_USE, __version__
from ..models import Model

HYPNOS_NS = "https://w3id.org/hypnos/terms#"


def citation_uris(model: Model, ds) -> List[str]:
    """Return resolvable DOI/PMID URIs for a model's primary citation."""
    uris: List[str] = []
    cit = ds.citation(model.primary_citation) if ds is not None else None
    if cit:
        if cit.get("doi"):
            uris.append(f"https://doi.org/{cit['doi']}")
        if cit.get("pmid"):
            uris.append(f"https://identifiers.org/pubmed:{cit['pmid']}")
    return uris


def provenance(model: Model, ds=None, tier: Optional[str] = None) -> Dict[str, Any]:
    """Structured provenance block reused across all exporters."""
    prov = {
        "hypnos:datasetVersion": __version__,
        "hypnos:modelId": model.id,
        "hypnos:confidenceTier": tier or model.tier,
        "hypnos:reviewStatus": model.review_status,
        "hypnos:clinicalUse": CLINICAL_USE,
        "bqmodel:isDerivedFrom": citation_uris(model, ds),
    }
    # Local-anesthetic records are doubly safety-critical (v0.6 §7/§9): the
    # concentration-vs-threshold view is exactly the artifact that tempts the
    # forbidden "max safe dose?" question, so a downstream consumer is twice-warned.
    if model.is_safety_critical:
        prov["hypnos:safetyCritical"] = "true"
    return prov


def banner(model: Model, tier: Optional[str] = None) -> str:
    """Human-readable banner for text-format exports (NONMEM, R, etc.)."""
    t = tier or model.tier
    safety = (
        "  SAFETY-CRITICAL (local anesthetic): systemic concentration only — NOT a\n"
        "  toxicity margin, max-dose, or 'is this safe?' tool. Thresholds are RANGES.\n"
        if model.is_safety_critical else ""
    )
    return (
        "============================================================\n"
        "  HYPNOS EXPORT — NOT FOR CLINICAL USE\n"
        f"  clinicalUse = {CLINICAL_USE}\n"
        f"  model: {model.id}  (tier {t}, review: {model.review_status})\n"
        f"  dataset version: {__version__}\n"
        + safety
        + "  Research / education / simulation only. Not a dosing tool.\n"
        "============================================================"
    )


def variability_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates describing the curated random-effects layer (v0.2 §8).

    SBML core cannot express population random effects, so the Ω/Σ ride as annotation
    metadata: a deterministic SBML consumer (COPASI/Tellurium) still sees only the
    typical patient (the parameters + rate rules), but the curated NLME object travels
    with the model and survives a round-trip. Returns the inner predicate lines (to be
    embedded inside the model's single ``rdf:RDF`` block) — or ``""`` when the model
    publishes no variability (the never-synthesize rule: absence is a true gap).
    """
    from ._variability import omega_correlations, omega_diagonal, residual_spec

    if not model.has_published_variability:
        return ""
    pad = f"{indent}    "
    lines = [
        f"{pad}<hypnos:variabilityStatus>{model.variability_status}</hypnos:variabilityStatus>",
        f"{pad}<hypnos:bandTier>{model.band_tier}</hypnos:bandTier>",
        f"{pad}<hypnos:betweenSubjectVariability>",
        f"{pad}  <rdf:Bag>",
    ]
    for sym, om2, cv in omega_diagonal(model):
        lines.append(f'{pad}    <rdf:li hypnos:parameter="{sym}" '
                     f'hypnos:omega2="{om2:.10g}" hypnos:cvPercent="{cv:.4g}"/>')
    lines += [f"{pad}  </rdf:Bag>", f"{pad}</hypnos:betweenSubjectVariability>"]
    for a, b, r in omega_correlations(model):
        lines.append(f'{pad}<hypnos:omegaCorrelation hypnos:between="{a} {b}" '
                     f'hypnos:correlation="{r:.6g}"/>')
    spec = residual_spec(model)
    if spec is not None:
        attrs = f'hypnos:model="{spec.model}"'
        if spec.log_sd is not None:
            attrs += f' hypnos:logSd="{spec.log_sd:.10g}"'
        if spec.prop_var is not None:
            attrs += f' hypnos:proportionalVariance="{spec.prop_var:.10g}"'
        if spec.add_sd is not None:
            attrs += f' hypnos:additiveSd="{spec.add_sd:.10g}"'
        lines.append(f"{pad}<hypnos:residualError {attrs}/>")
    return "\n".join(lines)


def estimation_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates describing the curated estimation uncertainty (v0.3 §8).

    The per-θ SE/RSE and any published estimate correlations ride as metadata beside the
    Ω/Σ of :func:`variability_rdf` — in their OWN predicates, because an estimation RSE is
    not a between-subject CV (v0.3 Trap 1). Each ``rdf:li`` carries the curated typical
    value the SE is on, which is *not* the patient-instantiated ``<parameter>`` value.
    Returns ``""`` when the model publishes no estimation uncertainty (never-synthesize).
    """
    from ._estimation import (
        covariance_step_note, estimate_correlations, independence_caveat, theta_estimates,
    )

    ests = theta_estimates(model)
    if not ests:
        return ""
    pad = f"{indent}    "
    lines = [
        f"{pad}<hypnos:uncertaintyStatus>{model.uncertainty_status}</hypnos:uncertaintyStatus>",
        f"{pad}<hypnos:estimationBandTier>{model.estimation_band_tier}</hypnos:estimationBandTier>",
        f"{pad}<hypnos:parameterStandardError>",
        f"{pad}  <rdf:Bag>",
    ]
    for e in ests:
        rse = f' hypnos:rsePercent="{e.rse_percent:.6g}"' if e.rse_percent is not None else ""
        ci = (f' hypnos:ci95Low="{e.ci95[0]:.10g}" hypnos:ci95High="{e.ci95[1]:.10g}"'
              if e.ci95 else "")
        lines.append(
            f'{pad}    <rdf:li hypnos:parameter="{e.symbol}" hypnos:estimate="{e.central:.10g}" '
            f'hypnos:se="{e.se:.10g}" hypnos:scale="{e.scale}"{rse}{ci} '
            f'hypnos:method="{e.method}" hypnos:tier="{e.tier}"/>')
    lines += [f"{pad}  </rdf:Bag>", f"{pad}</hypnos:parameterStandardError>"]
    for a, b, r in estimate_correlations(model):
        lines.append(f'{pad}<hypnos:estimateCorrelation hypnos:between="{a} {b}" '
                     f'hypnos:correlation="{r:.6g}"/>')
    for note in (covariance_step_note(model), independence_caveat(model)):
        if note:
            lines.append(f"{pad}<hypnos:estimationNote>{note}</hypnos:estimationNote>")
    return "\n".join(lines)


def la_rdf(model: Model, drug: Optional[Dict[str, Any]] = None, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the local-anesthetic blocks (v0.6 §9).

    Carries the site-absorption rank order, the toxicity-threshold *ranges*, and (v0.6
    LA2) the drug's ``cardiotoxicity_class`` as metadata so they travel with an
    SBML/PharmML model whose core is deterministic. Threshold ranges export **as
    ranges** — there is no format projection that collapses them to a single value
    (v0.6 §9). ``drug`` is the resolved drug record (for the cardiotoxicity class);
    pass ``None`` to omit it. Returns ``""`` for a non-LA model (and emits nothing for
    blocks the model does not carry)."""
    if not model.is_safety_critical:
        return ""
    pad = f"{indent}    "
    lines = [f"{pad}<hypnos:safetyCritical>true</hypnos:safetyCritical>"]
    absn = model.absorption or {}
    rates = absn.get("site_rates", [])
    if rates:
        lines.append(f"{pad}<hypnos:siteAbsorption>")
        lines.append(f"{pad}  <rdf:Seq>")
        for s in sorted(rates, key=lambda r: r.get("rank", 99)):
            ka = s.get("ka")
            ka_attr = f' hypnos:ka="{ka:.10g}"' if ka is not None else ""
            lines.append(f'{pad}    <rdf:li hypnos:site="{s.get("site")}" '
                         f'hypnos:rank="{s.get("rank")}"{ka_attr}/>')
        lines += [f"{pad}  </rdf:Seq>", f"{pad}</hypnos:siteAbsorption>"]
    for th in model.toxicity_thresholds:
        lines.append(
            f'{pad}<hypnos:toxicityThresholdRange hypnos:endpoint="{th.endpoint}" '
            f'hypnos:basis="{th.basis}" hypnos:low="{th.low:.10g}" hypnos:high="{th.high:.10g}" '
            f'hypnos:units="{th.units}" hypnos:tier="{th.tier}"/>')
    cc = (drug or {}).get("cardiotoxicity_class") if drug else None
    if cc:
        lines.append(
            f'{pad}<hypnos:cardiotoxicityClass hypnos:rank="{cc.get("rank")}" '
            f'hypnos:stereochemistry="{cc.get("stereochemistry")}" '
            f'hypnos:cnsToCvsMargin="{cc.get("cns_to_cvs_margin")}"/>')
    ffm = ((drug or {}).get("protein_binding") or {}).get("free_fraction_model") if drug else None
    if ffm:
        lines.append(
            f'{pad}<hypnos:freeFractionModel hypnos:type="{ffm.get("type")}" '
            f'hypnos:bindingCapacityUgMl="{ffm.get("binding_capacity_ug_ml")}" '
            f'hypnos:tier="{ffm.get("tier")}"/>')
    return "\n".join(lines)


def protein_binding_rdf(drug: Optional[Dict[str, Any]], indent: str = "  ") -> str:
    """``hypnos:proteinBinding`` for the drug-level binding block (v0.6 §9 / v0.5 §B3).

    The bound fraction and the two flags that say when a TOTAL concentration under-states
    free-drug exposure — ``saturable`` (the LA failure mode) and ``bindingSensitive`` (the
    hypoalbuminemia one). Carried as stated; never folded into a parameter. ``""`` when the
    drug curates no binding."""
    pb = (drug or {}).get("protein_binding")
    if not pb or pb.get("fraction_bound") is None:
        return ""
    attrs = f'hypnos:fractionBound="{pb["fraction_bound"]:.10g}"'
    if pb.get("binding_protein"):
        attrs += f' hypnos:bindingProtein="{pb["binding_protein"]}"'
    for key, name in (("saturable", "saturable"), ("binding_sensitive", "bindingSensitive")):
        if pb.get(key) is not None:
            attrs += f' hypnos:{name}="{str(bool(pb[key])).lower()}"'
    return f"{indent}    <hypnos:proteinBinding {attrs}/>"


def developmental_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the developmental extrapolation block (v0.8 §7).

    The allometric size scaling and the PMA maturation function ride as metadata so they
    travel with an SBML/PharmML model whose core is deterministic, carrying the Tier-D flag,
    the opt-in flag, and (for ``allometry_only``) the over-dose caveat. The maturation driver
    is emitted as PMA with a comment forbidding chronological-age substitution (Trap 1).
    Returns ``""`` when the model carries no developmental block."""
    dev = model.developmental_model
    if dev is None:
        return ""
    pad = f"{indent}    "
    lines = [
        f'{pad}<hypnos:developmentalModel hypnos:extrapolationBasis="{dev.extrapolation_basis}" '
        f'hypnos:evidenceTier="{dev.evidence_tier}" hypnos:appliedByDefault="{str(dev.applied_by_default).lower()}">',
    ]
    if dev.size is not None:
        s = dev.size
        desc = f' hypnos:sizeDescriptor="{s.size_descriptor}"' if s.size_descriptor else ""
        lines.append(
            f'{pad}  <hypnos:allometricSize hypnos:clExponent="{s.cl_exponent:g}" '
            f'hypnos:vExponent="{s.v_exponent:g}" hypnos:referenceWeightKg="{s.reference_weight_kg:g}" '
            f'hypnos:exponentBasis="{s.exponent_basis}"{desc} hypnos:tier="{s.tier}"/>')
    if dev.maturation is not None:
        m = dev.maturation
        lines.append(
            f'{pad}  <!-- maturation driver is PMA (post-menstrual age); a downstream consumer '
            f'MUST NOT substitute chronological age (v0.8 Trap 1) -->')
        lines.append(
            f'{pad}  <hypnos:maturationFunction hypnos:function="{m.function}" '
            f'hypnos:tm50Weeks="{m.tm50_weeks:g}" hypnos:hill="{m.hill:g}" hypnos:driver="{m.driver}" '
            f'hypnos:affectedParameter="{m.affected_parameter}" hypnos:tier="{m.tier}"/>')
    else:
        lines.append(f'{pad}  <hypnos:maturationCaveat>allometry_only — maturation un-modeled; '
                     f'neonatal clearance OVER-stated, exposure UNDER-predicted</hypnos:maturationCaveat>')
    lines.append(f'{pad}  <hypnos:developmentalExtrapolation>true</hypnos:developmentalExtrapolation>')
    lines.append(f'{pad}</hypnos:developmentalModel>')
    return "\n".join(lines)


def pharmacogenomics_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the pharmacogenomic blocks (v0.9 §7).

    Kinetic modifiers ride as ``hypnos:pharmacogenomicModifier`` metadata (genotype is a
    covariate the core formats describe only weakly); safety flags ride as
    ``hypnos:safetyCritical`` + ``hypnos:pharmacogenomicSafetyFlag``, **explicitly not** as
    any parameter change — the RDF carrier keeps the avoidance fact lossless and portable
    without ever encoding it as a number a deterministic consumer might apply. Returns ``""``
    when the model carries no pharmacogenomic blocks."""
    if not model.has_pharmacogenomics:
        return ""
    pad = f"{indent}    "
    lines: List[str] = []
    for mod in model.pharmacogenomic_modifiers:
        sf = mod.scale_factor
        sf_attr = f' hypnos:scaleFactor="{sf:.10g}"' if sf is not None else ""
        lines.append(
            f'{pad}<hypnos:pharmacogenomicModifier hypnos:gene="{mod.gene}" '
            f'hypnos:phenotype="{mod.phenotype_value}" hypnos:affectedParameter="{mod.affected_parameter}"'
            f'{sf_attr} hypnos:substrateScope="{" ".join(mod.substrate_scope)}" '
            f'hypnos:evidenceTier="{mod.evidence_tier}" hypnos:appliedByDefault="false"/>')
    for flag in model.pharmacogenomic_safety_flags:
        lines.append(
            f'{pad}<hypnos:pharmacogenomicSafetyFlag hypnos:gene="{flag.gene}" '
            f'hypnos:kind="{flag.kind}" hypnos:triggers="{" ".join(flag.trigger_agents)}" '
            f'hypnos:affectsKinetics="false" hypnos:evidenceTier="{flag.evidence_tier}"/>')
    return "\n".join(lines)


def _reversal_facts(rm: Dict[str, Any]) -> List[tuple]:
    """The mechanism-appropriate numbers of a ``reversal_mechanism`` as ``[(name, value)]``.

    Only what the record curates: a null binding constant is omitted, never zero-filled."""
    facts: List[tuple] = []
    for block, keys in (("stoichiometry", ("ratio", "basis")),
                        ("binding", ("k_on", "k_off", "kd", "approximation")),
                        ("antagonism", ("ki", "ki_units", "ce50_shift_model")),
                        ("indirect", ("ceiling", "depth_floor"))):
        for k in keys:
            v = (rm.get(block) or {}).get(k)
            if v is not None:
                facts.append((f"{block}.{k}", v))
    return facts


def _camel(dotted: str) -> str:
    """``binding.k_off`` -> ``bindingKOff`` (an XML attribute local name)."""
    parts = dotted.replace(".", "_").split("_")
    return parts[0] + "".join(p.capitalize() for p in parts[1:])


def reversal_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:`` RDF for the reversal-agent coupling (v0.5 §D).

    SBML/NONMEM core cannot express a cross-model binding reaction, so the mechanism,
    its stoichiometry/binding/antagonism/ceiling numbers, and the **link to the agonist
    model(s)** it acts on ride as metadata — antagonism is an interaction, and an export
    of the antagonist alone must still say what it reverses. Returns ``""`` for a model
    with no ``reversal_mechanism``."""
    rm = model.reversal_mechanism
    if not rm:
        return ""
    pad = f"{indent}    "
    attrs = "".join(f" hypnos:{_camel(k)}={quoteattr(str(v))}" for k, v in _reversal_facts(rm))
    lines = [f'{pad}<hypnos:reversalMechanism hypnos:type="{rm.get("type")}" '
             f'hypnos:tier="{rm.get("tier", "D")}"{attrs}>',
             f"{pad}  <hypnos:reversalTargets>", f"{pad}    <rdf:Bag>"]
    lines += [f'{pad}      <rdf:li hypnos:modelId="{t}"/>' for t in model.reversal_targets]
    lines += [f"{pad}    </rdf:Bag>", f"{pad}  </hypnos:reversalTargets>",
              f"{pad}</hypnos:reversalMechanism>"]
    return "\n".join(lines)


def _citation_uri(cid: Optional[str], ds) -> Optional[str]:
    cit = ds.citation(cid) if (ds is not None and cid) else None
    if cit and cit.get("doi"):
        return f"https://doi.org/{cit['doi']}"
    if cit and cit.get("pmid"):
        return f"https://identifiers.org/pubmed:{cit['pmid']}"
    return None


def _condition_text(dz) -> str:
    return f"{dz.dimension} < {dz.below:g}" if dz.below is not None else f"{dz.dimension} = {dz.value}"


def disease_modifier_rdf(model: Model, ds=None, indent: str = "  ") -> str:
    """``hypnos:diseaseStateModifier`` annotations (v0.5 §D).

    Each carries its citation, tier, and the opt-in flag, so a downstream consumer cannot
    apply one without seeing that it is an unvalidated, opt-in, Tier-D extrapolation. The
    exported parameters are **never** pre-adjusted: the modifier travels beside them.
    Returns ``""`` when the model carries none."""
    pad = f"{indent}    "
    lines: List[str] = []
    for dz in model.disease_state_modifiers:
        sf = dz.scale_factor
        sf_attr = f' hypnos:scaleFactor="{sf:.10g}"' if sf is not None else ""
        uri = _citation_uri(dz.primary_citation, ds)
        cite = f" hypnos:citation={quoteattr(uri)}" if uri else ""
        lines.append(
            f"{pad}<hypnos:diseaseStateModifier hypnos:condition={quoteattr(_condition_text(dz))} "
            f'hypnos:affectedParameter="{dz.affected_parameter}"{sf_attr} '
            f'hypnos:evidenceTier="{dz.evidence_tier}" '
            f'hypnos:appliedByDefault="{str(dz.applied_by_default).lower()}"{cite}/>')
    return "\n".join(lines)


def _validation_entries(model: Model) -> List[Dict[str, Any]]:
    """Hypnos-COMPUTED metric sets flattened for annotation: one dict per entry with the
    cohort, the stratification axis, and ``{metric name: value}`` (v0.4 §4.1)."""
    out = []
    for e in model.external_validation:
        cohort = e.get("cohort") or {}
        out.append({
            "dataset": e.get("dataset"), "mode": e.get("mode"), "target": e.get("target"),
            "n_subjects": cohort.get("n_subjects"), "in_envelope": cohort.get("in_envelope"),
            "metrics": [(m.get("name"), m.get("value"), m.get("units") or "")
                        for m in e.get("metrics", []) if m.get("value") is not None],
            "checksum": (e.get("provenance") or {}).get("manifest_checksum"),
        })
    return out


_NOT_CLEARANCE = ("Hypnos-COMPUTED and reproducible — distinct from publisher-reported "
                  "predictive_performance; good external performance is NOT a clinical clearance")


def external_validation_rdf(model: Model, indent: str = "  ") -> str:
    """``hypnos:externalValidation`` RDF (v0.4 §8).

    These formats describe the model, not its post-hoc validation, so the computed
    Varvel metrics + their cohort ride as annotation with ``hypnos:validationStatus``.
    Each entry keeps ``inEnvelope`` — a metric without its stratification axis would
    read as a blanket claim. Returns ``""`` when the model carries no computed metrics."""
    entries = _validation_entries(model)
    if not entries:
        return ""
    pad = f"{indent}    "
    lines = [f"{pad}<hypnos:validationStatus>{model.validation_status}</hypnos:validationStatus>",
             f"{pad}<!-- {_NOT_CLEARANCE} -->"]
    for e in entries:
        attrs = (f'hypnos:dataset={quoteattr(str(e["dataset"]))} hypnos:mode="{e["mode"]}" '
                 f'hypnos:target="{e["target"]}"')
        if e["n_subjects"] is not None:
            attrs += f' hypnos:nSubjects="{e["n_subjects"]}"'
        if e["in_envelope"] is not None:
            attrs += f' hypnos:inEnvelope="{str(e["in_envelope"]).lower()}"'
        if e["checksum"]:
            attrs += f' hypnos:manifestChecksum={quoteattr(str(e["checksum"]))}'
        lines.append(f"{pad}<hypnos:externalValidation {attrs}>")
        lines += [f'{pad}  <hypnos:metric hypnos:name="{n}" hypnos:value="{v:.10g}" '
                  f'hypnos:units={quoteattr(u)}/>' for n, v, u in e["metrics"]]
        lines.append(f"{pad}</hypnos:externalValidation>")
    return "\n".join(lines)


def breadth_comments(model: Model, ds=None) -> List[str]:
    """The v0.5 reversal + disease-modifier blocks, the v0.9 pharmacogenomic records, and the
    v0.4 computed-validation block as plain comment lines (no prefix).

    For the text formats (NONMEM / rxode2 / Pumas), whose core can express neither a
    cross-model binding reaction nor an opt-in adjustment: the same facts the RDF carries,
    as comments, so the annotation is never format-dependent (v0.5 §D)."""
    lines: List[str] = []
    rm = model.reversal_mechanism
    if rm:
        lines.append(f"hypnos:reversalMechanism = {rm.get('type')} (tier {rm.get('tier', 'D')}) — "
                     "an INTERACTION with the target model(s), not a solo trajectory")
        lines += [f"  hypnos:reversalTarget = {t}" for t in model.reversal_targets]
        lines += [f"  {k} = {v}" for k, v in _reversal_facts(rm)]
    for dz in model.disease_state_modifiers:
        sf = dz.scale_factor
        uri = _citation_uri(dz.primary_citation, ds)
        lines.append(
            f"hypnos:diseaseStateModifier = {_condition_text(dz)}: {dz.affected_parameter}"
            + (f" x {sf:g}" if sf is not None else "")
            + f" (evidenceTier {dz.evidence_tier}, appliedByDefault "
              f"{str(dz.applied_by_default).lower()}" + (f", {uri}" if uri else "") + ")")
        lines.append("  OPT-IN, UNVALIDATED EXTRAPOLATION — NOT applied to the parameters in this file")
    # v0.9 §7: pharmacogenomic records. A modifier is opt-in Tier-D; a safety flag is an
    # avoidance/awareness FACT and never a parameter change.
    for mod in model.pharmacogenomic_modifiers:
        sf = mod.scale_factor
        lines.append(
            f"hypnos:pharmacogenomicModifier = {mod.gene} ({mod.phenotype_value}): "
            f"{mod.affected_parameter}" + (f" x {sf:g}" if sf is not None else "")
            + f"; substrateScope {' '.join(mod.substrate_scope)} (evidenceTier {mod.evidence_tier}, "
              "appliedByDefault false)")
        lines.append("  OPT-IN, UNVALIDATED EXTRAPOLATION — NOT applied to the parameters in this file")
    for flag in model.pharmacogenomic_safety_flags:
        lines.append(
            f"hypnos:pharmacogenomicSafetyFlag = {flag.gene} ({flag.kind}): {flag.action} "
            f"(evidenceTier {flag.evidence_tier}, affectsKinetics false)")
        lines.append("  hypnos:safetyCritical = true — a contraindication/awareness fact, NOT a "
                     "kinetic number; nothing in this file is adjusted for it")
    # v0.4 §8: the Hypnos-computed external-validation metrics + their cohort.
    entries = _validation_entries(model)
    if entries:
        lines.append(f"hypnos:validationStatus = {model.validation_status}")
        for e in entries:
            where = {True: "in-envelope", False: "OUT-of-envelope", None: "unstratified"}[e["in_envelope"]]
            n = f", n={e['n_subjects']}" if e["n_subjects"] is not None else ""
            mets = ", ".join(f"{nm} {v:g}{u}" for nm, v, u in e["metrics"])
            lines.append(f"hypnos:externalValidation = {e['dataset']} ({e['mode']} -> {e['target']}, "
                         f"{where}{n}): {mets}")
        lines.append(f"  {_NOT_CLEARANCE}")
    return lines


def rdf_annotation_xml(model: Model, ds=None, tier: Optional[str] = None, indent: str = "  ",
                       extra_predicates: str = "") -> str:
    """A small RDF/MIRIAM block embeddable in SBML/PharmML <annotation> elements.

    ``extra_predicates`` (e.g. :func:`variability_rdf`) is injected inside the single
    ``rdf:RDF`` element — SBML allows only one ``<annotation>`` per element, so the
    random-effects layer rides alongside the provenance rather than in a second block.
    """
    t = tier or model.tier
    uris = citation_uris(model, ds)
    bag = "\n".join(
        f'{indent}      <rdf:li rdf:resource="{u}"/>' for u in uris
    ) or f'{indent}      <rdf:li rdf:resource="urn:hypnos:no-citation"/>'
    extra = (extra_predicates + "\n") if extra_predicates else ""
    return (
        f'{indent}<annotation>\n'
        f'{indent}  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n'
        f'{indent}           xmlns:bqmodel="http://biomodels.net/model-qualifiers/"\n'
        f'{indent}           xmlns:hypnos="{HYPNOS_NS}">\n'
        f'{indent}    <hypnos:clinicalUse>{CLINICAL_USE}</hypnos:clinicalUse>\n'
        f'{indent}    <hypnos:confidenceTier>{t}</hypnos:confidenceTier>\n'
        f'{indent}    <hypnos:datasetVersion>{__version__}</hypnos:datasetVersion>\n'
        f'{indent}    <bqmodel:isDerivedFrom>\n'
        f'{indent}      <rdf:Bag>\n{bag}\n{indent}      </rdf:Bag>\n'
        f'{indent}    </bqmodel:isDerivedFrom>\n'
        f'{extra}'
        f'{indent}  </rdf:RDF>\n'
        f'{indent}</annotation>'
    )
