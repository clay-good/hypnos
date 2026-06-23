"""Tests for the v0.9 pharmacogenomic-envelope subsystem.

The load-bearing discipline: a kinetic modifier scales a prediction; a safety flag forbids
a trigger and carries NO number. The schema separates the two, the validator enforces the
substrate guardrail, and nothing is ever inferred — a modifier/flag fires only on a
DECLARED genotype.
"""
import pytest

import hypnos
from hypnos.models import Model
from hypnos.pharmacogenomics import genotype_triggers, pgx_overlay
from hypnos.simulate import _apply_pgx_modifiers, evaluate_safety
from hypnos.validate import _check_pharmacogenomics
from hypnos.export import annotate

SUX = "nmb_agents.succinylcholine.roy_2004"
SEVO = "volatiles.sevoflurane.mac"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --- never-infer ------------------------------------------------------------ #
def test_no_genotype_no_effect(ds):
    pg = pgx_overlay(ds, genotype={})
    assert pg.safety_flags == [] and pg.modifiers == []
    # an undeclared dimension never triggers
    assert not genotype_triggers({}, "mh_susceptibility", True)
    assert not genotype_triggers({"bche_phenotype": "normal"}, "bche_phenotype", "homozygous_atypical")


def test_declared_atypical_triggers_but_normal_does_not():
    assert genotype_triggers({"bche_phenotype": "homozygous_atypical"},
                             "bche_phenotype", "homozygous_atypical")
    assert genotype_triggers({"bche_phenotype": "heterozygous_atypical"},
                             "bche_phenotype", "homozygous_atypical")   # any atypical triggers
    assert not genotype_triggers({"bche_phenotype": "normal"}, "bche_phenotype", "homozygous_atypical")


# --- the overlay: avoidance first, kinetics second, honest majority --------- #
def test_mh_overlay_avoidance_flags(ds):
    pg = pgx_overlay(ds, genotype={"mh_susceptibility": True})
    flagged_drugs = {f.drug for f in pg.avoidance_flags}
    assert {"sevoflurane", "isoflurane", "desflurane", "succinylcholine"} <= flagged_drugs
    assert all(f.record.gene == "RYR1" for f in pg.avoidance_flags)
    # a susceptibility carries no kinetic effect (Trap 2)
    assert all(f.record.affects_kinetics is False for f in pg.avoidance_flags)


def test_bche_overlay_modifier_and_awareness_and_guardrails(ds):
    pg = pgx_overlay(ds, genotype={"bche_phenotype": "homozygous_atypical"})
    assert any(f.drug == "succinylcholine" and f.gene == "BCHE" for f in pg.modifiers)
    assert any(f.drug == "succinylcholine" and f.record.kind == "awareness" for f in pg.awareness_flags)
    # substrate guardrail: a genetic effect never leaks to a non-substrate
    guarded = {g["drug"] for g in pg.substrate_guardrails}
    assert {"remifentanil", "lidocaine", "bupivacaine", "ropivacaine", "rocuronium"} <= guarded
    # the honest majority: drugs with no actionable evidence are named, not invented
    assert "propofol" in pg.no_evidence and "dexmedetomidine" in pg.no_evidence


def test_modifier_is_substrate_scoped(ds):
    mod = ds[SUX].pharmacogenomic_modifiers[0]
    assert "succinylcholine" in mod.substrate_scope
    assert "remifentanil" not in mod.substrate_scope     # NOT a BCHE substrate
    assert mod.evidence_tier == "D" and mod.applied_by_default is False


# --- safety-flag surfacing in simulate (a warning, never a number) ---------- #
def test_evaluate_safety_surfaces_declared_mh_flag(ds):
    floor, warns, _ = evaluate_safety(ds[SEVO], {"mh_susceptibility": True})
    assert any("PGx SAFETY [AVOID, RYR1]" in w for w in warns)
    # the flag does not change the tier floor (no kinetic effect)
    assert floor == ds[SEVO].tier


def test_evaluate_safety_silent_without_declaration(ds):
    _, warns, _ = evaluate_safety(ds[SEVO], {})
    assert not any("PGx SAFETY" in w for w in warns)


# --- kinetic modifier application (forward-only, substrate-scoped) ---------- #
def test_apply_pgx_modifier_scales_clearance(ds):
    from hypnos.reference import MicroParams
    params = MicroParams.from_volumes_clearances(V1=10.0, Cl1=2.0, V2=20.0, Cl2=1.0)
    out, warns = _apply_pgx_modifiers(ds[SUX], params, {"bche_phenotype": "homozygous_atypical"})
    # Cl scaled by the curated factor (0.2) -> reduced hydrolysis -> prolonged effect
    assert out.as_volumes_clearances()["Cl1"] == pytest.approx(2.0 * 0.2)
    assert any("PGx MODIFIER" in w for w in warns)


def test_apply_pgx_modifier_respects_substrate_guardrail(ds):
    """A BCHE modifier must not fire on a non-substrate drug even if its phenotype is declared."""
    from hypnos.reference import MicroParams
    # remifentanil is not in the modifier's substrate_scope; build a fake model that
    # (wrongly) carries the succinylcholine modifier but is a remifentanil drug record
    raw = dict(ds[SUX].raw)
    raw = {**raw, "drug": {"name": "remifentanil"}}
    fake = Model(raw)
    params = MicroParams.from_volumes_clearances(V1=10.0, Cl1=2.0)
    out, warns = _apply_pgx_modifiers(fake, params, {"bche_phenotype": "homozygous_atypical"})
    assert out.as_volumes_clearances()["Cl1"] == pytest.approx(2.0)   # unchanged (guardrail)


# --- validation ------------------------------------------------------------- #
def test_validate_accepts_dataset(ds):
    assert hypnos.validate_dataset(ds) == []


def test_validate_rejects_mis_scoped_modifier():
    bad = Model({"id": "opioids.remifentanil.x", "drug": {"name": "remifentanil"}, "purpose": "pk",
                 "tier": "D", "primary_citation": "minto-1997-remifentanil",
                 "extraction": {"review_status": "unverified"},
                 "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                 "parameters": [],
                 "pharmacogenomic_modifiers": [{
                     "gene": "BCHE", "phenotype": {"dimension": "bche_phenotype", "value": "homozygous_atypical"},
                     "affected_parameter": "Cl", "substrate_scope": ["succinylcholine"],
                     "adjustment": {"type": "scale_factor", "value": 0.2}, "evidence_tier": "D",
                     "applied_by_default": False, "primary_citation": "lockridge-2015-bche"}]})
    problems = _check_pharmacogenomics(bad, {"lockridge-2015-bche"})
    assert any("substrate_scope" in p for p in problems)


def test_validate_rejects_default_applied_modifier():
    bad = Model({"id": "nmb_agents.succinylcholine.x", "drug": {"name": "succinylcholine"},
                 "purpose": "pk", "tier": "D", "primary_citation": "lockridge-2015-bche",
                 "extraction": {"review_status": "unverified"},
                 "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                 "parameters": [],
                 "pharmacogenomic_modifiers": [{
                     "gene": "BCHE", "phenotype": {"dimension": "bche_phenotype", "value": "homozygous_atypical"},
                     "affected_parameter": "Cl", "substrate_scope": ["succinylcholine"],
                     "adjustment": {"type": "scale_factor", "value": 0.2}, "evidence_tier": "D",
                     "applied_by_default": True, "primary_citation": "lockridge-2015-bche"}]})
    assert any("applied_by_default" in p for p in _check_pharmacogenomics(bad, {"lockridge-2015-bche"}))


# --- export ----------------------------------------------------------------- #
def test_pgx_rdf_carries_flag_as_non_kinetic(ds):
    rdf = annotate.pharmacogenomics_rdf(ds[SUX])
    assert 'hypnos:affectsKinetics="false"' in rdf           # a flag is never a number
    assert "hypnos:pharmacogenomicModifier" in rdf
    assert "hypnos:pharmacogenomicSafetyFlag" in rdf


def test_safety_critical_propagates_to_export(ds):
    # a model carrying a pgx safety flag is safety-critical (twice-warned in exports)
    assert ds[SEVO].is_safety_critical
    assert ds[SUX].is_safety_critical
