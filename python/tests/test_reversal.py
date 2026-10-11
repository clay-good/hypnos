"""v0.5 Part A — the reversal subsystem.

Antagonism is an interaction, not a solo drug (§A1). Three mechanisms, three kernels:
sugammadex encapsulation (TOF recovers), naloxone competitive antagonism (the
renarcotization relapse is rendered forward), neostigmine indirect inhibition (the hard
deep-block ceiling). No reversal dose is ever computed (§A6).
"""
import numpy as np
import pytest

import hypnos
from hypnos.export.registry import rocuronium_wierda_1991
from hypnos.reversal import competitive_shift, encapsulation, indirect_inhibition

SUGAMMADEX = "reversal.sugammadex.encapsulation"
NALOXONE = "reversal.naloxone.competitive"
NEOSTIGMINE = "reversal.neostigmine.indirect"
ROC = "nmb_agents.rocuronium.wierda_1991"
TOF = "pd_effect.rocuronium.tof_sigmoid"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --- kernels --------------------------------------------------------------- #
def test_encapsulation_forms_complex_and_drops_free_drug():
    roc = rocuronium_wierda_1991({"weight": 70, "age": 40})
    t = np.linspace(0, 40, 400)
    no = encapsulation(roc, 42.0, 0.0, t, sugammadex_time_min=10.0)
    rev = encapsulation(roc, 42.0, 280.0, t, sugammadex_time_min=10.0)
    i = int(np.argmin(np.abs(t - 25)))
    assert rev.free_ce[i] < no.free_ce[i]            # reversal lowers free drug
    assert rev.complex_umol[i] > 0                   # an inert complex forms
    assert no.complex_umol[i] == pytest.approx(0.0)  # no sugammadex -> no complex


def test_competitive_shift_multiplier():
    # apparent Ce50 multiplier = 1 + C/Ki
    m = competitive_shift(np.array([0.0, 1.0, 3.0]), ki=1.0)
    np.testing.assert_allclose(m, [1.0, 2.0, 4.0])


def test_indirect_inhibition_ceiling_and_floor():
    # below the depth floor (deep block) -> unreversed; above -> partial reversal up to the ceiling
    block = np.array([1.0, 0.95, 0.5, 0.1])
    res = indirect_inhibition(block, ceiling=0.5, depth_floor=0.9)
    assert res[0] == pytest.approx(1.0)              # deep block: unchanged (ceiling)
    assert res[1] == pytest.approx(0.95)             # still below floor: unchanged
    assert res[2] == pytest.approx(0.25)             # 0.5 * (1 - 0.5)
    assert res[3] == pytest.approx(0.05)


# --- simulate_reversal end-to-end ----------------------------------------- #
def test_encapsulation_recovers_tof(ds):
    t = np.linspace(0, 40, 400)
    r = hypnos.simulate_reversal(ds, SUGAMMADEX, patient=dict(age=40, weight=70),
                                 reversal_dose="4 mg/kg", agonist_dose="0.6 mg/kg", t=t,
                                 reversal_time_min=10, pd_model=TOF)
    i30 = int(np.argmin(np.abs(t - 30)))
    # with reversal TOF (T1%) recovers far more than without
    assert r.effect[i30] > r.effect_no_reversal[i30] + 50
    assert any("never computes a reversal dose" in w for w in r.warnings)


def test_competitive_renarcotization_relapses(ds):
    # over 4 h naloxone (t1/2 ~65 min) clears faster than a long-acting opioid -> relapse
    t = np.linspace(0, 240, 600)
    r = hypnos.simulate_reversal(ds, NALOXONE, patient=dict(age=40, weight=70),
                                 reversal_dose="0.4 mg", agonist_dose="2 mcg/kg", t=t,
                                 reversal_time_min=10)
    m = r.ce50_multiplier
    assert m.max() > 1.5                              # naloxone shifts the apparent Ce50 up
    assert m[-1] < m.max()                            # then relapses toward 1 as it clears
    assert r.renarcotization is True
    assert any("RENARCOTIZATION" in w for w in r.warnings)


def test_enzyme_inhibition_deep_block_ceiling(ds):
    t = np.linspace(0, 40, 200)
    r = hypnos.simulate_reversal(ds, NEOSTIGMINE, patient=dict(age=40, weight=70),
                                 reversal_dose="2.5 mg", agonist_dose="0.6 mg/kg", t=t)
    assert r.effect[0] == pytest.approx(1.0)          # deep block unreversed (ceiling)
    assert r.effect[-1] < 0.5                         # shallow block reversed
    assert any("CEILING" in w for w in r.warnings)


# --- model + validate ------------------------------------------------------ #
def test_reversal_records_resolve_targets(ds):
    for rid in (SUGAMMADEX, NALOXONE, NEOSTIGMINE):
        m = ds[rid]
        assert m.is_reversal_agent
        for tgt in m.reversal_targets:
            assert tgt in ds                           # targets must resolve


def test_validate_requires_renarcotization_caveat(ds):
    from hypnos.models import Model
    from hypnos.validate import _check_reversal
    # a competitive_receptor agent with NO renarcotization failure mode is rejected (§A5).
    # the target (fentanyl) resolves in the real ds, so only the missing caveat should flag.
    bad = Model({"id": "reversal.x.y", "drug": {"name": "naloxone"}, "purpose": "reversal",
                 "tier": "C", "primary_citation": "clarke-2005-naloxone",
                 "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
                 "parameters": [], "extraction": {"review_status": "unverified"},
                 "reversal_mechanism": {"type": "competitive_receptor",
                                        "targets": ["opioids.fentanyl.shafer_1990"],
                                        "antagonism": {"ki": 1.0},
                                        "primary_citation": "clarke-2005-naloxone"}})
    probs = _check_reversal(bad, ds, {"clarke-2005-naloxone"})
    assert any("renarcotization" in p for p in probs)


def test_full_dataset_validates_clean():
    assert hypnos.validate_dataset() == []


# --------------------------------------------------------------------------- #
# v0.5 §D — the reversal coupling travels in every export; §C — the verify group
# --------------------------------------------------------------------------- #
import json  # noqa: E402
import xml.dom.minidom as minidom  # noqa: E402

from hypnos.export import export_model  # noqa: E402
from hypnos.verification import checklist_markdown, model_verification  # noqa: E402

SUGAMMADEX = "reversal.sugammadex.encapsulation"
NALOXONE = "reversal.naloxone.competitive"
NEOSTIGMINE = "reversal.neostigmine.indirect"
_AGENTS = [SUGAMMADEX, NALOXONE, NEOSTIGMINE]
_MODEL_FORMATS = ["nonmem", "pharmml", "sbml", "tci_json", "rxode2", "pumas"]


@pytest.mark.parametrize("mid", _AGENTS)
@pytest.mark.parametrize("fmt", _MODEL_FORMATS)
def test_every_export_links_the_antagonist_to_its_target(fmt, mid):
    # an export of the antagonist alone must still say WHAT it reverses and HOW
    ds = hypnos.load()
    m = ds[mid]
    text = export_model(fmt, m, ds)[1]
    assert m.reversal_mechanism["type"] in text
    for target in m.reversal_targets:
        assert target in text


def test_tci_json_passes_reversal_mechanism_verbatim():
    ds = hypnos.load()
    doc = json.loads(export_model("tci_json", ds[SUGAMMADEX], ds)[1])
    assert doc["reversal_mechanism"] == ds[SUGAMMADEX].raw["reversal_mechanism"]
    # a non-reversal model carries no such key
    roc = json.loads(export_model("tci_json", ds["nmb_agents.rocuronium.wierda_1991"], ds)[1])
    assert "reversal_mechanism" not in roc


def test_sbml_reversal_rdf_is_well_formed_and_carries_only_curated_numbers():
    ds = hypnos.load()
    dom = minidom.parseString(export_model("sbml", ds[SUGAMMADEX], ds)[1])
    (rm,) = dom.getElementsByTagName("hypnos:reversalMechanism")
    assert rm.getAttribute("hypnos:type") == "encapsulation"
    assert rm.getAttribute("hypnos:stoichiometryBasis") == "molar"      # molar, never mass
    assert float(rm.getAttribute("hypnos:stoichiometryRatio")) == 1.0
    assert rm.getAttribute("hypnos:bindingApproximation") == "instantaneous_1to1"
    # the binding constants are curated null -> omitted, never zero-filled
    assert not rm.hasAttribute("hypnos:bindingKd") and not rm.hasAttribute("hypnos:bindingKOn")
    targets = [li.getAttribute("hypnos:modelId") for li in rm.getElementsByTagName("rdf:li")]
    assert targets == ["nmb_agents.rocuronium.wierda_1991"]


@pytest.mark.parametrize("fmt", ["sbml", "pharmml"])
def test_xml_reversal_carries_mechanism_numbers(fmt):
    ds = hypnos.load()
    nal = export_model(fmt, ds[NALOXONE], ds)[1]
    minidom.parseString(nal)
    assert 'ntagonismKi="1.0"' in nal and "Ce50_apparent" in nal
    neo = export_model(fmt, ds[NEOSTIGMINE], ds)[1]
    minidom.parseString(neo)
    assert 'ndirectCeiling="0.5"' in neo


def test_text_exports_state_reversal_is_an_interaction():
    ds = hypnos.load()
    for fmt, prefix in (("nonmem", ";"), ("rxode2", "#"), ("pumas", "#")):
        text = export_model(fmt, ds[NEOSTIGMINE], ds)[1]
        line = next(ln for ln in text.splitlines() if "hypnos:reversalMechanism" in ln)
        assert line.startswith(prefix) and "not a solo trajectory" in line
        assert "indirect.ceiling = 0.5" in text


def test_non_reversal_model_exports_no_reversal_block():
    ds = hypnos.load()
    for fmt in _MODEL_FORMATS:
        assert "eversalMechanism" not in export_model(
            fmt, ds["hypnotics_iv.propofol.marsh_1991"], ds)[1]


def test_verify_gains_a_reversal_group():
    ds = hypnos.load()
    mv = model_verification(ds, NALOXONE)
    rev = [it for it in mv.checklist if it.group == "reversal"]
    labels = " | ".join(it.label for it in rev)
    assert "target model(s) resolve" in labels and "Ki value" in labels
    assert any("opioids.fentanyl.shafer_1990" in it.value for it in rev)
    # the renarcotization failure mode is found and quoted, not just asserted present
    fm = next(it for it in rev if "failure mode" in it.label)
    assert "renarcotization" in fm.value.lower() or "relapse" in fm.value.lower()
    assert "## Reversal coupling" in checklist_markdown(mv)


def test_verify_reversal_group_is_mechanism_specific_and_flags_gaps():
    ds = hypnos.load()
    neo = [it for it in model_verification(ds, NEOSTIGMINE).checklist if it.group == "reversal"]
    assert any("ceiling=0.5" in it.value for it in neo)
    sug = [it for it in model_verification(ds, SUGAMMADEX).checklist if it.group == "reversal"]
    assert any("molar" in it.value for it in sug)
    # sugammadex curates no failure mode yet: the checklist says so instead of staying silent
    assert any("none curated" in it.value for it in sug if "failure mode" in it.label)
    marsh = model_verification(ds, "hypnotics_iv.propofol.marsh_1991").checklist
    assert not [it for it in marsh if it.group == "reversal"]
