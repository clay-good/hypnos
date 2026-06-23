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
