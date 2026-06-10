"""Local-anesthetic systemic-absorption subsystem — v0.6 LA0.

The safety-first entry point: systemic absorption is SITE-driven, not
milligram-driven. LA0 curates disposition + site absorption + binding only — no
toxicity thresholds (LA1). These tests cover the Bateman kernel, the
site-dominance headline, the clean separation from the IV-disposition path, and
that the LA binding (α1-AG-driven) correctly does NOT trip the v0.5 albumin caveat.
"""
import math

import numpy as np
import pytest

import hypnos
from hypnos.la import absorption_pk, concentration_at_site, site_comparison
from hypnos.models import Model

LIDO = "local_anesthetics.lidocaine.systemic"


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


# --------------------------------------------------------------------------- #
# the Bateman absorption kernel
# --------------------------------------------------------------------------- #
def test_absorption_pk_zero_at_t0():
    c = absorption_pk(100.0, ka=0.2, V_L=73.0, k10=0.0064, t=np.array([0.0]))
    assert c[0] == pytest.approx(0.0, abs=1e-12)


def test_absorption_pk_peak_time_matches_closed_form():
    ka, k10, V, dose = 0.2, 0.0064, 73.0, 100.0
    t = np.linspace(0, 120, 4801)
    c = absorption_pk(dose, ka, V, k10, t)
    tmax_numeric = t[int(np.argmax(c))]
    tmax_closed = math.log(ka / k10) / (ka - k10)   # d/dt Bateman = 0
    assert tmax_numeric == pytest.approx(tmax_closed, abs=0.1)


def test_absorption_pk_degenerate_ka_equals_k10():
    # ka == k10 limit: C(t) = (dose/V)·k·t·exp(-k·t)
    k = 0.05
    t = np.array([10.0])
    c = absorption_pk(100.0, ka=k, V_L=50.0, k10=k, t=t)
    expected = (100.0 / 50.0) * k * 10.0 * math.exp(-k * 10.0)
    assert c[0] == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# site dominance — the headline
# --------------------------------------------------------------------------- #
def test_site_comparison_sorted_by_peak_and_rank_monotone(ds):
    rows = site_comparison(ds, LIDO, dose_mg=150.0)
    assert len(rows) == 5
    # sorted highest peak first == rank order (intercostal fastest -> highest peak)
    assert [r.rank for r in rows] == [1, 2, 3, 4, 5]
    assert rows[0].site == "intercostal" and rows[-1].site == "subcutaneous"
    cmaxes = [r.cmax for r in rows]
    assert cmaxes == sorted(cmaxes, reverse=True)        # strictly the site ordering
    # the same dose gives a materially different peak by site (the safety message)
    assert rows[0].cmax > rows[-1].cmax
    # faster-absorbing sites also peak earlier
    assert rows[0].tmax_min < rows[-1].tmax_min


def test_concentration_at_site_fields(ds):
    r = concentration_at_site(ds, LIDO, site="brachial_plexus", dose_mg=100.0)
    assert r.drug == "lidocaine" and r.site == "brachial_plexus" and r.rank == 4
    assert r.ka == pytest.approx(0.05)
    assert any("Tier-C" in w for w in r.warnings)        # honest about the soft magnitude


def test_unknown_site_raises(ds):
    with pytest.raises(ValueError):
        concentration_at_site(ds, LIDO, site="spinal", dose_mg=100.0)


# --------------------------------------------------------------------------- #
# curated chemistry + the binding/albumin separation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("drug, frac, sat", [
    ("lidocaine", 0.65, False), ("bupivacaine", 0.95, True), ("ropivacaine", 0.94, True)])
def test_la_protein_binding_curated(ds, drug, frac, sat):
    pb = ds.drug(drug)["protein_binding"]
    assert pb["fraction_bound"] == frac
    assert pb["saturable"] is sat
    # LA binding is α1-AG-driven, NOT the albumin axis -> must not trip the v0.5 caveat
    assert pb["binding_sensitive"] is False
    assert ds.citation(pb["citation"]) is not None


def test_la_does_not_trip_albumin_binding_caveat(ds):
    from hypnos.simulate import evaluate_safety
    m = ds[LIDO]
    _, warnings, _ = evaluate_safety(
        m, {"age": 55, "weight": 75, "albumin_g_dl": 2.4}, ds.drug(m.drug_name))
    assert not any(w.startswith("BINDING-SENSITIVE") for w in warnings)


# --------------------------------------------------------------------------- #
# clean separation from the IV-disposition simulate()/compare() path
# --------------------------------------------------------------------------- #
def test_simulate_refuses_la_model_cleanly(ds):
    from hypnos.simulate import simulate
    with pytest.raises(NotImplementedError) as e:
        simulate(ds, LIDO, patient={"age": 40, "weight": 70},
                 schedule=[("bolus", 0.0, "100 mg")], t=np.linspace(0, 10, 11))
    assert "hypnos.la" in str(e.value)


def test_compare_handles_la_drug_gracefully(ds):
    from hypnos.simulate import compare
    cmp = compare(ds, drug="lidocaine", patient={"age": 40, "weight": 70},
                  schedule=[("bolus", 0.0, "100 mg")], t=np.linspace(0, 10, 11))
    assert cmp.included == []
    assert any("hypnos.la" in u["reason"] for u in cmp.unavailable)


# --------------------------------------------------------------------------- #
# validate
# --------------------------------------------------------------------------- #
def test_full_dataset_validates_with_la(ds):
    assert hypnos.validate_dataset(ds) == []


class _FakeDS:
    def __init__(self, model):
        self._m = model
        self.citations = {"tucker-1979-la-pharmacokinetics": {}}
        self.drugs = {}

    def __iter__(self):
        return iter([self._m])


def _la_raw(site_rates, abs_cite="tucker-1979-la-pharmacokinetics"):
    return {
        "id": "local_anesthetics.x.systemic", "subsystem": "local_anesthetics",
        "drug": {"name": "x"}, "purpose": "pk", "tier": "C",
        "primary_citation": "tucker-1979-la-pharmacokinetics",
        "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
        "parameters": [{"symbol": "V1", "value": {"central": 70.0, "units": "L"},
                        "tier": "C", "extraction": {"review_status": "unverified"}}],
        "extraction": {"review_status": "unverified"},
        "absorption": {"model": "first_order", "site_rates": site_rates,
                       "primary_citation": abs_cite},
    }


def test_validate_flags_duplicate_site_rank():
    from hypnos.validate import validate_dataset
    raw = _la_raw([{"site": "a", "rank": 1}, {"site": "b", "rank": 1}])
    probs = validate_dataset(_FakeDS(Model(raw=raw)))
    assert any("duplicate site_rates rank" in p for p in probs)


def test_validate_flags_unknown_absorption_citation():
    from hypnos.validate import validate_dataset
    raw = _la_raw([{"site": "a", "rank": 1}], abs_cite="nope-not-real")
    probs = validate_dataset(_FakeDS(Model(raw=raw)))
    assert any("absorption cites unknown" in p for p in probs)
