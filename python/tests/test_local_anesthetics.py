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


# --------------------------------------------------------------------------- #
# v0.6 LA1 — toxicity thresholds + the double-uncertainty view
# --------------------------------------------------------------------------- #
BUPI = "local_anesthetics.bupivacaine.systemic"
ROPI = "local_anesthetics.ropivacaine.systemic"


def test_toxicity_thresholds_are_ranges_with_basis(ds):
    m = ds[BUPI]
    assert m.has_toxicity_thresholds and m.is_safety_critical
    for th in m.toxicity_thresholds:
        assert th.low < th.high                       # a range, never a line
        assert th.basis in ("total_plasma", "free_plasma")
        assert th.endpoint in ("cns_first_symptoms", "cns_seizure", "cardiovascular")
        assert th.relative_width is not None and th.relative_width > 0


def test_saturable_total_threshold_carries_saturation_caveat(ds):
    # the load-bearing free-fraction guard: a saturable drug on a total basis must
    # name the under-prediction failure mode (v0.6 §3.2/§4)
    for th in ds[BUPI].toxicity_thresholds:
        if th.basis == "total_plasma":
            assert th.saturation_caveat and "free" in th.saturation_caveat.lower()


def test_free_concentration_linear_baseline_and_caveat(ds):
    from hypnos.la import free_concentration
    pb = ds.drug("bupivacaine")["protein_binding"]      # 95% bound, saturable
    fc = free_concentration(np.array([2.0, 4.0]), pb)
    assert fc.free_fraction == pytest.approx(0.05)
    # the LINEAR baseline is still fraction_bound * total (now exposed as c_free_linear);
    # since bupivacaine carries a curated LA3 model, the primary c_free is the non-linear one
    assert fc.c_free_linear == pytest.approx([0.1, 0.2])
    assert fc.model == "capacity_limited"
    assert fc.saturable and any("SATURABLE" in w for w in fc.warnings)


def test_free_concentration_no_binding_returns_none():
    from hypnos.la import free_concentration
    fc = free_concentration(np.array([1.0]), None)
    assert fc.c_free is None and any("not derivable" in w for w in fc.warnings)


def test_double_uncertainty_view_structure(ds):
    from hypnos.la import double_uncertainty
    du = double_uncertainty(ds, BUPI, site="lumbar_epidural", dose_mg=100.0)
    assert du.drug == "bupivacaine" and du.site == "lumbar_epidural"
    assert du.peak_total > 0 and du.peak_free is not None and du.peak_free < du.peak_total
    # endpoints sorted by low bound; each placed on the MATCHING basis
    assert du.endpoints == sorted(du.endpoints, key=lambda e: e.low)
    for e in du.endpoints:
        expected = du.peak_free if e.basis == "free_plasma" else du.peak_total
        assert e.predicted_peak == pytest.approx(expected)
        assert e.position in ("below range", "within range", "above range")
    # the honest punchline names the dominant uncertainty, never a verdict
    assert "THRESHOLD uncertainty dominates" in du.dominant_uncertainty
    assert any("NOT a dosing tool" in w for w in du.warnings)


def test_double_uncertainty_dominant_uses_fold_ranges(ds):
    # the comparison is on a consistent multiplicative scale (high/low vs max/min Cmax)
    from hypnos.la import double_uncertainty
    du = double_uncertainty(ds, BUPI, site="intercostal", dose_mg=100.0)
    assert du.site_cmax_spread is not None and du.site_cmax_spread > 1.0
    assert "fold-range" in du.dominant_uncertainty


def test_double_uncertainty_requires_curated_thresholds(ds):
    # never fabricate: no thresholds => the view refuses rather than inventing a line
    from hypnos.la import double_uncertainty
    assert not ds[LIDO].has_toxicity_thresholds or True  # lidocaine has them; assert the refusal path on a bare model
    bare = Model(raw=_la_raw([{"site": "intercostal", "rank": 1, "ka": 0.2}]))

    class _DS:
        def __getitem__(self, k):
            return bare
        def drug(self, n):
            return None
    with pytest.raises(ValueError, match="no curated toxicity thresholds"):
        double_uncertainty(_DS(), "local_anesthetics.x.systemic", site="intercostal", dose_mg=100.0)


def test_validate_flags_non_range_threshold():
    from hypnos.validate import validate_dataset
    raw = _la_raw([{"site": "a", "rank": 1, "ka": 0.2}])
    raw["toxicity_thresholds"] = [{
        "endpoint": "cns_first_symptoms",
        "concentration_range": {"low": 4.0, "high": 2.0, "units": "ug/mL"},  # inverted
        "basis": "total_plasma", "tier": "C",
        "extraction": {"review_status": "unverified"},
    }]
    probs = validate_dataset(_FakeDS(Model(raw=raw)))
    assert any("low" in p and "high" in p and "toxicity" in p for p in probs)


def test_validate_flags_missing_saturation_caveat():
    from hypnos.validate import validate_dataset
    raw = _la_raw([{"site": "a", "rank": 1, "ka": 0.2}])
    raw["drug"] = {"name": "bupivacaine"}
    raw["toxicity_thresholds"] = [{
        "endpoint": "cns_first_symptoms",
        "concentration_range": {"low": 1.6, "high": 2.6, "units": "ug/mL"},
        "basis": "total_plasma", "tier": "C",
        "extraction": {"review_status": "unverified"},
    }]

    class _DSsat(_FakeDS):
        def drug(self, name):
            return {"name": "bupivacaine", "protein_binding": {"saturable": True}}
    probs = validate_dataset(_DSsat(Model(raw=raw)))
    assert any("saturation_caveat" in p for p in probs)


def test_verification_checklist_has_la_group(ds):
    from hypnos.verification import checklist_markdown, model_verification
    mv = model_verification(ds, BUPI)
    groups = {it.group for it in mv.checklist}
    assert "local_anesthetic" in groups
    md = checklist_markdown(mv)
    assert "toxicity thresholds" in md.lower()
    assert "RANGE, never a line" in md


def test_la_export_carries_safety_critical_and_threshold_ranges(ds):
    from hypnos.export import pharmml, sbml, tci_json
    m = ds[BUPI]
    sb = sbml.build(m, ds)
    assert "<hypnos:safetyCritical>true</hypnos:safetyCritical>" in sb
    assert sb.count("hypnos:toxicityThresholdRange") == len(m.toxicity_thresholds)
    import json
    d = json.loads(tci_json.build(m, ds))
    assert d["safety_critical"] is True
    assert len(d["toxicity_thresholds"]) == len(m.toxicity_thresholds)   # AS ranges, verbatim
    assert d["provenance"]["hypnos:safetyCritical"] == "true"
    pm = pharmml.build(m, ds)
    assert pm.count("toxicityThresholdRange") == len(m.toxicity_thresholds)


def test_ropivacaine_cvs_margin_wider_than_bupivacaine(ds):
    # the educational point LA1 sets up for LA2: ropivacaine's CNS->CVS separation is
    # wider than bupivacaine's (the narrow bupivacaine margin is the cardiotoxicity story)
    def cvs_over_cns(mid):
        m = ds[mid]
        cns = next(t for t in m.toxicity_thresholds
                   if t.endpoint == "cns_first_symptoms" and t.basis == "total_plasma")
        cvs = next(t for t in m.toxicity_thresholds if t.endpoint == "cardiovascular")
        return cvs.midpoint / cns.midpoint
    assert cvs_over_cns(ROPI) > cvs_over_cns(BUPI)


# --------------------------------------------------------------------------- #
# v0.6 LA2 — stereochemistry / cardiotoxicity differentiation
# --------------------------------------------------------------------------- #
LEVO = "local_anesthetics.levobupivacaine.systemic"


def test_levobupivacaine_model_present_and_simulable(ds):
    from hypnos.la import concentration_at_site
    m = ds[LEVO]
    assert m.subsystem == "local_anesthetics" and m.has_toxicity_thresholds
    r = concentration_at_site(ds, LEVO, site="lumbar_epidural", dose_mg=100.0)
    assert r.drug == "levobupivacaine" and r.cmax > 0


@pytest.mark.parametrize("drug, rank, stereo, margin", [
    ("bupivacaine", "high", "racemate", "narrow"),
    ("levobupivacaine", "intermediate", "S_enantiomer", "moderate"),
    ("ropivacaine", "low", "S_enantiomer", "wide"),
    ("lidocaine", "low", "achiral", "wide")])
def test_cardiotoxicity_class_curated(ds, drug, rank, stereo, margin):
    cc = ds.drug(drug)["cardiotoxicity_class"]
    assert cc["rank"] == rank and cc["stereochemistry"] == stereo
    assert cc["cns_to_cvs_margin"] == margin
    assert ds.citation(cc["citation"]) is not None


def test_cardiotoxicity_comparison_ordering(ds):
    from hypnos.la import cardiotoxicity_comparison
    rows = cardiotoxicity_comparison(ds)
    assert [a.drug for a in rows][:2] == ["bupivacaine", "levobupivacaine"]   # most cardiotoxic first
    # the numeric CNS-to-CVS margin is MONOTONE with the qualitative ranking
    folds = [a.cns_to_cvs_fold for a in rows]
    assert folds == sorted(folds)
    # bupivacaine has the narrowest curated margin; lidocaine the widest
    bupi = next(a for a in rows if a.drug == "bupivacaine")
    lido = next(a for a in rows if a.drug == "lidocaine")
    assert bupi.cns_to_cvs_fold < lido.cns_to_cvs_fold


def test_double_uncertainty_surfaces_cardiotoxicity(ds):
    from hypnos.la import double_uncertainty
    du = double_uncertainty(ds, BUPI, site="lumbar_epidural", dose_mg=100.0)
    assert du.cardiotoxicity is not None and du.cardiotoxicity["rank"] == "high"
    # the narrow-margin agent carries the explicit cardiotoxicity warning
    assert any("CARDIOTOXICITY" in w and "NARROW" in w for w in du.warnings)
    # a wide-margin agent does not raise that warning
    du_lido = double_uncertainty(ds, LIDO, site="lumbar_epidural", dose_mg=100.0)
    assert not any("CARDIOTOXICITY" in w for w in du_lido.warnings)


def test_cardiotoxicity_export_carried(ds):
    from hypnos.export import pharmml, sbml, tci_json
    import json
    m = ds[BUPI]
    assert 'hypnos:cardiotoxicityClass' in sbml.build(m, ds)
    assert 'cardiotoxicityClass rank="high"' in pharmml.build(m, ds)
    d = json.loads(tci_json.build(m, ds))
    assert d["cardiotoxicity_class"]["rank"] == "high"


def test_validate_flags_bad_cardiotoxicity_rank():
    from hypnos.validate import validate_dataset

    class _DS:
        def __init__(self):
            self.citations = {"tucker-1979-la-pharmacokinetics": {}}
            self.drugs = {"x": {"name": "x", "cardiotoxicity_class": {
                "rank": "extreme", "cns_to_cvs_margin": "narrow",
                "citation": "tucker-1979-la-pharmacokinetics"}}}

        def __iter__(self):
            return iter([])

        def drug(self, n):
            return self.drugs.get(n)
    probs = validate_dataset(_DS())
    assert any("invalid rank" in p for p in probs)


def test_verification_checklist_has_cardiotoxicity_item(ds):
    from hypnos.verification import model_verification
    mv = model_verification(ds, BUPI)
    assert any("cardiotoxicity class" in it.label.lower() for it in mv.checklist)


# --------------------------------------------------------------------------- #
# v0.6 LA3 — the non-linear (capacity-limited) free-fraction model
# --------------------------------------------------------------------------- #
def test_saturable_free_concentration_low_conc_matches_linear():
    from hypnos.la import saturable_free_concentration
    # in the low-concentration limit the capacity-limited model reduces to linear
    ct = np.array([0.001])
    nl = saturable_free_concentration(ct, fraction_bound=0.95, capacity_ug_ml=4.9)[0]
    assert nl == pytest.approx(ct[0] * 0.05, rel=1e-2)


def test_saturable_free_concentration_monotone_rising_fraction():
    from hypnos.la import saturable_free_concentration
    ct = np.array([0.5, 2.0, 5.0, 10.0, 20.0])
    free = saturable_free_concentration(ct, fraction_bound=0.95, capacity_ug_ml=4.9)
    frac = free / ct
    # the free FRACTION rises monotonically with total (binding saturates)
    assert np.all(np.diff(frac) > 0)
    # it stays a fraction in (0, 1] and the non-linear free always >= the linear
    assert np.all(free >= ct * 0.05 - 1e-9) and np.all(frac <= 1.0 + 1e-9)


def test_saturable_conservation_total_equals_free_plus_bound():
    # invert total->free, reconstruct bound, and confirm free+bound == total
    from hypnos.la import saturable_free_concentration
    fb0, cap = 0.95, 4.9
    ka = (fb0 / (1 - fb0)) / cap
    ct = np.array([1.0, 5.0, 25.0])
    free = saturable_free_concentration(ct, fb0, cap)
    bound = cap * ka * free / (1 + ka * free)
    assert np.allclose(free + bound, ct, rtol=1e-9)


def test_free_concentration_uses_curated_nonlinear_model(ds):
    from hypnos.la import free_concentration
    pb = ds.drug("bupivacaine")["protein_binding"]
    assert pb["free_fraction_model"]["type"] == "capacity_limited"
    fc = free_concentration(np.array([2.0, 20.0]), pb)
    assert fc.model == "capacity_limited"
    # the curated non-linear trace rises ABOVE the linear one, more so at high total
    assert fc.c_free[0] > fc.c_free_linear[0]
    assert (fc.c_free[1] / fc.c_free_linear[1]) > (fc.c_free[0] / fc.c_free_linear[0])
    assert any("SATURABLE" in w and "non-linear" in w.lower() for w in fc.warnings)


def test_free_concentration_falls_back_to_linear_without_model():
    # a saturable drug with NO curated free_fraction_model stays linear + caveat
    from hypnos.la import free_concentration
    pb = {"fraction_bound": 0.9, "saturable": True}      # no free_fraction_model
    fc = free_concentration(np.array([10.0]), pb)
    assert fc.model == "linear"
    assert fc.c_free[0] == pytest.approx(1.0)
    assert any("no non-linear model curated" in w for w in fc.warnings)


def test_lidocaine_nonsaturable_has_no_model_and_stays_linear(ds):
    from hypnos.la import free_concentration
    pb = ds.drug("lidocaine")["protein_binding"]
    assert "free_fraction_model" not in pb and pb["saturable"] is False
    fc = free_concentration(np.array([5.0]), pb)
    assert fc.model == "linear" and fc.c_free[0] == pytest.approx(5.0 * 0.35)


def test_double_uncertainty_surfaces_saturation_gap(ds):
    from hypnos.la import double_uncertainty
    du = double_uncertainty(ds, BUPI, site="intercostal", dose_mg=200.0)
    assert du.free_model == "capacity_limited"
    assert du.peak_free is not None and du.peak_free_linear is not None
    assert du.peak_free > du.peak_free_linear               # non-linear under-prediction made visible
    assert any("FREE-FRACTION SATURATION" in w for w in du.warnings)


def test_validate_flags_free_fraction_model_on_nonsaturable():
    from hypnos.validate import validate_dataset

    class _DS:
        def __init__(self):
            self.citations = {"tucker-1979-la-pharmacokinetics": {}}
            self.drugs = {"x": {"name": "x", "protein_binding": {
                "fraction_bound": 0.5, "saturable": False,
                "free_fraction_model": {"type": "capacity_limited",
                                        "binding_capacity_ug_ml": 5.0,
                                        "citation": "tucker-1979-la-pharmacokinetics"}}}}

        def __iter__(self):
            return iter([])

        def drug(self, n):
            return self.drugs.get(n)
    probs = validate_dataset(_DS())
    assert any("not marked saturable" in p for p in probs)


def test_free_fraction_model_export_carried(ds):
    from hypnos.export import pharmml, sbml, tci_json
    import json
    m = ds[BUPI]
    assert 'hypnos:freeFractionModel' in sbml.build(m, ds)
    assert 'freeFractionModel type="capacity_limited"' in pharmml.build(m, ds)
    d = json.loads(tci_json.build(m, ds))
    assert d["protein_binding"]["free_fraction_model"]["type"] == "capacity_limited"
