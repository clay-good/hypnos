import numpy as np
import pytest

import hypnos

DEX = "alpha2_agonists.dexmedetomidine.hannivoort_2015"
PAEDFUSOR = "hypnotics_iv.propofol.paedfusor_2005"
KATARIA = "hypnotics_iv.propofol.kataria_1994"
SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
MARSH = "hypnotics_iv.propofol.marsh_1991"
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
FENTANYL = "opioids.fentanyl.shafer_1990"

T = np.linspace(0, 60, 361)
PED_SCHED = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "10 mg/kg/h")]
CHILD = dict(age=6, weight=20, height=115, sex="M")
ADULT = dict(age=40, weight=70, height=175, sex="M")


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_dexmedetomidine_allometric_sane(ds):
    dex = [("infusion", 0.0, "6 mcg/kg/h"), ("infusion", 10.0, "0.5 mcg/kg/h")]
    res = hypnos.simulate(ds, DEX, patient=ADULT, schedule=dex, t=T)
    assert res.tier == "B"
    # plasma should be in the ng/mL range (ug/mL ~ 1e-3) for a clinical loading regimen
    assert 0.0005 < res.cp_peak < 0.01


def test_dexmedetomidine_allometric_scaling(ds):
    # central volume scales linearly with weight; double the weight -> ~half the conc
    light = hypnos.simulate(ds, DEX, patient=dict(age=40, weight=50, height=170, sex="M"),
                            schedule=[("bolus", 0.0, "50 mcg")], t=T)
    heavy = hypnos.simulate(ds, DEX, patient=dict(age=40, weight=100, height=180, sex="M"),
                            schedule=[("bolus", 0.0, "50 mcg")], t=T)
    assert heavy.cp[0] < light.cp[0]


def test_paedfusor_pediatric_bolus(ds):
    # 2 mg/kg into a 20 kg child = 40 mg; V1 = 0.4584*20 = 9.168 L -> Cp(0)=40/9.168
    res = hypnos.simulate(ds, PAEDFUSOR, patient=CHILD,
                          schedule=[("bolus", 0.0, "2 mg/kg")], t=T)
    assert res.tier == "B"
    assert res.warnings == []
    assert abs(res.cp[0] - 40.0 / (0.4584 * 20)) < 1e-6


def test_kataria_reference_child_parameters(ds):
    # the kernel reproduces the standard published Kataria parameters at a 23 kg, 7 y
    # reference child: V1 = 0.41*23 = 9.43 L, V2 = 0.78*23 + 3.1*7 - 16 = 23.64 L
    from hypnos.export.registry import propofol_kataria_1994
    p = propofol_kataria_1994(dict(age=7, weight=23))
    assert abs(p.V1 - 9.43) < 1e-9
    assert abs(p.derived["V2"] - 23.64) < 1e-9
    assert p.ke0 == 0.0          # PK-only (no published ke0)


def test_kataria_paedfusor_diverge_for_a_child(ds):
    # the named pediatric pair: both in-envelope for a child, but they disagree on
    # plasma — exactly the "Kataria vs Paedfusor" question the literature asks.
    child = dict(age=7, weight=23, height=122, sex="M")
    cmp = hypnos.compare(ds, drug="propofol", patient=child, schedule=PED_SCHED, t=T)
    included = {r.model_id for r in cmp.included}
    assert KATARIA in included and PAEDFUSOR in included
    kat = next(r for r in cmp.included if r.model_id == KATARIA)
    pae = next(r for r in cmp.included if r.model_id == PAEDFUSOR)
    assert abs(kat.cp_peak - pae.cp_peak) > 0.1     # a real, measurable divergence
    assert cmp.divergence["cp"]["max_rel"] > 0      # plasma spread is reported


def test_kataria_in_adult_is_extrapolation(ds):
    res = hypnos.simulate(ds, KATARIA, patient=ADULT, schedule=PED_SCHED, t=T)
    assert res.tier == "D"
    assert res.excluded
    assert any("EXTRAPOLATION" in w for w in res.warnings)


def test_adult_model_in_child_is_pediatric_extrapolation(ds):
    res = hypnos.simulate(ds, SCHNIDER, patient=CHILD, schedule=PED_SCHED, t=T)
    assert res.tier == "D"
    assert res.excluded
    assert any("PEDIATRIC EXTRAPOLATION" in w for w in res.warnings)


def test_pediatric_model_in_adult_is_extrapolation(ds):
    res = hypnos.simulate(ds, PAEDFUSOR, patient=ADULT, schedule=PED_SCHED, t=T)
    assert res.tier == "D"
    assert any("EXTRAPOLATION" in w for w in res.warnings)


def test_geriatric_extrapolation_label(ds):
    # Dexmedetomidine envelope tops out at 70 y; a 90 y patient -> geriatric extrapolation
    old = dict(age=90, weight=70, height=170, sex="M")
    res = hypnos.simulate(ds, DEX, patient=old, schedule=[("bolus", 0.0, "50 mcg")], t=T)
    assert res.tier == "D"
    assert any("GERIATRIC EXTRAPOLATION" in w for w in res.warnings)


def test_pediatric_compare_greys_adult_models(ds):
    cmp = hypnos.compare(ds, drug="propofol", patient=CHILD, schedule=PED_SCHED, t=T)
    included = {r.model_id for r in cmp.included}
    # Paedfusor (pediatric) and Eleveld (broad-envelope, covers children) are in-envelope;
    # the adult-only Marsh and Schnider are greyed out as pediatric extrapolations.
    assert PAEDFUSOR in included and ELEVELD in included
    excluded = {e["model_id"] for e in cmp.excluded}
    assert SCHNIDER in excluded and MARSH in excluded


def test_fentanyl_kernel_pending_refuses(ds):
    with pytest.raises(NotImplementedError):
        hypnos.simulate(ds, FENTANYL, patient=ADULT, schedule=PED_SCHED, t=T)


def test_new_subsystems_present(ds):
    s = hypnos.summary(ds)
    assert s["by_subsystem"].get("alpha2_agonists") == 1
    assert s["n_models"] >= 9
    assert s["kernels_implemented"] >= 7
