"""Hypnos dashboard — browse models + the model-divergence view.

    streamlit run dashboard/app.py

The killer feature (spec §6): pick a drug, a virtual patient, and a dosing
schedule; overlay the predicted plasma/effect-site curves from *every* eligible
model, grey out the ones whose envelope the patient violates, and report the
quantitative divergence. It answers "how much do published models disagree?" — a
research question — never "what should I give this patient?"

All compute is the tested package logic (compare, time_to_peak_effect,
simulate_interaction) over `hypnos.presets` doses, so the dashboard never drifts
from the CLI. NOT FOR CLINICAL USE. Research / education / simulation only.
"""
from __future__ import annotations

import numpy as np

try:
    import pandas as pd
    import streamlit as st
except ImportError:  # pragma: no cover
    raise SystemExit("Install the dashboard extra: pip install -e '.[dashboard]'")

import hypnos
from hypnos.filter import pk_drugs
from hypnos.presets import default_schedule_for
from hypnos.reference import alveolar_washin, alveolar_washout, mac_age_corrected

st.set_page_config(page_title="Hypnos — model-divergence view", layout="wide")
ds = hypnos.load()

st.title("Hypnos — anesthetic PK/PD model divergence")
st.error(
    "**NOT FOR CLINICAL USE.** Research / education / simulation only. "
    "Not a dosing tool, not a TCI pump driver. This view answers "
    "*how much do published models disagree?* — never *what should I give?*"
)

drugs = pk_drugs(ds)

with st.sidebar:
    st.header("Drug")
    drug = st.selectbox("Drug", drugs, index=drugs.index("propofol") if "propofol" in drugs else 0)
    st.header("Virtual patient")
    age = st.slider("Age (years)", 1, 95, 72)
    weight = st.slider("Weight (kg)", 5, 180, 60)
    height = st.slider("Height (cm)", 90, 210, 162)
    sex = st.radio("Sex", ["F", "M"], horizontal=True)
    st.header("Dosing")
    # default doses are drug-appropriate (a 2 mg/kg propofol dose would be a
    # ~1000x overdose for remifentanil) — shared with the CLI via hypnos.presets
    preset = default_schedule_for(drug)
    bolus = next((spec for kind, _, spec in preset if kind == "bolus"), "")
    infusion = next((spec for kind, _, spec in preset if kind == "infusion"), "")
    bolus = st.text_input("Bolus (blank = none)", bolus)
    infusion = st.text_input("Infusion (blank = none)", infusion)
    tmax = st.slider("Horizon (min)", 10, 180, 60)

patient = dict(age=age, weight=weight, height=height, sex=sex)
schedule = []
if bolus.strip():
    schedule.append(("bolus", 0.0, bolus.strip()))
if infusion.strip():
    schedule.append(("infusion", 0.0, infusion.strip()))
t = np.linspace(0, tmax, 6 * tmax + 1)

cmp = hypnos.compare(ds, drug=drug, patient=patient, schedule=schedule, t=t)
unit = cmp.concentration_unit
f = cmp.conc_factor

col1, col2 = st.columns([3, 2])

with col1:
    st.subheader(f"Effect-site concentration ({unit})")
    df = pd.DataFrame({"t (min)": t})
    any_ce = False
    for r in cmp.included:
        if float(np.max(r.ce)) > 0:  # only models with a ke0 link have an effect-site curve
            df[r.model_id.split(".")[-1] + f" (tier {r.tier})"] = r.ce * f
            any_ce = True
    if any_ce:
        st.line_chart(df, x="t (min)")
    st.subheader(f"Plasma concentration ({unit})")
    dfp = pd.DataFrame({"t (min)": t})
    for r in cmp.included:
        dfp[r.model_id.split(".")[-1]] = r.cp * f
    st.line_chart(dfp, x="t (min)")

with col2:
    st.subheader("Model status")
    for r in cmp.included:
        mdape = sorted({e["value"] for e in ds[r.model_id].predictive_mdape})
        acc = (f", published MDAPE {mdape[0]:g}%" if len(mdape) == 1
               else f", published MDAPE {mdape[0]:g}–{mdape[-1]:g}%" if mdape else "")
        st.success(f"✓ {r.model_id} — tier {r.tier}, Cp peak {r.cp_peak_display:.2f} {unit}{acc}")
    for e in cmp.excluded:
        st.warning(f"⬚ {e['model_id']} — greyed out (tier {e['tier']}): "
                   + "; ".join(e["reasons"])[:300])
    for u in cmp.unavailable:
        st.info(f"… {u['model_id']} — {u['reason']}")

    d = cmp.divergence.get("ce") or cmp.divergence.get("cp") or {}
    if d:
        st.metric("Peak divergence across models", f"{d['max_abs'] * f:.2f} {unit}",
                  f"{100 * d['max_rel']:.0f}% peak relative spread")
        drv = d.get("driver")
        if drv:
            st.caption(f"Driver of the disagreement: **{drv['high'].split('.')[-1]}** vs "
                       f"**{drv['low'].split('.')[-1]}** (furthest apart at the peak instant).")

    # Onset (time-to-peak-effect) for models with a ke0 link
    onset = []
    for r in cmp.included:
        try:
            pe = hypnos.time_to_peak_effect(ds, r.model_id, patient=patient)
            onset.append({"model": r.model_id.split(".")[-1], "tpeak (min)": round(pe.tpeak_min, 2),
                          "ke0 (1/min)": round(pe.ke0, 3)})
        except (ValueError, NotImplementedError):
            pass
    if onset:
        st.subheader("Onset (time to peak effect)")
        st.dataframe(pd.DataFrame(onset), hide_index=True)

if drug == "propofol":
    st.divider()
    st.subheader("Propofol–remifentanil hypnotic synergy")
    st.caption("How much does adding remifentanil deepen hypnosis at the same propofol dose? "
               "Greco response surface — coefficients illustrative (Tier C). NOT a dosing tool.")
    if st.checkbox("Show interaction (propofol Schnider + remifentanil Minto)"):
        prop_sched = schedule or default_schedule_for("propofol")
        remi_sched = default_schedule_for("remifentanil")
        ir = hypnos.simulate_interaction(
            ds, "interactions.propofol_remifentanil.greco_bis",
            pk_a="hypnotics_iv.propofol.schnider_1998",
            pk_b="opioids.remifentanil.minto_1997",
            patient=patient, schedule_a=prop_sched, schedule_b=remi_sched, t=t,
        )
        alone = hypnos.simulate(ds, "hypnotics_iv.propofol.schnider_1998", patient=patient,
                                schedule=prop_sched, t=t, pd_model="pd_effect.propofol.bis_sigmoid")
        bis = pd.DataFrame({"t (min)": t, "propofol alone": alone.effect,
                            "propofol + remifentanil": ir.effect})
        st.line_chart(bis, x="t (min)")
        c1, c2, c3 = st.columns(3)
        c1.metric("BIS min — propofol alone", f"{float(alone.effect.min()):.1f}")
        c2.metric("BIS min — + remifentanil", f"{ir.effect_min:.1f}")
        c3.metric("Composed tier", ir.tier)

# --- Inhalational agents: a different (non-compartmental) parameter convention ---
# Volatiles are characterized by MAC + its age correction + the blood:gas-driven
# wash-in, not compartmental PK — so they live outside the divergence view above.
# All numbers come from the same tested package functions the CLI uses.
volatiles = sorted(m.id for m in ds if m.purpose == "physicochemical")
if volatiles:
    st.divider()
    st.subheader("Inhalational agents (volatiles)")
    st.caption("Not compartmental: characterized by MAC (minimum alveolar concentration), "
               "its age correction, and the blood:gas-driven wash-in / wash-out. "
               "Research / education only.")
    vc1, vc2 = st.columns(2)
    with vc1:
        st.markdown(f"**Age-corrected MAC at age {age} y**")
        rows = []
        for vid in volatiles:
            r = hypnos.mac(ds, vid, age=age)
            rows.append({"agent": vid.split(".")[1], "MAC (vol%)": round(r.mac_age, 2),
                         "MAC-awake": round(r.mac_awake_age, 2), "blood:gas λ": r.blood_gas})
        st.dataframe(pd.DataFrame(rows), hide_index=True)
        ages = np.arange(1, 91)
        mac_curve = pd.DataFrame({"age (y)": ages})
        for vid in volatiles:
            mac40 = next(p.central for p in ds[vid].parameters if p.symbol == "MAC40")
            mac_curve[vid.split(".")[1]] = [mac_age_corrected(mac40, a) for a in ages]
        st.caption("MAC vs age (Mapleson/Nickalls 2003, ~6% per decade)")
        st.line_chart(mac_curve, x="age (y)")
    with vc2:
        tw = np.linspace(0, 8, 160)
        st.markdown("**Wash-in FA/FI — lower blood:gas λ washes in faster**")
        washin = pd.DataFrame({"t (min)": tw})
        for vid in volatiles:
            lam = next(p.central for p in ds[vid].parameters if p.symbol == "blood_gas")
            washin[vid.split(".")[1]] = alveolar_washin(lam, tw)[0]
        st.line_chart(washin, x="t (min)")
        st.markdown("**Wash-out FA/FA₀ (emergence) — lower λ washes out faster**")
        washout = pd.DataFrame({"t (min)": tw})
        for vid in volatiles:
            lam = next(p.central for p in ds[vid].parameters if p.symbol == "blood_gas")
            washout[vid.split(".")[1]] = alveolar_washout(lam, tw)[0]
        st.line_chart(washout, x="t (min)")
        st.caption("Single-compartment alveolar model, standard ventilation. "
                   "Comparative only — NOT a per-patient predictor.")
