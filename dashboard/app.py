"""Hypnos dashboard — browse models + the model-divergence view.

    streamlit run dashboard/app.py

The killer feature (spec §6): pick a virtual patient and a dosing schedule, and
overlay the predicted plasma and effect-site curves from *every* eligible model,
greying out the ones whose envelope the patient violates and reporting the
quantitative divergence between them. This makes model-selection risk visible
and measurable. It is framed as "how much do published models disagree?" — a
research question — never "what should I give this patient?"

NOT FOR CLINICAL USE. Research / education / simulation only.
"""
from __future__ import annotations

import numpy as np

try:
    import streamlit as st
except ImportError:  # pragma: no cover
    raise SystemExit("Install the dashboard extra: pip install -e '.[dashboard]'")

import hypnos

st.set_page_config(page_title="Hypnos — model-divergence view", layout="wide")
ds = hypnos.load()

st.title("Hypnos — anesthetic PK/PD model divergence")
st.error(
    "**NOT FOR CLINICAL USE.** Research / education / simulation only. "
    "Not a dosing tool, not a TCI pump driver. This view answers "
    "*how much do published models disagree?* — never *what should I give?*"
)

with st.sidebar:
    st.header("Virtual patient")
    age = st.slider("Age (years)", 18, 95, 72)
    weight = st.slider("Weight (kg)", 30, 180, 60)
    height = st.slider("Height (cm)", 140, 210, 162)
    sex = st.radio("Sex", ["F", "M"], horizontal=True)
    st.header("Dosing")
    bolus = st.text_input("Bolus", "2 mg/kg")
    infusion = st.text_input("Infusion", "6 mg/kg/h")
    tmax = st.slider("Horizon (min)", 10, 180, 60)

patient = dict(age=age, weight=weight, height=height, sex=sex)
schedule = [("bolus", 0.0, bolus), ("infusion", 0.0, infusion)]
t = np.linspace(0, tmax, 6 * tmax + 1)

cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=schedule, t=t)

col1, col2 = st.columns([3, 2])

with col1:
    st.subheader("Effect-site concentration (µg/mL)")
    import pandas as pd

    df = pd.DataFrame({"t (min)": t})
    for r in cmp.included:
        df[r.model_id.split(".")[-1] + f" (tier {r.tier})"] = r.ce
    st.line_chart(df, x="t (min)")

with col2:
    st.subheader("Model status")
    for r in cmp.included:
        st.success(f"✓ {r.model_id}  — tier {r.tier}, Ce peak {r.ce_peak:.2f} µg/mL")
    for e in cmp.excluded:
        st.warning(f"⬚ {e['model_id']} — greyed out (tier {e['tier']}): "
                   + "; ".join(e["reasons"])[:300])
    for u in cmp.unavailable:
        st.info(f"… {u['model_id']} — {u['reason']}")

    d = cmp.divergence.get("ce", {})
    if d:
        st.metric("Peak effect-site divergence", f"{d['max_abs']:.2f} µg/mL",
                  f"{100 * d['max_rel']:.0f}% peak relative spread")

st.divider()
st.subheader("Propofol–remifentanil hypnotic synergy")
st.caption("How much does adding remifentanil deepen hypnosis at the same propofol dose? "
           "Greco response surface — coefficients illustrative (Tier C). NOT a dosing tool.")
if st.checkbox("Show interaction (propofol Schnider + remifentanil Minto)"):
    remi_sched = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]
    ir = hypnos.simulate_interaction(
        ds, "interactions.propofol_remifentanil.greco_bis",
        pk_a="hypnotics_iv.propofol.schnider_1998",
        pk_b="opioids.remifentanil.minto_1997",
        patient=patient, schedule_a=schedule, schedule_b=remi_sched, t=t,
    )
    alone = hypnos.simulate(ds, "hypnotics_iv.propofol.schnider_1998", patient=patient,
                            schedule=schedule, t=t, pd_model="pd_effect.propofol.bis_sigmoid")
    import pandas as pd

    bis = pd.DataFrame({"t (min)": t, "propofol alone": alone.effect,
                        "propofol + remifentanil": ir.effect})
    st.line_chart(bis, x="t (min)")
    c1, c2, c3 = st.columns(3)
    c1.metric("BIS min — propofol alone", f"{float(alone.effect.min()):.1f}")
    c2.metric("BIS min — + remifentanil", f"{ir.effect_min:.1f}")
    c3.metric("Composed tier", ir.tier)
