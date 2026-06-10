"""Hypnos dashboard — browse models + the model-divergence view.

    streamlit run dashboard/app.py

The killer feature (spec §6): pick a drug, a virtual patient, and a dosing
schedule; overlay the predicted plasma/effect-site curves from *every* eligible
model, grey out the ones whose envelope the patient violates, and report the
quantitative divergence. It answers "how much do published models disagree?" — a
research question — never "what should I give this patient?"

With "prediction bands" enabled (v0.2 §7), each band-eligible model also draws a
seeded 5–95% ribbon, and the view answers the sharper question — *are the models
distinguishable beyond their own stated variability?* (separation index) and
*what dominates the uncertainty here?* (variance decomposition). Models with no
published BSV stay bare lines and are named, never given a fabricated band.

All compute is the tested package logic (compare, time_to_peak_effect,
simulate_interaction) over `hypnos.presets` doses, so the dashboard never drifts
from the CLI. NOT FOR CLINICAL USE. Research / education / simulation only.
"""
from __future__ import annotations

import numpy as np

try:
    import altair as alt
    import pandas as pd
    import streamlit as st
except ImportError:  # pragma: no cover
    raise SystemExit("Install the dashboard extra: pip install -e '.[dashboard]'")

import hypnos
from hypnos.filter import pk_drugs
from hypnos.presets import default_schedule_for
from hypnos.reference import alveolar_washin, alveolar_washout, mac_age_corrected


def band_chart(t, results, key, f, unit, percentile=(5, 95)):
    """Layered Altair chart of per-model curves with seeded prediction ribbons (v0.2 §7).

    Band-eligible models (those that publish between-subject variability) get a
    shaded percentile ribbon + a solid median line; models with no curated Ω are
    drawn as a dashed point-estimate line and named — never given a fabricated
    band (the never-synthesize rule, carried into the view). Returns ``None`` when
    no model has a curve for ``key`` (e.g. effect-site for a PK-only comparison).
    """
    lo, hi = percentile
    band_rows, line_rows = [], []
    for r in results:
        name = r.model_id.split(".")[-1] + f" (tier {r.tier})"
        q = getattr(r, f"{key}_quantiles")
        point = getattr(r, key)
        if q is not None:  # band-eligible: ribbon + median
            for i, tt in enumerate(t):
                band_rows.append({"t (min)": float(tt), "model": name,
                                  "lo": float(q[lo][i] * f), "hi": float(q[hi][i] * f)})
                line_rows.append({"t (min)": float(tt), "model": name,
                                  unit: float(q[50][i] * f), "kind": "median"})
        elif float(point.max()) > 0:  # point-estimate line only
            for i, tt in enumerate(t):
                line_rows.append({"t (min)": float(tt), "model": name,
                                  unit: float(point[i] * f), "kind": "point (no published BSV)"})
    if not line_rows:
        return None
    layers = []
    if band_rows:
        layers.append(
            alt.Chart(pd.DataFrame(band_rows)).mark_area(opacity=0.18).encode(
                x=alt.X("t (min):Q"),
                y=alt.Y("lo:Q", title=unit), y2="hi:Q",
                color=alt.Color("model:N", legend=alt.Legend(title=f"model (band = {lo}–{hi}%)")),
            )
        )
    layers.append(
        alt.Chart(pd.DataFrame(line_rows)).mark_line().encode(
            x=alt.X("t (min):Q"), y=alt.Y(f"{unit}:Q", title=unit), color="model:N",
            strokeDash=alt.StrokeDash("kind:N", legend=alt.Legend(title="line")),
        )
    )
    return alt.layer(*layers).resolve_scale(color="shared").properties(height=340)

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
    st.header("Organ function (v0.5)")
    st.caption("Declare organ impairment to watch the eligible-model set shrink. A model "
               "with no cited standing in that state is greyed + Tier-D (extrapolation); "
               "remifentanil keeps standing (esterase clearance). Defaults = normal.")
    child_pugh = st.selectbox("Child-Pugh (hepatic)", ["— none —", "A", "B", "C"], index=0)
    crcl = st.slider("CrCl (mL/min, renal)", 5, 120, 100,
                     help="renal impairment (KDIGO stage ≥3) below 60")
    albumin = st.slider("Albumin (g/dL)", 1.5, 5.0, 4.0, step=0.1,
                        help="hypoalbuminemia below 3.5 → raised free fraction")
    ejection_fraction = st.slider("Ejection fraction (%)", 10, 70, 60,
                                  help="low cardiac output below 40")
    st.header("Dosing")
    # default doses are drug-appropriate (a 2 mg/kg propofol dose would be a
    # ~1000x overdose for remifentanil) — shared with the CLI via hypnos.presets
    preset = default_schedule_for(drug)
    bolus = next((spec for kind, _, spec in preset if kind == "bolus"), "")
    infusion = next((spec for kind, _, spec in preset if kind == "infusion"), "")
    bolus = st.text_input("Bolus (blank = none)", bolus)
    infusion = st.text_input("Infusion (blank = none)", infusion)
    tmax = st.slider("Horizon (min)", 10, 180, 60)
    st.header("Prediction bands (v0.2)")
    show_bands = st.checkbox(
        "Seeded 5–95% prediction bands", value=False,
        help="Monte-Carlo over each model's curated between-subject variability (Ω) "
             "and residual error (Σ). Only models that publish random effects earn a "
             "band; the rest stay bare lines (never a fabricated ribbon). Seeded, so "
             "the quantiles are byte-reproducible. A band is a statement about the "
             "model's stated uncertainty — NOT a claim about a real patient.")
    seed = st.number_input("Seed", 0, 99999, 7, disabled=not show_bands)
    samples = st.select_slider("Monte-Carlo samples", [500, 1000, 2000, 4000], 2000,
                               disabled=not show_bands)

patient = dict(age=age, weight=weight, height=height, sex=sex,
               crcl_ml_min=crcl, albumin_g_dl=albumin, ejection_fraction_pct=ejection_fraction)
if child_pugh != "— none —":
    patient["child_pugh"] = child_pugh
schedule = []
if bolus.strip():
    schedule.append(("bolus", 0.0, bolus.strip()))
if infusion.strip():
    schedule.append(("infusion", 0.0, infusion.strip()))
t = np.linspace(0, tmax, 6 * tmax + 1)

cmp = hypnos.compare(ds, drug=drug, patient=patient, schedule=schedule, t=t,
                     bands=bool(show_bands), percentile=(5, 95),
                     samples=int(samples), seed=int(seed) if show_bands else None)
unit = cmp.concentration_unit
f = cmp.conc_factor

col1, col2 = st.columns([3, 2])

with col1:
    st.subheader(f"Effect-site concentration ({unit})")
    if show_bands:
        # shaded 5–95% ribbons (band-eligible models) + bare lines (the rest)
        ce_chart = band_chart(t, cmp.included, "ce", f, unit, percentile=(5, 95))
        if ce_chart is not None:
            st.altair_chart(ce_chart, use_container_width=True)
    else:
        df = pd.DataFrame({"t (min)": t})
        any_ce = False
        for r in cmp.included:
            if float(np.max(r.ce)) > 0:  # only models with a ke0 link have an effect-site curve
                df[r.model_id.split(".")[-1] + f" (tier {r.tier})"] = r.ce * f
                any_ce = True
        if any_ce:
            st.line_chart(df, x="t (min)")
    st.subheader(f"Plasma concentration ({unit})")
    if show_bands:
        cp_chart = band_chart(t, cmp.included, "cp", f, unit, percentile=(5, 95))
        if cp_chart is not None:
            st.altair_chart(cp_chart, use_container_width=True)
    else:
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

    # v0.2 uncertainty-aware divergence: does the gap survive the models' own BSV?
    if cmp.bands:
        st.subheader("Are the models distinguishable? (v0.2)")
        sep = d.get("separation")
        if sep and sep.get("value") is not None:
            verdict = "bands DISJOINT" if sep["bands_disjoint_at_tstar"] else "bands overlap"
            st.metric(f"Separation index @ t*  ({verdict})", f"{sep['value']:+.2f}",
                      f"{100 * sep['fraction_trajectory_disjoint']:.0f}% of trajectory disjoint")
            st.caption(
                f"Driver pair **{sep['driver_high'].split('.')[-1]}** vs "
                f"**{sep['driver_low'].split('.')[-1]}**, band-tier {sep['band_tier']}. "
                "separation > 0 ⇒ a structural disagreement neither model's stated "
                "variability explains away — *model-selection risk you cannot variability-away.*")
        vs = d.get("variance_share")
        if vs:
            st.caption(f"**Variance decomposition @ t* = {vs['t_star_min']:.1f} min** — "
                       "what dominates the predictive uncertainty here?")
            st.bar_chart(
                pd.DataFrame({"share": [vs["structural"], vs["bsv"], vs["residual"]]},
                             index=["structural (which model)", "BSV (which patient)",
                                    "residual (Σ)"]),
                horizontal=True)
        for e in cmp.excluded_from_bands:
            st.warning(f"⬚ {e['model_id'].split('.')[-1]} — excluded from band math: {e['reason']}")
        if not any(getattr(r, "ce_quantiles", None) is not None
                   or getattr(r, "cp_quantiles", None) is not None for r in cmp.included):
            st.info("No band-eligible model for this drug yet — only models that publish "
                    "between-subject variability (Ω) earn a band. Today that is Eleveld propofol.")

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

# --- PD effect (BIS) prediction band: PK BSV propagated through the Hill link (v0.2 §14) ---
# Compose the band-eligible PK model with a matching PD model: the same curated Ω that
# scatters the concentration scatters the predicted effect. Quantiles are taken on the
# effect draws directly (correct under the monotone Hill transform). PD-parameter BSV
# (Ce50, γ) is not curated, so the ribbon is an honest LOWER BOUND on true effect spread.
if show_bands:
    pk_band = next((r for r in cmp.included
                    if ds[r.model_id].has_published_variability), None)
    pds = [m.id for m in ds if m.purpose == "pd" and f".{drug}." in m.id]
    if pk_band is not None and pds:
        author = pk_band.model_id.split(".")[-1].split("_")[0]
        pd_id = next((p for p in pds if author in p), pds[0])
        eff = hypnos.simulate(ds, pk_band.model_id, patient=patient, schedule=schedule, t=t,
                              pd_model=pd_id, bands=True, percentile=(5, 95),
                              samples=int(samples), seed=int(seed))
        if eff.effect_quantiles is not None:
            st.divider()
            st.subheader(f"PD effect band — {eff.effect_label}")
            st.caption(
                f"PK between-subject variability from **{pk_band.model_id.split('.')[-1]}** "
                f"pushed through the (fixed) Hill link. PD-parameter BSV (Ce50, γ) is not "
                "curated, so this 5–95% ribbon is an honest **lower bound** on true effect "
                "spread (never-invent, carried into effect space). NOT a per-patient range.")
            eq = eff.effect_quantiles
            band_df = pd.DataFrame({"t (min)": t, "lo": eq[5], "hi": eq[95]})
            med_df = pd.DataFrame({"t (min)": t, "BIS": eq[50]})
            ribbon = alt.Chart(band_df).mark_area(opacity=0.20, color="#9467bd").encode(
                x="t (min):Q", y=alt.Y("lo:Q", title="effect (BIS — lower = deeper)"), y2="hi:Q")
            line = alt.Chart(med_df).mark_line(color="#6a3d9a").encode(
                x="t (min):Q", y="BIS:Q")
            st.altair_chart(alt.layer(ribbon, line).properties(height=300),
                            use_container_width=True)
            j = int(np.argmin(eq[50]))   # peak effect = minimum median BIS
            st.caption(f"Peak-effect BIS band: **{eq[50][j]:.1f}** "
                       f"[{eq[5][j]:.1f}, {eq[95][j]:.1f}]  (band-tier {eff.band_tier}, "
                       f"PK-BSV only — lower bound).")

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
