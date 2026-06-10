#!/usr/bin/env python3
"""Deterministically regenerate all generated artifacts: exports + figures.

Exports are pure projections of the dataset (no plotting deps). Figures require
matplotlib (an optional dev dependency) and are skipped with a notice if it is
not installed, so this script always runs the export half in CI.

    python scripts/regenerate.py            # exports + figures (if matplotlib)
    python scripts/regenerate.py --exports  # exports only

The dataset is the single source of truth; everything this script writes is a
deterministic function of it (and of the dataset version). Re-running on the
same dataset produces byte-identical exports, and every README figure is
regenerated from the live dataset/kernels — so a figure can never silently rot
out of sync with the data (e.g. when a kernel is implemented or a model added).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import hypnos
from hypnos.export import FORMATS, bibtex, combine, csv_flat, export_model
from hypnos.reference import (
    alveolar_washin,
    alveolar_washout,
    greco_response_surface,
    mac_age_corrected,
)

ROOT = Path(__file__).resolve().parent.parent
EXPORTS = ROOT / "exports"
# Figures ARE deterministic projections of the dataset, exactly like the exports,
# so they regenerate in place into the committed docs/images/ (the README's
# figures) — driven by the live compare()/simulate()/kernels, never hand-drawn.
IMAGES = ROOT / "docs" / "images"

# Stable per-model colours so a model keeps its colour across figures.
_COLORS = {
    "eleveld_2018": "#2ca02c", "marsh_1991": "#1f77b4", "schnider_1998": "#d62728",
    "paedfusor_2005": "#9467bd", "kataria_1994": "#ff7f0e",
}
_DISCLAIMER = "NOT FOR CLINICAL USE"


def _params(model) -> dict:
    return {p.symbol: p.central for p in model.parameters}


def regenerate_exports(ds) -> int:
    EXPORTS.mkdir(parents=True, exist_ok=True)
    pk = [m for m in ds if m.purpose == "pk"]
    count = 0
    for fmt in [f for f in FORMATS if f not in ("omex", "bibtex", "csv")]:
        d = EXPORTS / fmt
        d.mkdir(parents=True, exist_ok=True)
        for m in pk:
            fname, text = export_model(fmt, m, ds)
            (d / fname).write_text(text, encoding="utf-8")
            count += 1
    # dataset-level text exports + per-model omex + a single dataset omex
    (EXPORTS / "citations.bib").write_text(bibtex.build(ds), encoding="utf-8")
    (EXPORTS / "parameters.csv").write_text(csv_flat.build(ds), encoding="utf-8")
    omex_dir = EXPORTS / "omex"
    omex_dir.mkdir(parents=True, exist_ok=True)
    for m in pk:
        (omex_dir / combine.model_filename(m)).write_bytes(combine.build_model_archive(m, ds))
    (EXPORTS / "hypnos.omex").write_bytes(combine.build_dataset_archive(ds))
    print(f"exports: {count} text file(s) + {len(pk)} omex + dataset omex + citations.bib + parameters.csv -> {EXPORTS}/")
    return count


def _fig_divergence(ds, plt) -> None:
    """The headline: propofol effect-site, every eligible model, envelope-greyed."""
    t = np.linspace(0, 60, 361)
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for ax, patient, title in [
        (axes[0], dict(age=72, weight=60, height=162, sex="F"), "Elderly patient (72 y, 60 kg, 162 cm, F)"),
        (axes[1], dict(age=40, weight=140, height=172, sex="M"), "Obese patient (40 y, 140 kg, 172 cm, M)"),
    ]:
        cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=sched, t=t)
        # effect-site panel: only models that actually carry a ke0 link (a PK-only
        # model has ce ≡ 0 and would draw a misleading flat line at zero).
        for r in sorted(cmp.included, key=lambda r: r.model_id):
            if float(np.max(r.ce)) <= 0:
                continue
            name = r.model_id.split(".")[-1]
            ax.plot(t, r.ce, color=_COLORS.get(name, "#333"), lw=2.3, label=f"{name}  (tier {r.tier})")
        for e in cmp.excluded:
            if float(np.max(e["result"].ce)) <= 0:
                continue
            name = e["model_id"].split(".")[-1]
            reason = "out of envelope" if any("ENVELOPE" in x for x in e["reasons"]) else "greyed"
            ax.plot(t, e["result"].ce, color="#aaa", lw=2, ls="--",
                    label=f"{name}  (tier {e['tier']} — greyed: {reason})")
        d = cmp.divergence.get("ce", {})
        sub = (f"peak divergence {d['max_abs']:.1f} µg/mL ({100*d['max_rel']:.0f}%)"
               if d else "single model in envelope")
        ax.set_title(f"{title}\n{sub}", fontsize=11)
        ax.set_xlabel("time (min)"); ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
    axes[0].set_ylabel("propofol effect-site conc. (µg/mL)")
    fig.suptitle("Hypnos model-divergence view — propofol effect-site (Marsh · Schnider · Eleveld)"
                 f"  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "divergence.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_synergy(ds, plt) -> None:
    """Propofol-remifentanil hypnotic synergy: the time course + the Greco surface."""
    t = np.linspace(0, 30, 181)
    patient = dict(age=40, weight=70, height=170, sex="M")
    prop = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    remi = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]
    alone = hypnos.simulate(ds, "hypnotics_iv.propofol.schnider_1998", patient=patient,
                            schedule=prop, t=t, pd_model="pd_effect.propofol.bis_sigmoid")
    ir = hypnos.simulate_interaction(ds, "interactions.propofol_remifentanil.greco_bis",
                                     pk_a="hypnotics_iv.propofol.schnider_1998",
                                     pk_b="opioids.remifentanil.minto_1997", patient=patient,
                                     schedule_a=prop, schedule_b=remi, t=t)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    axL.plot(t, alone.effect, color="#1f77b4", lw=2.3, label=f"propofol alone (tier {alone.tier})")
    axL.plot(t, ir.effect, color="#d62728", lw=2.3, label=f"propofol + remifentanil (tier {ir.tier})")
    axL.axhline(50, color="#888", ls=":", lw=1); axL.text(t[-1], 51, "BIS 50", ha="right", color="#888")
    axL.set_title("Hypnotic synergy — same propofol dose, ± remifentanil", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("BIS index"); axL.set_ylim(0, 95)
    axL.legend(fontsize=9); axL.grid(alpha=0.25)

    sp = _params(ds["interactions.propofol_remifentanil.greco_bis"])
    cp = np.linspace(0, 10, 160)          # propofol Ce, ug/mL
    cr_ng = np.linspace(0, 40, 160)       # remifentanil Ce, ng/mL
    CP, CR = np.meshgrid(cp, cr_ng)
    bis = greco_response_surface(CP, CR / 1000.0, E0=sp["E0"], Emax=sp["Emax"],
                                 Ce50_a=sp["Ce50_prop"], Ce50_b=sp["Ce50_remi"],
                                 alpha=sp["alpha"], gamma=sp["gamma"])
    cf = axR.contourf(CP, CR, bis, levels=12, cmap="viridis_r")
    cs = axR.contour(CP, CR, bis, levels=[40, 50, 60], colors="white", linewidths=1.3)
    axR.clabel(cs, fmt=lambda v: f"BIS {v:.0f}", fontsize=8)
    fig.colorbar(cf, ax=axR, label="BIS index")
    axR.set_title("Greco response surface (illustrative, tier C)", fontsize=11)
    axR.set_xlabel("propofol Ce (µg/mL)"); axR.set_ylabel("remifentanil Ce (ng/mL)")
    fig.suptitle(f"Hypnos Phase B — propofol–remifentanil interaction  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "synergy.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_pediatric(ds, plt) -> None:
    """Pediatric propofol: who is in-envelope for a child, who is a Tier-D extrapolation."""
    t = np.linspace(0, 30, 181)
    patient = dict(age=6, weight=20, height=115, sex="M")
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=sched, t=t)
    fig, ax = plt.subplots(figsize=(11, 5.4))
    for r in sorted(cmp.included, key=lambda r: r.model_id):
        name = r.model_id.split(".")[-1]
        ax.plot(t, r.cp, color=_COLORS.get(name, "#2ca02c"), lw=2.6,
                label=f"{name}  (tier {r.tier}, IN envelope)")
    for e in cmp.excluded:
        name = e["model_id"].split(".")[-1]
        ped = any("PEDIATRIC" in x for x in e["reasons"])
        tag = "pediatric extrapolation" if ped else "out of envelope"
        ax.plot(t, e["result"].cp, lw=2, ls="--",
                label=f"{name}  (tier {e['tier']} — greyed: {tag})")
    ax.set_title("Hypnos Phase C — pediatric propofol (6 y, 20 kg child)\n"
                 "The pediatric pair Kataria & Paedfusor + the broad-envelope Eleveld cover the "
                 "child; the adult-only models are explicit Tier-D extrapolations", fontsize=11)
    ax.set_xlabel("time (min)"); ax.set_ylabel("propofol plasma conc. (µg/mL)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    fig.text(0.5, -0.02, f"{_DISCLAIMER} — research / education / simulation only",
             ha="center", color="#c0392b", fontsize=9)
    fig.tight_layout(); fig.savefig(IMAGES / "pediatric.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_mac_age(ds, plt) -> None:
    """Volatile MAC age-correction: absolute and the one universal normalized curve."""
    ages = np.linspace(1, 90, 200)
    agents = [("desflurane", "#d62728"), ("sevoflurane", "#1f77b4"), ("isoflurane", "#2ca02c")]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    for name, color in agents:
        mac40 = _params(ds[f"volatiles.{name}.mac"])["MAC40"]
        mac = np.array([mac_age_corrected(mac40, a) for a in ages])
        axL.plot(ages, mac, color=color, lw=2.3, label=f"{name} (MAC40 {mac40:g}%)")
        axR.plot(ages, mac / mac40, color=color, lw=2.3, label=name)
    for ax in (axL, axR):
        ax.axvline(40, color="#999", ls=":", lw=1); ax.set_xlabel("age (years)"); ax.grid(alpha=0.25)
    axL.text(41, axL.get_ylim()[1] * 0.92, "age 40\n(anchor)", color="#777", fontsize=8)
    axL.set_ylabel("MAC (vol%)"); axL.set_title("Age-corrected MAC (Mapleson/Nickalls 2003)", fontsize=11)
    axL.legend(fontsize=9)
    axR.axhline(1.0, color="#999", ls=":", lw=1); axR.set_ylabel("MAC / MAC40")
    axR.set_title("Normalized MAC/MAC40 — one universal curve\n"
                  "MAC(age) = MAC40 · 10^(−0.00269·(age−40))  (~6%/decade)", fontsize=11)
    axR.legend(fontsize=9, title="all agents overlap")
    fig.suptitle(f"Hypnos Phase D — volatile-agent MAC  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "mac_age.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_washin(ds, plt) -> None:
    """Inhalational wash-in (FA/FI): solubility sets the speed, straight from blood:gas."""
    t = np.linspace(0, 10, 240)
    agents = [("desflurane", "#d62728"), ("nitrous_oxide", "#9467bd"),
              ("sevoflurane", "#1f77b4"), ("isoflurane", "#2ca02c")]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    lam_pts, plat_pts = [], []
    for name, color in agents:
        lam = _params(ds[f"volatiles.{name}.mac"])["blood_gas"]
        fa_fi, plateau, tau = alveolar_washin(lam, t)
        axL.plot(t, fa_fi, color=color, lw=2.3, label=f"{name} (λ {lam:g})")
        axL.axhline(plateau, color=color, ls=":", lw=1, alpha=0.6)
        lam_pts.append(lam); plat_pts.append(plateau)
        axR.scatter([lam], [plateau], color=color, s=70, zorder=3)
        axR.annotate(name, (lam, plateau), textcoords="offset points", xytext=(8, 4), fontsize=9, color=color)
    axL.set_title("Alveolar wash-in FA/FI(t) — single-compartment, standard ventilation", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("FA / FI"); axL.set_ylim(0, 0.8)
    axL.legend(fontsize=9, title="lower λ = faster"); axL.grid(alpha=0.25)
    order = np.argsort(lam_pts)
    axR.plot(np.array(lam_pts)[order], np.array(plat_pts)[order], color="#888", lw=1, ls="--", zorder=1)
    axR.set_title("Solubility sets the wash-in 'knee'\nhigher blood:gas λ → lower early FA/FI plateau", fontsize=11)
    axR.set_xlabel("blood:gas partition coefficient λ"); axR.set_ylabel("early FA/FI plateau")
    axR.grid(alpha=0.25)
    fig.suptitle(f"Hypnos Phase D — inhalational wash-in (FA/FI uptake)  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "washin.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_washout(ds, plt) -> None:
    """Inhalational wash-out (FA/FA₀): the offset mirror — solubility sets emergence speed."""
    t = np.linspace(0, 10, 240)
    agents = [("desflurane", "#d62728"), ("nitrous_oxide", "#9467bd"),
              ("sevoflurane", "#1f77b4"), ("isoflurane", "#2ca02c")]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    lam_pts, floor_pts = [], []
    for name, color in agents:
        lam = _params(ds[f"volatiles.{name}.mac"])["blood_gas"]
        fa, floor, tau = alveolar_washout(lam, t)
        axL.plot(t, fa, color=color, lw=2.3, label=f"{name} (λ {lam:g})")
        axL.axhline(floor, color=color, ls=":", lw=1, alpha=0.6)
        lam_pts.append(lam); floor_pts.append(floor)
        axR.scatter([lam], [floor], color=color, s=70, zorder=3)
        axR.annotate(name, (lam, floor), textcoords="offset points", xytext=(8, 4), fontsize=9, color=color)
    axL.set_title("Alveolar wash-out FA/FA₀(t) — single-compartment, standard ventilation", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("FA / FA₀"); axL.set_ylim(0, 1.0)
    axL.legend(fontsize=9, title="lower λ = faster emergence"); axL.grid(alpha=0.25)
    order = np.argsort(lam_pts)
    axR.plot(np.array(lam_pts)[order], np.array(floor_pts)[order], color="#888", lw=1, ls="--", zorder=1)
    axR.set_title("Solubility sets the wash-out floor\nhigher blood:gas λ → higher residual FA/FA₀", fontsize=11)
    axR.set_xlabel("blood:gas partition coefficient λ"); axR.set_ylabel("early elimination floor (1 − plateau)")
    axR.grid(alpha=0.25)
    fig.suptitle(f"Hypnos Phase D — inhalational wash-out (FA/FA₀ emergence)  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "washout.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_variability(ds, plt) -> None:
    """v0.2 headline: the prediction band + variance decomposition.

    Of the three propofol models eligible for the elderly patient, only Eleveld
    publishes the random-effects structure, so only it earns a band; Marsh and
    Schnider are drawn as bare median lines and named as excluded from the band
    math — the never-synthesize rule made visible (a missing band is honest, a
    borrowed one is a lie with error bars).
    """
    t = np.linspace(0, 60, 361)
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    patient = dict(age=72, weight=60, height=162, sex="F")
    cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=sched, t=t,
                         bands=True, percentile=(5, 95), samples=2000, seed=7)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    for r in sorted(cmp.included, key=lambda r: r.model_id):
        if float(np.max(r.ce)) <= 0:
            continue
        name = r.model_id.split(".")[-1]
        color = _COLORS.get(name, "#333")
        if r.ce_quantiles is not None:
            lo, hi = r.band_percentile
            axL.fill_between(t, r.ce_quantiles[lo], r.ce_quantiles[hi], color=color, alpha=0.18,
                             label=f"{name} {lo}–{hi}% band (band-tier {r.band_tier})")
            axL.plot(t, r.ce_quantiles[50], color=color, lw=2.4)
        else:
            axL.plot(t, r.ce, color=color, lw=2.0, ls="--",
                     label=f"{name} (tier {r.tier} — no published BSV; line only)")
    axL.set_title("Propofol effect-site — only Eleveld publishes between-subject\n"
                  "variability, so only Eleveld earns a band (never-synthesize rule)", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("propofol Ce (µg/mL)")
    axL.legend(fontsize=8.5); axL.grid(alpha=0.25)

    # right: time-resolved variance decomposition for the band-eligible model
    elev = next(r for r in cmp.included if r.ce_quantiles is not None)
    bsv = elev.ce_bsv_var
    res = elev.ce_resid_var
    total = bsv + res
    with np.errstate(divide="ignore", invalid="ignore"):
        f_bsv = np.where(total > 0, bsv / total, 0.0)
        f_res = np.where(total > 0, res / total, 0.0)
    axR.stackplot(t, f_bsv, f_res, labels=["between-subject (η / Ω)", "residual (ε / Σ)"],
                  colors=["#2ca02c", "#bbbbbb"], alpha=0.85)
    vs = cmp.divergence["ce"]["variance_share"]
    axR.axvline(vs["t_star_min"], color="#333", ls=":", lw=1)
    axR.set_ylim(0, 1); axR.set_xlim(t.min(), t.max())
    axR.set_title("Where does the uncertainty come from? (Eleveld)\n"
                  f"at t*={vs['t_star_min']:g} min: BSV {vs['bsv']:.0%} vs residual {vs['residual']:.0%}"
                  "  — the patient, not the assay", fontsize=11)
    axR.set_xlabel("time (min)"); axR.set_ylabel("share of within-model predictive variance")
    axR.legend(fontsize=9, loc="lower right")
    fig.suptitle("Hypnos v0.2 — population-variability layer (seeded, reproducible bands)"
                 f"  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "variability.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_effect_band(ds, plt) -> None:
    """The PD effect band: PK between-subject variability propagated through the link.

    The same curated Eleveld Ω that scatters the effect-site concentration also
    scatters the predicted BIS — each virtual individual's true Ce, pushed through
    the (deterministic) two-slope Hill model. Because PD-parameter BSV (Ce50, γ) is
    not curated, the effect ribbon is an honest LOWER BOUND on true effect spread.
    """
    t = np.linspace(0, 30, 301)
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    patient = dict(age=72, weight=60, height=162, sex="F")
    r = hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018", patient=patient,
                        schedule=sched, t=t, pd_model="pd_effect.propofol.eleveld_bis",
                        bands=True, percentile=(5, 95), samples=2000, seed=7)
    lo, hi = r.band_percentile
    color = _COLORS.get("eleveld_2018", "#1f77b4")

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))
    axL.fill_between(t, r.ce_quantiles[lo], r.ce_quantiles[hi], color=color, alpha=0.18,
                     label=f"{lo}–{hi}% band (band-tier {r.band_tier})")
    axL.plot(t, r.ce_quantiles[50], color=color, lw=2.4, label="population median")
    axL.plot(t, r.ce, color="#333", lw=1.2, ls="--", label="typical individual")
    axL.set_title("Propofol effect-site concentration (Eleveld, with BSV)", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("propofol Ce (µg/mL)")
    axL.legend(fontsize=8.5); axL.grid(alpha=0.25)

    axR.fill_between(t, r.effect_quantiles[lo], r.effect_quantiles[hi], color="#9467bd", alpha=0.20,
                     label=f"{lo}–{hi}% effect band (lower bound)")
    axR.plot(t, r.effect_quantiles[50], color="#6a3d9a", lw=2.4, label="population median BIS")
    axR.plot(t, r.effect, color="#333", lw=1.2, ls="--", label="typical individual")
    axR.set_ylim(0, 100); axR.invert_yaxis()  # lower BIS (deeper hypnosis) plotted downward
    axR.set_title("Predicted BIS — the SAME PK variability, pushed through the\n"
                  "PD link; PD-parameter BSV uncurated ⇒ a lower bound", fontsize=11)
    axR.set_xlabel("time (min)"); axR.set_ylabel("BIS  (lower = deeper hypnosis)")
    axR.legend(fontsize=8.5, loc="lower right"); axR.grid(alpha=0.25)

    fig.suptitle("Hypnos v0.2 — PD effect band (PK BSV propagated through the Hill link)"
                 f"  ·  {_DISCLAIMER}", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "effect_band.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_la_double_uncertainty(ds, plt) -> None:
    """v0.6 LA1 — the double-uncertainty view: predicted concentration vs the
    threshold *bands*, with the honest punchline that the threshold uncertainty
    dwarfs the PK uncertainty (so no single safe-concentration line is defensible).

    LEFT: the same bupivacaine dose at every site (the site-dominance trajectory),
    with the CNS and cardiovascular threshold RANGES drawn as shaded bands — never
    lines. RIGHT: the total vs free trace at one site, showing the binding-saturation
    gap. Framed throughout as research/education, never a dosing margin (v0.6 §6/§7).
    """
    from hypnos.la import concentration_at_site, double_uncertainty, free_concentration

    mid = "local_anesthetics.bupivacaine.systemic"
    dose, horizon = 100.0, 90.0
    du = double_uncertainty(ds, mid, site="lumbar_epidural", dose_mg=dose, t_min=horizon)
    sites = ["intercostal", "lumbar_epidural", "brachial_plexus", "subcutaneous"]
    site_colors = {"intercostal": "#d62728", "lumbar_epidural": "#1f77b4",
                   "brachial_plexus": "#2ca02c", "subcutaneous": "#888"}

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))

    # LEFT — total-plasma trajectory by site + total-basis threshold BANDS
    for s in sites:
        r = concentration_at_site(ds, mid, site=s, dose_mg=dose, t_min=horizon)
        axL.plot(r.t, r.cp, color=site_colors[s], lw=2.2, label=f"{s} (Cmax {r.cmax:.2f})")
    band_styles = {"cns_first_symptoms": ("#f0ad4e", "CNS first symptoms"),
                   "cardiovascular": ("#c0392b", "cardiovascular")}
    for e in du.endpoints:
        if e.basis != "total_plasma":
            continue
        color, lab = band_styles.get(e.endpoint, ("#999", e.endpoint))
        axL.axhspan(e.low, e.high, color=color, alpha=0.18,
                    label=f"{lab} threshold RANGE [{e.low:g}–{e.high:g}]")
    axL.set_title(f"Bupivacaine {dose:g} mg — total plasma by site\n"
                  "vs toxicity-threshold RANGES (never lines)", fontsize=11)
    axL.set_xlabel("time (min)"); axL.set_ylabel("total plasma conc. (µg/mL)")
    axL.legend(fontsize=8); axL.grid(alpha=0.25)

    # RIGHT — total vs free trace at one site + the free-basis threshold band
    r = concentration_at_site(ds, mid, site="lumbar_epidural", dose_mg=dose, t_min=horizon)
    pb = ds.drug("bupivacaine")["protein_binding"]
    fc = free_concentration(r.cp, pb)
    axR.plot(r.t, r.cp, color="#1f77b4", lw=2.3, label="total plasma")
    axR.plot(r.t, fc.c_free, color="#6a3d9a", lw=2.3,
             label=f"free (linear, ff={fc.free_fraction:g}) — saturable ⇒ a lower bound")
    for e in du.endpoints:
        if e.basis == "free_plasma":
            axR.axhspan(e.low, e.high, color="#6a3d9a", alpha=0.16,
                        label=f"free CNS threshold RANGE [{e.low:g}–{e.high:g}]")
    axR.set_title("Total vs FREE concentration (lumbar epidural)\n"
                  "toxicity tracks free drug; binding saturation widens the gap", fontsize=11)
    axR.set_xlabel("time (min)"); axR.set_ylabel("conc. (µg/mL)")
    axR.legend(fontsize=8); axR.grid(alpha=0.25)

    fig.suptitle("Hypnos v0.6 LA1 — the double-uncertainty view: threshold uncertainty dwarfs "
                 f"the PK spread ⇒ no single safe line  ·  {_DISCLAIMER}", y=1.02, fontsize=10.5)
    fig.text(0.5, -0.03, "RESEARCH / EDUCATION ONLY — not a dosing tool, no maximum dose, no "
             "margin-as-guarantee, no 'is this safe?' answer (v0.6 §7)",
             ha="center", color="#c0392b", fontsize=9)
    fig.tight_layout(); fig.savefig(IMAGES / "la_double_uncertainty.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_la_cardiotoxicity(ds, plt) -> None:
    """v0.6 LA2 — the agent-choice cardiotoxicity comparison: why a similar CNS
    threshold can hide a very different cardiovascular margin.

    For each amide, the CNS-first-symptoms range and the cardiovascular range are
    drawn as bands on a shared (log) total-plasma axis; the GAP between them is the
    margin. It widens monotonically bupivacaine -> levobupivacaine -> ropivacaine ->
    lidocaine, mirroring the curated cardiotoxicity class. Comparative / educational,
    never a dose recommendation (v0.6 §7).
    """
    from hypnos.la import cardiotoxicity_comparison

    rows = cardiotoxicity_comparison(ds)
    fig, ax = plt.subplots(figsize=(11, 4.8))
    cns_c, cvs_c = "#f0ad4e", "#c0392b"
    for i, a in enumerate(rows):
        m = ds[a.model_id]
        cns = next((t for t in m.toxicity_thresholds
                    if t.endpoint == "cns_first_symptoms" and t.basis == "total_plasma"), None)
        cvs = next((t for t in m.toxicity_thresholds if t.endpoint == "cardiovascular"), None)
        if cns is None or cvs is None:
            continue
        y = len(rows) - 1 - i
        ax.barh(y, cns.high - cns.low, left=cns.low, height=0.32, color=cns_c,
                label="CNS first symptoms" if i == 0 else None)
        ax.barh(y, cvs.high - cvs.low, left=cvs.low, height=0.32, color=cvs_c,
                label="cardiovascular" if i == 0 else None)
        fold = f"{a.cns_to_cvs_fold:.1f}x" if a.cns_to_cvs_fold else "?"
        ax.text(cvs.high * 1.05, y, f"margin {a.cns_to_cvs_margin} ({fold})  ·  {a.rank}",
                va="center", fontsize=9, color="#333")
        ax.text(cns.low * 0.95, y, a.drug, va="center", ha="right", fontsize=9.5, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(0.8, 60)
    ax.set_yticks([])
    ax.set_xlabel("total plasma concentration (µg/mL, log scale) — toxicity-threshold RANGES, never lines")
    ax.set_title("Hypnos v0.6 LA2 — stereochemistry & the agent-choice cardiotoxicity margin\n"
                 "the CNS→cardiovascular gap widens bupivacaine → levobupivacaine → ropivacaine → lidocaine "
                 f"·  {_DISCLAIMER}", fontsize=10.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), fontsize=9, frameon=True)
    ax.grid(alpha=0.25, axis="x")
    fig.text(0.5, -0.02, "RESEARCH / EDUCATION ONLY — comparative agent differences, not a dose "
             "recommendation; thresholds are wide, individual, Tier-C (v0.6 §7)",
             ha="center", color="#c0392b", fontsize=9)
    fig.tight_layout(); fig.savefig(IMAGES / "la_cardiotoxicity.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fig_la_saturation(ds, plt) -> None:
    """v0.6 LA3 — the binding-saturation failure mode made visible: total
    concentration under-predicts free-drug (toxic) concentration exactly when risk
    is highest, because plasma-protein binding saturates and the free fraction rises.

    LEFT: the free FRACTION vs total plasma concentration — the curated capacity-
    limited (Langmuir, Tier-D illustrative) curves climb from the low-concentration
    fraction toward 1 as binding saturates, while the linear assumption (dashed) stays
    flat. RIGHT: bupivacaine's free CONCENTRATION, linear vs non-linear, with the
    under-prediction gap shaded. Comparative / educational, never a dose (v0.6 §7).
    """
    from hypnos.la import saturable_free_concentration

    agents = [("bupivacaine", "#c0392b"), ("levobupivacaine", "#8e44ad"), ("ropivacaine", "#2980b9")]
    ct = np.logspace(np.log10(0.2), np.log10(40.0), 240)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.8))

    for name, color in agents:
        pb = ds.drug(name)["protein_binding"]
        fb = pb["fraction_bound"]
        cap = pb["free_fraction_model"]["binding_capacity_ug_ml"]
        free_nl = saturable_free_concentration(ct, fb, cap)
        axL.plot(ct, 100.0 * free_nl / ct, color=color, lw=2.3,
                 label=f"{name} (non-linear, cap {cap:g} µg/mL)")
        axL.axhline(100.0 * (1.0 - fb), color=color, lw=1.1, ls=":", alpha=0.7)
    axL.set_xscale("log")
    axL.set_xlabel("total plasma concentration (µg/mL, log scale)")
    axL.set_ylabel("free fraction (%)")
    axL.set_title("Free fraction RISES as binding saturates\n(dotted = the flat linear assumption)", fontsize=11)
    axL.legend(fontsize=8.5, loc="upper left"); axL.grid(alpha=0.25)

    pb = ds.drug("bupivacaine")["protein_binding"]
    fb, cap = pb["fraction_bound"], pb["free_fraction_model"]["binding_capacity_ug_ml"]
    free_nl = saturable_free_concentration(ct, fb, cap)
    free_lin = ct * (1.0 - fb)
    axR.fill_between(ct, free_lin, free_nl, color="#c0392b", alpha=0.15,
                     label="under-prediction gap")
    axR.plot(ct, free_nl, color="#c0392b", lw=2.4, label="free (non-linear, Tier-D)")
    axR.plot(ct, free_lin, color="#333", lw=1.6, ls="--", label="free (linear assumption)")
    axR.set_xlabel("total plasma concentration (µg/mL)")
    axR.set_ylabel("free (unbound) concentration (µg/mL)")
    axR.set_xlim(0, 20); axR.set_ylim(0, 8)
    axR.set_title("Bupivacaine: the linear free trace UNDER-predicts\nfree-drug risk at high total", fontsize=11)
    axR.legend(fontsize=8.5, loc="upper left"); axR.grid(alpha=0.25)

    fig.suptitle("Hypnos v0.6 LA3 — the binding-saturation failure mode (total under-predicts free "
                 f"when risk is highest)  ·  {_DISCLAIMER}", y=1.02, fontsize=10.5)
    fig.text(0.5, -0.03, "RESEARCH / EDUCATION ONLY — capacity-limited model is Tier-D illustrative; "
             "the QUALITATIVE rise is the documented fact, not the magnitude (v0.6 §7)",
             ha="center", color="#c0392b", fontsize=9)
    fig.tight_layout(); fig.savefig(IMAGES / "la_saturation.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def regenerate_figures(ds) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("figures: matplotlib not installed; skipped (pip install matplotlib)")
        return False
    IMAGES.mkdir(parents=True, exist_ok=True)
    _fig_divergence(ds, plt)
    _fig_synergy(ds, plt)
    _fig_pediatric(ds, plt)
    _fig_mac_age(ds, plt)
    _fig_washin(ds, plt)
    _fig_washout(ds, plt)
    _fig_variability(ds, plt)
    _fig_effect_band(ds, plt)
    _fig_la_double_uncertainty(ds, plt)
    _fig_la_cardiotoxicity(ds, plt)
    _fig_la_saturation(ds, plt)
    print(f"figures: regenerated divergence, synergy, pediatric, mac_age, washin, washout, "
          f"variability, effect_band, la_double_uncertainty, la_cardiotoxicity, la_saturation -> {IMAGES}/")
    return True


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ds = hypnos.load()
    problems = hypnos.validate_dataset(ds)
    if problems:
        print(f"refusing to regenerate: dataset has {len(problems)} validation problem(s)", file=sys.stderr)
        return 1
    regenerate_exports(ds)
    if "--exports" not in argv:
        regenerate_figures(ds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
