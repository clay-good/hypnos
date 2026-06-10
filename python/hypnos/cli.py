"""Hypnos command-line interface.

    hypnos version
    hypnos validate
    hypnos info
    hypnos simulate <model_id> --age .. --weight .. --height .. --sex .. [--pd <id>]
    hypnos compare  --drug propofol --age .. --weight .. ...
    hypnos export --format <fmt> --output <dir> [--model <id>]

Forward simulation only. No inverse control, no dosing recommendation (spec §10).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

from . import CLINICAL_USE, __version__
from .analysis import decrement_time, time_to_peak_effect
from .export import FORMATS, bibtex, combine, csv_flat, export_model
from .filter import performance_table, pk_drugs, select, summary
from .inhalational import mac as mac_eval
from .inhalational import washin, washin_comparison, washout, washout_comparison
from .load import load
from .presets import default_schedule_for
from .simulate import compare as compare_models
from .simulate import simulate, simulate_interaction
from .validate import validate_dataset
from .verification import checklist_markdown, model_verification, verification_summary

# Drug-appropriate default dose schedules live in hypnos.presets (shared with the
# dashboard so the two never drift). Aliased here for the existing call sites.
_default_schedule_for = default_schedule_for


def _patient_from_args(args) -> dict:
    patient = {"age": args.age, "weight": args.weight, "height": args.height, "sex": args.sex}
    # Organ-function covariates (v0.5 §B) — present only when supplied, so a normal
    # simulation is unaffected. Each one a model has no standing in greys it to Tier D.
    for key, attr in (("child_pugh", "child_pugh"), ("crcl_ml_min", "crcl"),
                      ("albumin_g_dl", "albumin"), ("ejection_fraction_pct", "ejection_fraction")):
        val = getattr(args, attr, None)
        if val is not None:
            patient[key] = val
    return patient


def _add_patient_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--age", type=float, default=50.0)
    p.add_argument("--weight", type=float, default=77.0)
    p.add_argument("--height", type=float, default=177.0)
    p.add_argument("--sex", default="M")
    # Organ-function covariates (optional; v0.5 organ-failure envelope).
    p.add_argument("--child-pugh", dest="child_pugh", choices=["A", "B", "C"], default=None,
                   help="Child-Pugh class (hepatic impairment)")
    p.add_argument("--crcl", type=float, default=None,
                   help="creatinine clearance mL/min (renal impairment if <60)")
    p.add_argument("--albumin", type=float, default=None,
                   help="serum albumin g/dL (hypoalbuminemia if <3.5)")
    p.add_argument("--ejection-fraction", dest="ejection_fraction", type=float, default=None,
                   help="ejection fraction %% (low cardiac output if <40)")
    p.add_argument("--tmax", type=float, default=60.0, help="simulation horizon (min)")
    p.add_argument("--n", type=int, default=361, help="number of time samples")


def cmd_version(args) -> int:
    print(f"hypnos {__version__}")
    print(f"clinicalUse = {CLINICAL_USE}")
    return 0


def cmd_validate(args) -> int:
    problems = validate_dataset()
    if not problems:
        ds = load()
        print(f"OK — {len(ds)} model(s) valid against schema + integrity checks.")
        return 0
    print(f"FAILED — {len(problems)} problem(s):", file=sys.stderr)
    for p in problems:
        print("  " + p, file=sys.stderr)
    return 1


def cmd_info(args) -> int:
    print(json.dumps(summary(load()), indent=2))
    return 0


def _run_sim(args):
    ds = load()
    patient = _patient_from_args(args)
    t = np.linspace(0.0, args.tmax, args.n)
    schedule = _default_schedule_for(ds[args.model].drug_name)
    bands = bool(getattr(args, "bands", False))
    res = simulate(ds, args.model, patient=patient, schedule=schedule, t=t, pd_model=args.pd,
                   bands=bands, percentile=_parse_percentile(getattr(args, "percentile", "5,95")),
                   samples=getattr(args, "samples", 2000),
                   seed=getattr(args, "seed", 7) if bands else None)
    return res


def cmd_simulate(args) -> int:
    res = _run_sim(args)
    print(f"model: {res.model_id}")
    print(f"tier (propagated): {res.tier}")
    u = res.concentration_unit
    cf = res.conc_factor
    print(f"Cp peak: {res.cp_peak_display:.3f} {u}   Ce peak: {res.ce_peak_display:.3f} {u}")
    if res.ce_quantiles is not None:
        lo, hi = res.band_percentile
        q = res.ce_quantiles if res.ce_peak > 0 else res.cp_quantiles
        i = int(np.argmax(q[50]))
        print(f"  band-tier {res.band_tier}  Ce peak {q[50][i] * cf:.3f} "
              f"[{q[lo][i] * cf:.3f}, {q[hi][i] * cf:.3f}] {u}  ({lo}-{hi}%, seeded)")
    if res.effect is not None:
        print(f"effect ({res.effect_label}): min {float(np.min(res.effect)):.1f}")
        if res.effect_quantiles is not None:
            lo, hi = res.band_percentile
            j = int(np.argmin(res.effect_quantiles[50]))  # peak effect (min median BIS)
            print(f"  effect band @ peak: {res.effect_quantiles[50][j]:.1f} "
                  f"[{res.effect_quantiles[lo][j]:.1f}, {res.effect_quantiles[hi][j]:.1f}]  "
                  f"({lo}-{hi}%, PK-BSV propagated — lower bound on true effect spread)")
    if res.warnings:
        print("warnings:")
        for w in res.warnings:
            print("  - " + w)
    else:
        print("warnings: none (request within envelope)")
    return 0


def _parse_percentile(spec: str) -> tuple:
    """Parse a '5,95' percentile pair into a (lo, hi) int tuple."""
    parts = [int(x) for x in str(spec).split(",")]
    if len(parts) != 2:
        raise ValueError(f"--percentile expects 'lo,hi' (e.g. '5,95'), got {spec!r}")
    return (parts[0], parts[1])


def _mdape_str(model) -> str:
    """Compact in-envelope MDAPE (inaccuracy) badge for the divergence view."""
    vals = sorted({e["value"] for e in model.predictive_mdape})
    if not vals:
        return "MDAPE n/a"
    if len(vals) == 1:
        return f"MDAPE {vals[0]:g}%"
    return f"MDAPE {vals[0]:g}-{vals[-1]:g}%"


def cmd_compare(args) -> int:
    ds = load()
    patient = _patient_from_args(args)
    t = np.linspace(0.0, args.tmax, args.n)
    pct = _parse_percentile(args.percentile)
    cmp = compare_models(ds, drug=args.drug, patient=patient,
                         schedule=_default_schedule_for(args.drug), t=t,
                         bands=args.bands, percentile=pct, samples=args.samples,
                         seed=args.seed if args.bands else None)
    cu = cmp.concentration_unit
    cf = cmp.conc_factor
    print(f"drug: {cmp.drug}   patient: {patient}   (concentrations in {cu})")
    print(f"included ({len(cmp.included)}):")
    for r in cmp.included:
        peak = r.ce_peak_display if r.ce_peak > 0 else r.cp_peak_display
        kind = "Ce" if r.ce_peak > 0 else "Cp"
        acc = _mdape_str(ds[r.model_id])
        bt = ""
        if args.bands:
            q = r.ce_quantiles if r.ce_peak > 0 else r.cp_quantiles
            if q is not None:
                lo, hi = r.band_percentile
                i = int(np.argmax(q[50]))   # band reported at the median-peak instant
                bt = (f"  band-tier {r.band_tier}  {kind} {q[50][i]*cf:.2f} "
                      f"[{q[lo][i]*cf:.2f}, {q[hi][i]*cf:.2f}]")
            else:
                bt = "  band-tier — (no published BSV; line only)"
        print(f"  - {r.model_id:42s} tier {r.tier}  {kind} peak {peak:.3f} {cu}   {acc}{bt}")
    if cmp.excluded:
        print(f"excluded for envelope ({len(cmp.excluded)}):")
        for e in cmp.excluded:
            print(f"  - {e['model_id']:42s} tier {e['tier']}  ({'; '.join(e['reasons'])})")
    if cmp.unavailable:
        print(f"unavailable ({len(cmp.unavailable)}):")
        for u in cmp.unavailable:
            print(f"  - {u['model_id']:42s} ({u['reason']})")
    for key, label in [("cp", "plasma"), ("ce", "effect-site")]:
        d = cmp.divergence.get(key) or {}
        if not d:
            continue
        line = (f"{label} divergence across included models: "
                f"peak abs {d['max_abs'] * cf:.3f} {cu}, peak rel {100*d['max_rel']:.0f}%")
        drv = d.get("driver")
        if drv:
            line += (f"  (driver: {drv['high'].split('.')[-1]} vs {drv['low'].split('.')[-1]})")
        print(line)
        sep = d.get("separation")
        if sep is not None:
            disj = "DISJOINT" if sep["bands_disjoint_at_tstar"] else "overlapping"
            print(f"  separation@t* {sep['value']:+g} (bands {disj}); "
                  f"{100*sep['fraction_trajectory_disjoint']:.0f}% of trajectory disjoint "
                  f"[{sep['percentile'][0]}-{sep['percentile'][1]}%, band-tier {sep['band_tier']}]")
        vs = d.get("variance_share")
        if vs is not None:
            print(f"  variance share @ t*={vs['t_star_min']:g} min: "
                  f"structural {vs['structural']:.2f} | BSV {vs['bsv']:.2f} | residual {vs['residual']:.2f}")
    if args.bands and cmp.excluded_from_bands:
        for e in cmp.excluded_from_bands:
            print(f"  note: {e['model_id'].split('.')[-1]} excluded from band metrics ({e['reason']})")
    return 0


def cmd_interact(args) -> int:
    ds = load()
    patient = _patient_from_args(args)
    t = np.linspace(0.0, args.tmax, args.n)
    prop_sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
    remi_sched = [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")]
    ir = simulate_interaction(
        ds, args.surface,
        pk_a=args.propofol, pk_b=args.remifentanil,
        patient=patient, schedule_a=prop_sched, schedule_b=remi_sched, t=t,
    )
    # propofol-alone BIS via the single-drug sigmoid, for the synergy contrast
    alone = simulate(ds, args.propofol, patient=patient, schedule=prop_sched, t=t,
                     pd_model="pd_effect.propofol.bis_sigmoid")
    print(f"surface: {ir.surface_id}   patient: {patient}")
    print(f"tier (propagated): {ir.tier}")
    print(f"BIS min — propofol alone:        {float(np.min(alone.effect)):.1f}")
    print(f"BIS min — propofol + remifentanil: {ir.effect_min:.1f}  (synergy)")
    if ir.warnings:
        print("warnings:")
        for w in ir.warnings:
            print("  - " + w)
    return 0


def cmd_status(args) -> int:
    ds = load()
    s = verification_summary(ds)
    bs = s["by_review_status"]
    print(f"verification coverage: {bs.get('verified', 0)}/{s['n_models']} verified "
          f"({100 * s['verified_fraction']:.0f}%)   "
          f"unverified {bs.get('unverified', 0)}   contested {bs.get('contested', 0)}")
    print("\nstart here (highest-leverage unverified models — implemented kernel + best tier first):")
    for item in s["next_to_verify"]:
        kern = "kernel" if item["kernel"] else "no-kernel"
        print(f"  - {item['model_id']:48s} tier {item['tier']}  {kern:9s}  cite {item['citation']}")
    print("\nrun `hypnos verify <model_id>` for the field-by-field checklist.")
    return 0


def cmd_verify(args) -> int:
    ds = load()
    if args.model not in ds:
        print(f"unknown model {args.model!r}", file=sys.stderr)
        return 2
    mv = model_verification(ds, args.model)
    if args.markdown:
        sys.stdout.write(checklist_markdown(mv))
        return 0
    print(f"model: {mv.model_id}   tier {mv.tier}   status: {mv.review_status}   "
          f"kernel: {'implemented' if mv.kernel_implemented else 'pending'}")
    print(f"source: doi:{mv.doi or '?'}  PMID:{mv.pmid or '?'}"
          + (f"  ({mv.source_locator})" if mv.source_locator else ""))
    if mv.blocking:
        print("blocking before 'verified':")
        for b in mv.blocking:
            print("  - " + b)
    print(f"\nchecklist ({mv.n_items} items to confirm against the PDF):")
    n_missing = 0
    for it in mv.checklist:
        if it.locator:
            where = f"   @ {it.locator}"
        elif it.group == "structural":
            where = "   [! no source locator curated — add one]"
            n_missing += 1
        else:
            where = ""
        print(f"  [ ] ({it.group}) {it.label} = {it.value}{where}")
    if n_missing:
        print(f"\n{n_missing} structural parameter(s) have no curated source locator — "
              "add one (e.g. 'Author YYYY, Table N') while confirming against the PDF.")
    print("\nLLMs do not promote: a human confirms these, then edits review_status -> verified.")
    return 0


def cmd_models(args) -> int:
    ds = load()
    models = select(ds, drug=args.drug) if args.drug else list(ds)
    models = sorted(models, key=lambda m: m.id)
    if not models:
        print(f"no models{' for drug ' + args.drug if args.drug else ''}", file=sys.stderr)
        return 2
    print(f"{'model_id':46s} {'purpose':14s} {'tier':4s} {'kernel':8s} review")
    for m in models:
        kern = "yes" if m.kernel_implemented else "pending"
        print(f"{m.id:46s} {m.purpose:14s} {m.tier:4s} {kern:8s} {m.review_status}")
    if not args.drug:
        print(f"\nsimulatable drugs (>=1 PK kernel): {', '.join(pk_drugs(ds))}")
    return 0


def cmd_performance(args) -> int:
    ds = load()
    rows = performance_table(ds, drug=args.drug)
    if not rows:
        print(f"no published predictive-performance metrics{' for drug ' + args.drug if args.drug else ''}",
              file=sys.stderr)
        return 2
    print("Published predictive performance — the numeric counterpart to the tier "
          "(spec §5; MDPE = bias, MDAPE = inaccuracy).")
    print(f"{'model_id':46s} {'tier':4s} {'metric':10s} {'value':>8s}  population (citation)")
    for r in rows:
        val = f"{r['value']:g} {r['units']}".strip()
        src = r["citation"] or "?"
        pop = r["population"] or ""
        print(f"{r['model_id']:46s} {r['tier']:4s} {r['metric']:10s} {val:>8s}  {pop} [{src}]")
    return 0


def cmd_tpeak(args) -> int:
    ds = load()
    patient = _patient_from_args(args)
    try:
        pe = time_to_peak_effect(ds, args.model, patient=patient)
    except (ValueError, NotImplementedError) as e:
        print(str(e), file=sys.stderr)
        return 2
    print(f"model: {pe.model_id}   tier {pe.tier}")
    print(f"time to peak effect (after a bolus): {pe.tpeak_min:.2f} min   (ke0 {pe.ke0:.3f} /min)")
    print(f"Ce/Cp at peak: {pe.ce_cp_ratio_at_peak:.3f}  (1.0 by definition: dCe/dt=0 when Ce=Cp)")
    for w in pe.warnings:
        print("  - " + w)
    return 0


def cmd_decrement(args) -> int:
    ds = load()
    patient = _patient_from_args(args)
    drug = ds[args.model].drug_name
    infusion = args.infusion
    if infusion is None:  # default to the drug's preset infusion rate
        infusion = next((spec for kind, _, spec in default_schedule_for(drug) if kind == "infusion"), None)
        if infusion is None:
            print(f"no preset infusion for {drug}; pass --infusion (e.g. '6 mg/kg/h')", file=sys.stderr)
            return 2
    try:
        dt = decrement_time(ds, args.model, patient=patient, infusion=infusion,
                            duration=args.duration, fraction=args.fraction)
    except (ValueError, NotImplementedError) as e:
        print(str(e), file=sys.stderr)
        return 2
    u = (ds.drug(drug) or {}).get("concentration_unit", "ug/mL")
    from .models import concentration_factor
    cstop = dt.conc_at_stop * concentration_factor(u)
    print(f"model: {dt.model_id}   tier {dt.tier}")
    print(f"infusion: {dt.infusion} for {dt.duration_min:g} min  ->  plasma at stop {cstop:.3f} {u}")
    val = "not reached" if dt.decrement_min == float("inf") else f"{dt.decrement_min:.1f} min"
    print(f"{100*dt.fraction:g}% plasma decrement time (constant-rate; not classic CSHT): {val}")
    for w in dt.warnings:
        print("  - " + w)
    return 0


def cmd_mac(args) -> int:
    ds = load()
    agent = args.agent
    if "." not in agent:  # short name -> volatiles.<agent>.mac
        agent = f"volatiles.{agent.replace(' ', '_')}.mac"
    res = mac_eval(ds, agent, age=args.age, end_tidal_pct=args.end_tidal,
                   n2o_end_tidal_pct=args.n2o)
    print(f"agent: {res.agent_id}   age: {res.age:g} y   tier: {res.tier}")
    print(f"MAC (age-corrected): {res.mac_age:.2f} vol%   (MAC40 {res.mac40:g})")
    print(f"MAC-awake (age-corrected): {res.mac_awake_age:.2f} vol%")
    print(f"blood:gas {res.blood_gas:g}   oil:gas {res.oil_gas:g}")
    if res.mac_fraction is not None:
        print(f"end-tidal {res.end_tidal_pct:g} vol% -> MAC fraction {res.mac_fraction:.2f}")
    if res.combined_mac_fraction is not None and args.n2o is not None:
        print(f"+ N2O {args.n2o:g} vol% -> combined MAC fraction {res.combined_mac_fraction:.2f}")
    for w in res.warnings:
        print("  - " + w)
    return 0


def cmd_la(args) -> int:
    from .la import cardiotoxicity_comparison, double_uncertainty, site_comparison
    ds = load()
    if args.agents:
        rows = cardiotoxicity_comparison(ds)
        if not rows:
            print("no curated cardiotoxicity classes", file=sys.stderr)
            return 2
        print("local-anesthetic cardiotoxicity comparison (v0.6 LA2) — most cardiotoxic first")
        print("  the agent-choice axis: why a similar CNS threshold can hide a very different "
              "cardiovascular margin")
        print(f"{'drug':16s} {'class':>13s} {'stereochem':>13s} {'CNS->CVS margin':>16s} {'fold':>6s}")
        for a in rows:
            fold = f"{a.cns_to_cvs_fold:.1f}x" if a.cns_to_cvs_fold is not None else "n/a"
            print(f"{a.drug:16s} {a.rank:>13s} {a.stereochemistry:>13s} {a.cns_to_cvs_margin:>16s} {fold:>6s}")
        amides = [a for a in rows if a.rank in ("high", "intermediate")
                  or a.stereochemistry == "S_enantiomer"]
        if len(amides) >= 2:
            worst, best = amides[0], amides[-1]
            print(f"agent choice (long-acting amides): {worst.drug} (margin {worst.cns_to_cvs_margin}, "
                  f"{worst.cns_to_cvs_fold:.1f}x) vs {best.drug} (margin {best.cns_to_cvs_margin}, "
                  f"{best.cns_to_cvs_fold:.1f}x) — at SIMILAR local-anesthetic potency the cardiovascular "
                  "margin differs by stereochemistry, which is why the S-enantiomers displaced racemic "
                  "bupivacaine where cardiotoxicity risk dominates.")
        print("RESEARCH/EDUCATION — comparative only; no dose is ranked or recommended (v0.6 §7).")
        return 0
    if not args.drug:
        print("provide --drug <name> (or --agents for the cross-agent cardiotoxicity comparison)",
              file=sys.stderr)
        return 2
    mid = args.drug if "." in args.drug else f"local_anesthetics.{args.drug.replace(' ', '_')}.systemic"
    try:
        rows = site_comparison(ds, mid, dose_mg=args.dose, t_min=args.tmax)
    except (KeyError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 2
    if not rows:
        print("no simulable sites (ranks curated but ka magnitudes absent)", file=sys.stderr)
        return 2
    model = ds[mid]
    has_thr = model.has_toxicity_thresholds
    tail = "thresholds shown as RANGES via the double-uncertainty view (v0.6 LA1)" if has_thr else \
        "NO toxicity threshold is drawn (v0.6 LA0)"
    print(f"{rows[0].drug} — systemic plasma concentration after {args.dose:g} mg, by injection site")
    print(f"  tier {rows[0].tier}  (systemic concentration only — NOT block efficacy; {tail})")
    print(f"{'site':18s} {'rank':>4s} {'ka(1/min)':>9s} {'Cmax(ug/mL)':>12s} {'Tmax(min)':>10s}")
    for r in rows:
        print(f"{r.site:18s} {r.rank:>4d} {r.ka:>9.2f} {r.cmax:>12.3f} {r.tmax_min:>10.1f}")
    hi, lo = rows[0], rows[-1]
    print(f"site dominance: the SAME {args.dose:g} mg gives Cmax {hi.cmax:.2f} at {hi.site} vs "
          f"{lo.cmax:.2f} ug/mL at {lo.site} ({hi.cmax / lo.cmax:.1f}x), peaking {lo.tmax_min - hi.tmax_min:.0f} "
          "min later — so the mg/kg ceiling is the wrong mental model (absorption is site-driven).")

    if not has_thr:
        return 0

    # v0.6 LA1 — the double-uncertainty view. Always a RANGE, never a line.
    print("\ntoxicity thresholds (RANGES, never lines — v0.6 LA1):")
    for th in sorted(model.toxicity_thresholds, key=lambda t: t.low):
        print(f"  {th.endpoint:18s} {th.basis:12s} [{th.low:g}, {th.high:g}] {th.units}"
              f"  (tier {th.tier}, individual variability {th.individual_variability or '?'})")
    site = args.site
    if site is None:
        print("\npass --site <site> for the double-uncertainty view (predicted concentration vs "
              "the threshold ranges, free trace, dominant uncertainty). RESEARCH/EDUCATION — NOT a dosing tool.")
        return 0
    try:
        du = double_uncertainty(ds, mid, site=site, dose_mg=args.dose, t_min=args.tmax)
    except (KeyError, ValueError) as e:
        print(str(e), file=sys.stderr)
        return 2
    print(f"\ndouble-uncertainty view — {du.drug} {du.dose_mg:g} mg at {du.site} (tier {du.tier}):")
    if du.cardiotoxicity:
        cc = du.cardiotoxicity
        print(f"  cardiotoxicity: {cc.get('rank')} ({cc.get('stereochemistry')}) · "
              f"CNS-to-CVS margin {cc.get('cns_to_cvs_margin')} (v0.6 LA2; `hypnos la --agents` to compare)")
    fp = f"{du.peak_free:.4f}" if du.peak_free is not None else "n/a"
    print(f"  predicted peak: total {du.peak_total:.3f} ug/mL  ·  free {fp} ug/mL  ·  Tmax {du.tmax_min:.0f} min")
    print(f"  {'endpoint':18s} {'basis':12s} {'threshold range':>18s} {'predicted peak':>15s}  position")
    for e in du.endpoints:
        rng = f"[{e.low:g}, {e.high:g}] {e.units}"
        pk = f"{e.predicted_peak:.3f}" if e.predicted_peak is not None else "n/a"
        print(f"  {e.endpoint:18s} {e.basis:12s} {rng:>18s} {pk:>15s}  {e.position}")
    print(f"  >> {du.dominant_uncertainty}")
    for w in du.warnings:
        print("  - " + w)
    return 0


def cmd_washin(args) -> int:
    ds = load()
    vent = dict(t_min=args.t, alveolar_ventilation=args.valv, frc=args.frc, cardiac_output=args.co)
    if args.agent:
        agent = args.agent if "." in args.agent else f"volatiles.{args.agent.replace(' ', '_')}.mac"
        try:
            r = washin(ds, agent, **vent)
        except (KeyError, ValueError) as e:
            print(str(e), file=sys.stderr)
            return 2
        print(f"agent: {r.agent_id}   tier {r.tier}   blood:gas lambda = {r.blood_gas:g}")
        print(f"single-compartment alveolar wash-in (V_A {r.alveolar_ventilation:g}, "
              f"FRC {r.frc:g}, Q {r.cardiac_output:g} L/min):")
        print(f"  early FA/FI plateau (the wash-in knee): {r.plateau:.3f}")
        print(f"  time constant tau: {r.tau_min:.2f} min")
        print(f"  FA/FI at {r.t_min:g} min: {r.fa_fi:.3f}")
        print("comparative/education-grade; lower lambda -> faster wash-in. NOT a per-patient predictor.")
        return 0
    rows = washin_comparison(ds, **vent)
    print("Inhalational wash-in (FA/FI) -- single-compartment alveolar model, sorted fastest-first.")
    print(f"{'agent':14s} {'blood:gas':>10s} {'plateau':>8s} {'tau(min)':>9s} {f'FA/FI@{args.t:g}min':>12s}")
    for r in rows:
        name = r.agent_id.split(".")[1]
        print(f"{name:14s} {r.blood_gas:10g} {r.plateau:8.3f} {r.tau_min:9.2f} {r.fa_fi:12.3f}")
    print("Lower blood:gas solubility -> higher plateau -> faster wash-in (des/N2O fast, iso slow).")
    return 0


def cmd_washout(args) -> int:
    ds = load()
    vent = dict(t_min=args.t, alveolar_ventilation=args.valv, frc=args.frc, cardiac_output=args.co)
    if args.agent:
        agent = args.agent if "." in args.agent else f"volatiles.{args.agent.replace(' ', '_')}.mac"
        try:
            r = washout(ds, agent, **vent)
        except (KeyError, ValueError) as e:
            print(str(e), file=sys.stderr)
            return 2
        print(f"agent: {r.agent_id}   tier {r.tier}   blood:gas lambda = {r.blood_gas:g}")
        print(f"single-compartment alveolar wash-out (V_A {r.alveolar_ventilation:g}, "
              f"FRC {r.frc:g}, Q {r.cardiac_output:g} L/min):")
        print(f"  early elimination floor (= 1 - wash-in plateau): {r.floor:.3f}")
        print(f"  time constant tau: {r.tau_min:.2f} min")
        print(f"  FA/FA0 at {r.t_min:g} min: {r.fa_fa0:.3f}")
        print("comparative/education-grade; lower lambda -> lower floor -> faster emergence. NOT a per-patient predictor.")
        return 0
    rows = washout_comparison(ds, **vent)
    print("Inhalational wash-out (FA/FA0) -- single-compartment alveolar model, sorted fastest-first.")
    print(f"{'agent':14s} {'blood:gas':>10s} {'floor':>8s} {'tau(min)':>9s} {f'FA/FA0@{args.t:g}min':>13s}")
    for r in rows:
        name = r.agent_id.split(".")[1]
        print(f"{name:14s} {r.blood_gas:10g} {r.floor:8.3f} {r.tau_min:9.2f} {r.fa_fa0:13.3f}")
    print("Lower blood:gas solubility -> lower floor -> faster wash-out (des/N2O fast, iso slow).")
    return 0


def cmd_export(args) -> int:
    ds = load()
    if args.format not in FORMATS:
        print(f"unknown format {args.format!r}; choose from {FORMATS}", file=sys.stderr)
        return 2
    out = Path(args.output)
    models = [ds[args.model]] if args.model else list(ds)
    pk_models = [m for m in models if m.purpose == "pk"]

    # COMBINE .omex archive (binary): single archive if --output ends in .omex, else per-model.
    if args.format == "omex":
        if str(out).endswith(".omex"):
            out.parent.mkdir(parents=True, exist_ok=True)
            sel = [m for m in pk_models if m.kernel_implemented] if not args.model else pk_models
            out.write_bytes(combine.build_dataset_archive(ds, sel))
            print(f"wrote dataset archive {out} ({len(sel)} model(s))")
            return 0
        out.mkdir(parents=True, exist_ok=True)
        written = 0
        for m in pk_models:
            (out / combine.model_filename(m)).write_bytes(combine.build_model_archive(m, ds))
            written += 1
        print(f"wrote {written} omex archive(s) to {out}/")
        return 0

    # Dataset-level text exports (whole-dataset by default, per-model with --model):
    # BibTeX citation library and the flat-parameter CSV.
    if args.format in ("bibtex", "csv") and not args.model:
        if args.format == "bibtex":
            default_name, text = "citations.bib", bibtex.build(ds, models)
            target = out if str(out).endswith(".bib") else out / default_name
        else:
            default_name, text = "parameters.csv", csv_flat.build(ds, models)
            target = out if str(out).endswith(".csv") else out / default_name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        print(f"wrote {target}")
        return 0

    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for m in pk_models:
        fname, text = export_model(args.format, m, ds)
        (out / fname).write_text(text, encoding="utf-8")
        written += 1
    print(f"wrote {written} {args.format} file(s) to {out}/")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hypnos", description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("version").set_defaults(func=cmd_version)
    sub.add_parser("validate").set_defaults(func=cmd_validate)
    sub.add_parser("info").set_defaults(func=cmd_info)
    sub.add_parser("status", help="verification-coverage report + what to verify next").set_defaults(func=cmd_status)

    mlp = sub.add_parser("models", help="list models (optionally filtered by drug)")
    mlp.add_argument("--drug", default=None)
    mlp.set_defaults(func=cmd_models)

    pp = sub.add_parser("performance", help="published predictive-performance metrics (MDPE/MDAPE/…)")
    pp.add_argument("--drug", default=None)
    pp.set_defaults(func=cmd_performance)

    vp = sub.add_parser("verify", help="field-by-field verification checklist for a model")
    vp.add_argument("model")
    vp.add_argument("--markdown", action="store_true", help="emit a copy-pasteable markdown checklist")
    vp.set_defaults(func=cmd_verify)

    sp = sub.add_parser("simulate", help="forward-simulate one model")
    sp.add_argument("model")
    sp.add_argument("--pd", default=None, help="optional PD model id (e.g. pd_effect.propofol.bis_sigmoid)")
    sp.add_argument("--bands", action="store_true",
                    help="draw a seeded prediction band (Cp/Ce, and the effect band when --pd is set)")
    sp.add_argument("--percentile", default="5,95", help="band percentile pair, e.g. '5,95'")
    sp.add_argument("--samples", type=int, default=2000, help="Monte-Carlo draws")
    sp.add_argument("--seed", type=int, default=7, help="RNG seed (bands are byte-reproducible)")
    _add_patient_args(sp)
    sp.set_defaults(func=cmd_simulate)

    cp = sub.add_parser("compare", help="model-divergence comparison")
    cp.add_argument("--drug", required=True)
    cp.add_argument("--bands", action="store_true",
                    help="draw seeded prediction bands + uncertainty-aware divergence (v0.2)")
    cp.add_argument("--percentile", default="5,95", help="band percentile pair, e.g. '5,95'")
    cp.add_argument("--samples", type=int, default=2000, help="Monte-Carlo draws per model")
    cp.add_argument("--seed", type=int, default=7, help="RNG seed (bands are byte-reproducible)")
    _add_patient_args(cp)
    cp.set_defaults(func=cmd_compare)

    ip = sub.add_parser("interact", help="propofol-remifentanil response-surface synergy")
    ip.add_argument("--propofol", default="hypnotics_iv.propofol.schnider_1998")
    ip.add_argument("--remifentanil", default="opioids.remifentanil.minto_1997")
    ip.add_argument("--surface", default="interactions.propofol_remifentanil.greco_bis")
    _add_patient_args(ip)
    ip.set_defaults(func=cmd_interact)

    tp = sub.add_parser("tpeak", help="time to peak effect (onset) for a PK model with a ke0 link")
    tp.add_argument("model")
    _add_patient_args(tp)
    tp.set_defaults(func=cmd_tpeak)

    dp = sub.add_parser("decrement", help="plasma decrement time after a constant-rate infusion (offset)")
    dp.add_argument("model")
    dp.add_argument("--infusion", default=None, help="infusion rate, e.g. '6 mg/kg/h' (default: drug preset)")
    dp.add_argument("--duration", type=float, default=60.0, help="infusion duration (min)")
    dp.add_argument("--fraction", type=float, default=0.5, help="fractional decrement (0-1)")
    _add_patient_args(dp)
    dp.set_defaults(func=cmd_decrement)

    mp = sub.add_parser("mac", help="age-corrected MAC for a volatile agent")
    mp.add_argument("--agent", required=True, help="agent name (sevoflurane) or full model id")
    mp.add_argument("--age", type=float, required=True)
    mp.add_argument("--end-tidal", type=float, default=None, dest="end_tidal",
                    help="end-tidal concentration (vol%%) -> MAC fraction")
    mp.add_argument("--n2o", type=float, default=None, help="nitrous oxide end-tidal (vol%%), additive")
    mp.set_defaults(func=cmd_mac)

    wp = sub.add_parser("washin", help="inhalational wash-in (FA/FI) — solubility-driven uptake")
    wp.add_argument("--agent", default=None, help="agent name or model id (omit for a comparison table)")
    wp.add_argument("--t", type=float, default=3.0, help="time point for FA/FI (min)")
    wp.add_argument("--valv", type=float, default=4.0, help="alveolar ventilation (L/min)")
    wp.add_argument("--frc", type=float, default=2.5, help="functional residual capacity (L)")
    wp.add_argument("--co", type=float, default=5.0, help="cardiac output (L/min)")
    wp.set_defaults(func=cmd_washin)

    lap = sub.add_parser("la", help="local-anesthetic systemic concentration by injection site + "
                                    "the double-uncertainty threshold view + the cardiotoxicity "
                                    "agent comparison (v0.6 LA0/LA1/LA2)")
    lap.add_argument("--drug", default=None, help="lidocaine | bupivacaine | levobupivacaine | "
                     "ropivacaine (or full model id)")
    lap.add_argument("--dose", type=float, default=150.0, help="dose in mg")
    lap.add_argument("--tmax", type=float, default=90.0, help="horizon (min)")
    lap.add_argument("--site", default=None, help="injection site for the double-uncertainty view "
                     "(e.g. lumbar_epidural); omit to list the site table + threshold ranges")
    lap.add_argument("--agents", action="store_true", help="cross-agent cardiotoxicity comparison "
                     "(stereochemistry / CNS-to-CVS margin); v0.6 LA2")
    lap.set_defaults(func=cmd_la)

    op = sub.add_parser("washout", help="inhalational wash-out (FA/FA0) — solubility-driven emergence")
    op.add_argument("--agent", default=None, help="agent name or model id (omit for a comparison table)")
    op.add_argument("--t", type=float, default=3.0, help="time point for FA/FA0 (min)")
    op.add_argument("--valv", type=float, default=4.0, help="alveolar ventilation (L/min)")
    op.add_argument("--frc", type=float, default=2.5, help="functional residual capacity (L)")
    op.add_argument("--co", type=float, default=5.0, help="cardiac output (L/min)")
    op.set_defaults(func=cmd_washout)

    ep = sub.add_parser("export", help="export models to a pharmacometric format")
    ep.add_argument("--format", required=True, choices=FORMATS)
    ep.add_argument("--output", required=True)
    ep.add_argument("--model", default=None, help="single model id (default: all PK models)")
    ep.set_defaults(func=cmd_export)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
