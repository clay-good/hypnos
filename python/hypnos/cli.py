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
from .inhalational import washin, washin_comparison
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
    return {"age": args.age, "weight": args.weight, "height": args.height, "sex": args.sex}


def _add_patient_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--age", type=float, default=50.0)
    p.add_argument("--weight", type=float, default=77.0)
    p.add_argument("--height", type=float, default=177.0)
    p.add_argument("--sex", default="M")
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
    res = simulate(ds, args.model, patient=patient, schedule=schedule, t=t, pd_model=args.pd)
    return res


def cmd_simulate(args) -> int:
    res = _run_sim(args)
    print(f"model: {res.model_id}")
    print(f"tier (propagated): {res.tier}")
    u = res.concentration_unit
    print(f"Cp peak: {res.cp_peak_display:.3f} {u}   Ce peak: {res.ce_peak_display:.3f} {u}")
    if res.effect is not None:
        print(f"effect ({res.effect_label}): min {float(np.min(res.effect)):.1f}")
    if res.warnings:
        print("warnings:")
        for w in res.warnings:
            print("  - " + w)
    else:
        print("warnings: none (request within envelope)")
    return 0


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
    cmp = compare_models(ds, drug=args.drug, patient=patient,
                         schedule=_default_schedule_for(args.drug), t=t)
    cu = cmp.concentration_unit
    print(f"drug: {cmp.drug}   patient: {patient}   (concentrations in {cu})")
    print(f"included ({len(cmp.included)}):")
    for r in cmp.included:
        peak = r.ce_peak_display if r.ce_peak > 0 else r.cp_peak_display
        kind = "Ce" if r.ce_peak > 0 else "Cp"
        acc = _mdape_str(ds[r.model_id])
        print(f"  - {r.model_id:42s} tier {r.tier}  {kind} peak {peak:.3f} {cu}   {acc}")
    if cmp.excluded:
        print(f"excluded for envelope ({len(cmp.excluded)}):")
        for e in cmp.excluded:
            print(f"  - {e['model_id']:42s} tier {e['tier']}  ({'; '.join(e['reasons'])})")
    if cmp.unavailable:
        print(f"unavailable ({len(cmp.unavailable)}):")
        for u in cmp.unavailable:
            print(f"  - {u['model_id']:42s} ({u['reason']})")
    cf = cmp.conc_factor
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
    for it in mv.checklist:
        print(f"  [ ] ({it.group}) {it.label} = {it.value}")
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
    _add_patient_args(sp)
    sp.set_defaults(func=cmd_simulate)

    cp = sub.add_parser("compare", help="model-divergence comparison")
    cp.add_argument("--drug", required=True)
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
