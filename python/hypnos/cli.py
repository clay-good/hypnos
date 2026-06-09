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
from .export import FORMATS, bibtex, combine, csv_flat, export_model
from .filter import summary
from .load import load
from .inhalational import mac as mac_eval
from .simulate import compare as compare_models
from .simulate import simulate, simulate_interaction
from .validate import validate_dataset
from .verification import checklist_markdown, model_verification, verification_summary

# Drug-appropriate default dose schedules for the CLI. A propofol regimen
# (2 mg/kg) applied to remifentanil would be a ~1000x overdose, so each drug
# gets a clinically sensible default; users can still vary the patient.
_DEFAULT_SCHEDULE = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
_DEFAULT_SCHEDULES = {
    "propofol": [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")],
    "remifentanil": [("bolus", 0.0, "1 mcg/kg"), ("infusion", 0.0, "0.25 mcg/kg/min")],
    "dexmedetomidine": [("infusion", 0.0, "6 mcg/kg/h"), ("infusion", 10.0, "0.5 mcg/kg/h")],
    "fentanyl": [("bolus", 0.0, "2 mcg/kg")],
    "rocuronium": [("bolus", 0.0, "0.6 mg/kg")],
}


def _default_schedule_for(drug: str):
    return _DEFAULT_SCHEDULES.get(drug, _DEFAULT_SCHEDULE)


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
        print(f"  - {r.model_id:42s} tier {r.tier}  {kind} peak {peak:.3f} {cu}")
    if cmp.excluded:
        print(f"excluded for envelope ({len(cmp.excluded)}):")
        for e in cmp.excluded:
            print(f"  - {e['model_id']:42s} tier {e['tier']}  ({'; '.join(e['reasons'])})")
    if cmp.unavailable:
        print(f"unavailable ({len(cmp.unavailable)}):")
        for u in cmp.unavailable:
            print(f"  - {u['model_id']:42s} ({u['reason']})")
    d = cmp.divergence.get("ce", {})
    if d:
        print(f"effect-site divergence across included models: "
              f"peak abs {d['max_abs'] * cmp.conc_factor:.3f} {cu}, peak rel {100*d['max_rel']:.0f}%")
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

    mp = sub.add_parser("mac", help="age-corrected MAC for a volatile agent")
    mp.add_argument("--agent", required=True, help="agent name (sevoflurane) or full model id")
    mp.add_argument("--age", type=float, required=True)
    mp.add_argument("--end-tidal", type=float, default=None, dest="end_tidal",
                    help="end-tidal concentration (vol%%) -> MAC fraction")
    mp.add_argument("--n2o", type=float, default=None, help="nitrous oxide end-tidal (vol%%), additive")
    mp.set_defaults(func=cmd_mac)

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
