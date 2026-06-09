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
from .export import FORMATS, export_model
from .filter import summary
from .load import load
from .simulate import compare as compare_models
from .simulate import simulate
from .validate import validate_dataset

_DEFAULT_SCHEDULE = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]


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
    res = simulate(ds, args.model, patient=patient, schedule=_DEFAULT_SCHEDULE, t=t,
                   pd_model=args.pd)
    return res


def cmd_simulate(args) -> int:
    res = _run_sim(args)
    print(f"model: {res.model_id}")
    print(f"tier (propagated): {res.tier}")
    print(f"Cp peak: {res.cp_peak:.3f} ug/mL   Ce peak: {res.ce_peak:.3f} ug/mL")
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
    cmp = compare_models(ds, drug=args.drug, patient=patient, schedule=_DEFAULT_SCHEDULE, t=t)
    print(f"drug: {cmp.drug}   patient: {patient}")
    print(f"included ({len(cmp.included)}):")
    for r in cmp.included:
        print(f"  - {r.model_id:42s} tier {r.tier}  Ce peak {r.ce_peak:.3f}")
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
              f"peak abs {d['max_abs']:.3f} ug/mL, peak rel {100*d['max_rel']:.0f}%")
    return 0


def cmd_export(args) -> int:
    ds = load()
    if args.format not in FORMATS:
        print(f"unknown format {args.format!r}; choose from {FORMATS}", file=sys.stderr)
        return 2
    models = [ds[args.model]] if args.model else list(ds)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for m in models:
        if m.purpose != "pk":
            continue
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

    sp = sub.add_parser("simulate", help="forward-simulate one model")
    sp.add_argument("model")
    sp.add_argument("--pd", default=None, help="optional PD model id (e.g. pd_effect.propofol.bis_sigmoid)")
    _add_patient_args(sp)
    sp.set_defaults(func=cmd_simulate)

    cp = sub.add_parser("compare", help="model-divergence comparison")
    cp.add_argument("--drug", required=True)
    _add_patient_args(cp)
    cp.set_defaults(func=cmd_compare)

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
