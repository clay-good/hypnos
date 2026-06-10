#!/usr/bin/env python3
"""Fetch a VitalDB cohort and compute Hypnos external-validation metrics (v0.4 VE1).

VitalDB (https://vitaldb.net, open) is a high-fidelity intra-operative vital-signs
database; cases carry the TCI pump's propofol delivery and the *measured* BIS. This
script pulls a cohort **locally**, maps it through :func:`hypnos.subjects_from_vitaldb`,
and runs the tested Varvel engine (:func:`hypnos.validate_against_cohort`) to produce
the reproducible ``external_validation[]`` record + a manifest.

Run locally — NOT in CI, NOT committed with raw data:

    pip install vitaldb pandas
    python scripts/fetch_vitaldb.py --n 30 --model hypnotics_iv.propofol.eleveld_2018 \
        --pd-model pd_effect.propofol.eleveld_bis

Per spec §3 the raw records never enter the repo: outputs land in ``data/vitaldb/``
(gitignored); only the *derived metrics + manifest* are committable, and only after the
maintainer confirms VitalDB's data-use terms. The metrics validate a propofol-only
PK→BIS stack against measured BIS; remifentanil's synergistic BIS deepening is NOT
modelled here, so expect a systematic over-prediction of BIS (under-prediction of depth)
in balanced-anaesthesia cases — a real, interpretable caveat, not a bug. Treat every
number as pending domain review (the clinical track/unit choices live in the adapter).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "vitaldb"   # gitignored — raw cohort + metrics stay local

# VitalDB track names (DeviceName/ParameterName). VERIFY against vitaldb.get_track_names()
# at runtime — device/param strings can vary by case/era.
TRK_BIS = "BIS/BIS"
TRK_PPF_RATE = "Orchestra/PPF20_RATE"   # propofol delivery rate (mL/h); PPF20 = 20 mg/mL


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fetch a VitalDB cohort and compute Varvel metrics.")
    ap.add_argument("--n", type=int, default=20, help="number of cases to pull")
    ap.add_argument("--model", default="hypnotics_iv.propofol.eleveld_2018")
    ap.add_argument("--pd-model", dest="pd_model", default="pd_effect.propofol.eleveld_bis")
    ap.add_argument("--interval", type=float, default=10.0, help="track resample interval (s)")
    ap.add_argument("--mg-per-ml", type=float, default=20.0, help="propofol vial conc (PPF20 = 20)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    try:
        import vitaldb
    except ImportError:
        print("vitaldb not installed — run `pip install vitaldb pandas` (local only).", file=sys.stderr)
        return 2

    import numpy as np

    import hypnos
    from hypnos import subjects_from_vitaldb, validate_against_cohort

    ds = hypnos.load()
    OUT.mkdir(parents=True, exist_ok=True)

    # cases that carry BOTH a measured BIS and the propofol delivery track
    case_ids = vitaldb.find_cases([TRK_BIS, TRK_PPF_RATE])[: args.n]
    clinical = {int(r["caseid"]): r for r in _clinical_info(vitaldb)}
    cases = []
    for cid in case_ids:
        recs = vitaldb.load_case(cid, [TRK_BIS, TRK_PPF_RATE], args.interval)  # (n, 2) array
        if recs is None or len(recs) == 0:
            continue
        t = np.arange(len(recs)) * args.interval
        info = clinical.get(int(cid), {})
        cases.append({
            "id": str(cid),
            "age": info.get("age"), "sex": info.get("sex"),
            "height": info.get("height"), "weight": info.get("weight"),
            "bis": [(tt, v) for tt, v in zip(t, recs[:, 0]) if v == v],
            "infusion_ml_h": [(tt, v) for tt, v in zip(t, recs[:, 1]) if v == v],
        })

    subjects = subjects_from_vitaldb(cases, propofol_mg_per_ml=args.mg_per_ml)
    if not subjects:
        print("no usable cases (missing BIS or infusion after filtering).", file=sys.stderr)
        return 1
    cv = validate_against_cohort(ds, args.model, subjects, target="bis",
                                 pd_model=args.pd_model, dataset="vitaldb", seed=args.seed)
    record = cv.to_record()
    manifest = {
        "source": "vitaldb", "n_cases_requested": args.n, "n_subjects_scored": cv.n_subjects,
        "tracks": [TRK_BIS, TRK_PPF_RATE], "interval_s": args.interval,
        "propofol_mg_per_ml": args.mg_per_ml, "model": args.model, "pd_model": args.pd_model,
        "seed": args.seed,
        "caveats": [
            "Observation = measured BIS (independent); pump-predicted Ce is NOT used.",
            "Propofol-only PK->BIS stack: remifentanil synergy is NOT modelled (expect BIS over-prediction).",
            "Track names + PPF20 vial conc are domain-review items; verify per cohort.",
            "Raw VitalDB records are NOT committed (spec §3); only these derived metrics + manifest are.",
        ],
    }
    (OUT / "metrics.json").write_text(json.dumps({"manifest": manifest, "record": record}, indent=2))
    p = cv.population
    print(f"VitalDB external validation — {cv.model_id} (+ {args.pd_model}) vs measured BIS")
    print(f"  {cv.n_subjects} subject(s) · target bis · seed {cv.seed}")
    print(f"  MDPE {p.mdpe:.2f}%   MDAPE {p.mdape:.2f}%   wobble {p.wobble:.2f}%   "
          f"divergence {p.divergence:.2f}%/h")
    print(f"  wrote {OUT / 'metrics.json'} (raw data stayed local; review caveats before committing)")
    return 0


def _clinical_info(vitaldb):
    """Per-case demographics from VitalDB's open clinical-information table."""
    import pandas as pd
    return pd.read_csv("https://api.vitaldb.net/cases").to_dict("records")


if __name__ == "__main__":
    raise SystemExit(main())
