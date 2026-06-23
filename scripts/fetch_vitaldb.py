#!/usr/bin/env python3
"""Fetch a VitalDB cohort and compute the Hypnos external-validation leaderboard (v0.4 VE1–VE3).

VitalDB (https://vitaldb.net) is an **open** high-fidelity intra-operative database; its
cases carry the TCI pump's propofol delivery and the *measured* bispectral index (BIS). This
harness pulls a cohort **locally**, caches the *derived* per-subject cohort (demographics +
reconstructed infusion schedule + measured BIS — never the raw waveforms), and runs the tested
Varvel engine across **every eligible propofol PK→BIS stack at once**, **envelope-stratified**,
to produce a reproducible, apples-to-apples **leaderboard** plus a committable report.

End to end (run locally — NOT in CI):

    pip install vitaldb pandas
    python scripts/fetch_vitaldb.py --n 40 --seed 7
    python scripts/fetch_vitaldb.py --use-cache --seed 7      # offline re-run, byte-identical

What is and isn't committed (spec §3):
  * the raw cohort and the derived per-subject cache land in ``data/vitaldb/`` (gitignored);
  * only the **aggregate, derived leaderboard** (``docs/validation/vitaldb_bis_leaderboard.{json,md}``)
    is committable — no patient-level data, only MDPE/MDAPE/wobble/divergence per model + a manifest.

Scientific caveats (domain-review items, surfaced in the report, NOT verified facts):
  * the independent observation is **measured BIS**, never the pump's own predicted Ce;
  * a **propofol-only** PK→BIS stack ignores remifentanil's synergistic BIS deepening, so expect a
    systematic **over-prediction of BIS** (under-prediction of depth) in balanced-anaesthesia cases —
    a real, interpretable bias, not a bug; the leaderboard *ranking* across PK models is still
    apples-to-apples because every model shares the same cohort and the same PD link;
  * the track names + the PPF20 = 20 mg/mL vial concentration + the shared PD-BIS link are
    domain-review choices recorded in the manifest.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "vitaldb"                       # gitignored — raw + derived cohort stay local
REPORT_DIR = ROOT / "docs" / "validation"              # committable — aggregate metrics only
CACHE = DATA / "cohort.json"

# VitalDB track names (DeviceName/ParameterName). VERIFY against vitaldb.get_track_names()
# at runtime — device/param strings can vary by case/era.
TRK_BIS = "BIS/BIS"
TRK_PPF_RATE = "Orchestra/PPF20_RATE"   # propofol delivery rate (mL/h); PPF20 = 20 mg/mL
CLINICAL_URL = "https://api.vitaldb.net/cases"

# The eligible propofol PK→BIS stacks: every adult propofol PK model with an effect
# compartment (ke0 > 0), each composed with one SHARED PD-BIS link so the leaderboard
# isolates the PK model (apples-to-apples). PK-only models (Kataria/Paedfusor) cannot
# predict an effect-site BIS and are excluded by construction.
PK_MODELS = [
    "hypnotics_iv.propofol.eleveld_2018",
    "hypnotics_iv.propofol.marsh_1991",
    "hypnotics_iv.propofol.schnider_1998",
]


def _read_clinical(url: str = CLINICAL_URL) -> dict:
    """Per-case demographics from VitalDB's open clinical-information table (gzip-robust)."""
    import pandas as pd
    raw = urllib.request.urlopen(url, timeout=60).read()
    if raw[:2] == b"\x1f\x8b":                          # gzip magic
        raw = gzip.decompress(raw)
    df = pd.read_csv(io.BytesIO(raw))
    return {int(r["caseid"]): r for r in df.to_dict("records")}


def _fetch_cohort(n: int, interval: float) -> list:
    """Pull `n` VitalDB cases carrying BOTH measured BIS and propofol delivery; build case dicts."""
    import numpy as np
    import vitaldb

    case_ids = vitaldb.find_cases([TRK_BIS, TRK_PPF_RATE])[:n]
    clinical = _read_clinical()
    cases = []
    for cid in case_ids:
        recs = vitaldb.load_case(int(cid), [TRK_BIS, TRK_PPF_RATE], interval)
        if recs is None or len(recs) == 0:
            continue
        t = (np.arange(len(recs)) * interval).tolist()
        info = clinical.get(int(cid), {})

        def _num(key):
            v = info.get(key)
            return float(v) if v is not None and v == v else None

        cases.append({
            "id": str(int(cid)),
            "age": _num("age"), "sex": info.get("sex"),
            "height": _num("height"), "weight": _num("weight"),
            "bis": [(tt, float(v)) for tt, v in zip(t, recs[:, 0]) if v == v],
            "infusion_ml_h": [(tt, float(r)) for tt, r in zip(t, recs[:, 1]) if r == r],
        })
    return cases


def _leaderboard_markdown(lb, manifest) -> str:
    """A human-readable, citable leaderboard report (derived metrics only)."""
    lines = [
        "# Hypnos external validation — propofol PK→BIS on VitalDB",
        "",
        "> **NOT FOR CLINICAL USE.** Reproducible, envelope-stratified external validation of the "
        "propofol depth-of-anaesthesia stack against **measured BIS** on the open VitalDB cohort "
        "(v0.4 VE2/VE3). Computed by Hypnos's own reference kernels; kept strictly separate from any "
        "publisher-reported metric.",
        "",
        f"- **cohort:** {manifest['source']} · {lb.n_subjects} subjects scored "
        f"(of {manifest.get('n_cases_requested', '?')} requested) · seed {lb.seed}",
        f"- **target:** measured BIS · **shared PD link:** `{manifest['pd_model']}`",
        f"- **tracks:** {', '.join(manifest['tracks'])} · PPF20 = {manifest['propofol_mg_per_ml']} mg/mL",
        "",
        "## Leaderboard (ranked by overall MDAPE — inaccuracy; lower is better)",
        "",
        "| Rank | Model | Tier | MDPE % (bias) | MDAPE % | wobble % | in-env MDAPE % (n) | out-env MDAPE % (n) |",
        "| ---: | --- | :---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    def cell(cv):
        if cv is None:
            return "—"
        return f"{cv.population.mdape:.1f} ({cv.n_subjects})"

    for i, e in enumerate(lb.entries):
        p = e.overall.population
        lines.append(
            f"| {i + 1} | `{e.model_id.split('.')[-1]}` | {e.declared_tier} | "
            f"{p.mdpe:+.1f} | {p.mdape:.1f} | {p.wobble:.1f} | "
            f"{cell(e.in_envelope)} | {cell(e.out_envelope)} |")
    lines += [
        "",
        "## Caveats (domain-review items — these numbers are evidence, not yet verified facts)",
        "",
    ]
    for c in manifest["caveats"]:
        lines.append(f"- {c}")
    lines += [
        "",
        "## Reproduce",
        "",
        "```bash",
        "pip install vitaldb pandas",
        f"python scripts/fetch_vitaldb.py --use-cache --seed {lb.seed}   # byte-identical from the local cache",
        "```",
        "",
        "*See the companion `vitaldb_bis_leaderboard.json` for the full machine-readable "
        "record + manifest (the pinned cohort, tracks, seed, and Hypnos version).*",
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fetch a VitalDB cohort and compute the Varvel leaderboard.")
    ap.add_argument("--n", type=int, default=40, help="number of cases to pull")
    ap.add_argument("--pd-model", dest="pd_model", default="pd_effect.propofol.bis_sigmoid",
                    help="shared PD-BIS link composed with every PK model (apples-to-apples)")
    ap.add_argument("--interval", type=float, default=10.0, help="track resample interval (s)")
    ap.add_argument("--mg-per-ml", type=float, default=20.0, help="propofol vial conc (PPF20 = 20)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--use-cache", action="store_true",
                    help="score the cached cohort offline (no network); fails if no cache exists")
    ap.add_argument("--no-report", action="store_true", help="print only; do not write the report")
    args = ap.parse_args(argv)

    import hypnos
    from hypnos import cross_model_leaderboard, subjects_from_vitaldb

    ds = hypnos.load()
    DATA.mkdir(parents=True, exist_ok=True)

    # ---- 1. cohort: from the local cache (offline, reproducible) or a fresh fetch ----
    if args.use_cache:
        if not CACHE.is_file():
            print(f"no cache at {CACHE}; run once without --use-cache to fetch.", file=sys.stderr)
            return 2
        cases = json.loads(CACHE.read_text())
        print(f"loaded {len(cases)} cached case(s) from {CACHE}")
    else:
        try:
            import vitaldb  # noqa: F401
        except ImportError:
            print("vitaldb not installed — run `pip install vitaldb pandas` (local only).", file=sys.stderr)
            return 2
        print(f"fetching {args.n} VitalDB case(s) carrying {TRK_BIS} + {TRK_PPF_RATE} …")
        cases = _fetch_cohort(args.n, args.interval)
        CACHE.write_text(json.dumps(cases))
        print(f"cached {len(cases)} derived case(s) -> {CACHE} (gitignored; raw waveforms not kept)")

    subjects = subjects_from_vitaldb(cases, propofol_mg_per_ml=args.mg_per_ml)
    if not subjects:
        print("no usable cases (missing BIS or infusion after filtering).", file=sys.stderr)
        return 1

    # ---- 2. the cross-model, envelope-stratified leaderboard (VE2/VE3) ----
    candidates = [(pk, args.pd_model) for pk in PK_MODELS]
    manifest = {
        "source": "vitaldb", "n_cases_requested": args.n, "tracks": [TRK_BIS, TRK_PPF_RATE],
        "interval_s": args.interval, "propofol_mg_per_ml": args.mg_per_ml,
        "pd_model": args.pd_model, "seed": args.seed, "hypnos_version": hypnos.__version__,
        "caveats": [
            "Observation = MEASURED BIS (independent); the pump's predicted Ce is never used.",
            "Propofol-only PK->BIS stack: remifentanil synergy is NOT modelled, so BIS is "
            "systematically OVER-predicted (depth under-predicted) in balanced anaesthesia — a real, "
            "interpretable bias. The cross-model RANKING is still apples-to-apples (shared cohort + PD link).",
            "All PK models share one PD-BIS link, so differences reflect the PK model + that shared link.",
            "Track names + PPF20 vial conc + the shared PD link are domain-review items; verify per cohort.",
            "Raw VitalDB records are NOT committed (spec §3); only this aggregate leaderboard + manifest are.",
        ],
    }
    lb = cross_model_leaderboard(ds, subjects, candidates, target="bis", stratify_by_envelope=True,
                                 dataset="vitaldb", seed=args.seed, manifest=manifest)

    # ---- 3. print + write the committable report ----
    print(f"\nVitalDB external-validation leaderboard — propofol PK→BIS ({lb.n_subjects} subjects, "
          f"shared PD {args.pd_model}, seed {args.seed})")
    print(f"  {'rank model':34s} tier  MDPE%   MDAPE%   in-env  out-env")
    for i, e in enumerate(lb.entries):
        p = e.overall.population
        ie = f"{e.in_envelope.population.mdape:.1f}({e.in_envelope.n_subjects})" if e.in_envelope else "—"
        oe = f"{e.out_envelope.population.mdape:.1f}({e.out_envelope.n_subjects})" if e.out_envelope else "—"
        print(f"  {i + 1}. {e.model_id.split('.')[-1]:30s} {e.declared_tier}   "
              f"{p.mdpe:+6.1f}  {p.mdape:6.1f}   {ie:9s}  {oe}")

    if not args.no_report:
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
        (REPORT_DIR / "vitaldb_bis_leaderboard.json").write_text(json.dumps(lb.to_record(), indent=2))
        (REPORT_DIR / "vitaldb_bis_leaderboard.md").write_text(_leaderboard_markdown(lb, manifest))
        print(f"\nwrote committable report -> {REPORT_DIR}/vitaldb_bis_leaderboard.{{json,md}}")
        print("(aggregate metrics only; the per-subject cohort stayed in the gitignored data/ cache)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
