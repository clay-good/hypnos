#!/usr/bin/env python3
"""Deterministically regenerate all generated artifacts: exports + figures.

Exports are pure projections of the dataset (no plotting deps). Figures require
matplotlib (an optional dev dependency) and are skipped with a notice if it is
not installed, so this script always runs the export half in CI.

    python scripts/regenerate.py            # exports + figures (if matplotlib)
    python scripts/regenerate.py --exports  # exports only

The dataset is the single source of truth; everything this script writes is a
deterministic function of it (and of the dataset version). Re-running on the
same dataset produces byte-identical exports.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import hypnos
from hypnos.export import FORMATS, bibtex, combine, csv_flat, export_model

ROOT = Path(__file__).resolve().parent.parent
EXPORTS = ROOT / "exports"
# Figures regenerate into exports/ (gitignored) so a run never clobbers the
# curated, committed figures in docs/images/. Exports themselves are the
# byte-deterministic reproducibility guarantee.
IMAGES = EXPORTS / "figures"


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


def regenerate_figures(ds) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("figures: matplotlib not installed; skipped (pip install matplotlib)")
        return False
    IMAGES.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 60, 361)
    sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]

    # 1) divergence (elderly vs obese)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    colors = {"marsh_1991": "#1f77b4", "schnider_1998": "#d62728"}
    for ax, patient, title in [
        (axes[0], dict(age=72, weight=60, height=162, sex="F"), "Elderly (72 y, 60 kg, F)"),
        (axes[1], dict(age=40, weight=140, height=172, sex="M"), "Obese (40 y, 140 kg, M)"),
    ]:
        cmp = hypnos.compare(ds, drug="propofol", patient=patient, schedule=sched, t=t)
        for r in cmp.included:
            ax.plot(t, r.ce, color=colors.get(r.model_id.split(".")[-1], "#333"), lw=2.2,
                    label=f"{r.model_id.split('.')[-1]} (tier {r.tier})")
        for e in cmp.excluded:
            ax.plot(t, e["result"].ce, color="#bbb", lw=2, ls="--",
                    label=f"{e['model_id'].split('.')[-1]} (tier {e['tier']}, greyed)")
        d = cmp.divergence.get("ce", {})
        sub = f"peak div {d.get('max_abs',0):.1f} ug/mL" if d else "single model"
        ax.set_title(f"{title}\n{sub}", fontsize=11)
        ax.set_xlabel("time (min)"); ax.legend(fontsize=8); ax.grid(alpha=0.25)
    axes[0].set_ylabel("propofol effect-site (ug/mL)")
    fig.suptitle("Hypnos model-divergence view  ·  NOT FOR CLINICAL USE", y=1.02)
    fig.tight_layout(); fig.savefig(IMAGES / "divergence.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    print(f"figures: regenerated divergence.png in {IMAGES}/ (curated docs/images left untouched)")
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
