#!/usr/bin/env python3
"""Generate notebooks/03_local_anesthetics.ipynb — the v0.6 reference notebook.

Hand-writing ipynb JSON is error-prone, so the notebook is *built* from this
script (run once; the committed artifact is the .ipynb). Every code cell calls
only the shipped, tested ``hypnos.la`` API, so the notebook is a deterministic
projection of the dataset/kernels — exactly like the exports and figures — and
is executed in CI via ``nbmake`` so it cannot rot. Plotting is optional and
guarded (skipped if matplotlib is absent), matching notebooks 01/02.
"""
from __future__ import annotations

import json
from pathlib import Path

CELLS = []


def md(text: str) -> None:
    CELLS.append({"cell_type": "markdown", "id": f"c{len(CELLS)}",
                  "metadata": {}, "source": _lines(text)})


def code(text: str) -> None:
    CELLS.append({"cell_type": "code", "id": f"c{len(CELLS)}", "metadata": {},
                  "execution_count": None, "outputs": [], "source": _lines(text)})


def _lines(text: str):
    lines = text.strip("\n").split("\n")
    return [ln + "\n" for ln in lines[:-1]] + [lines[-1]]


md(r"""
# Hypnos reference notebook: the local-anesthetic subsystem (v0.6)

**NOT FOR CLINICAL USE.** Research / education / simulation only. This view computes **no**
dose, ceiling, margin, or "is this safe?" answer — it is a forward, comparative research
instrument (v0.6 §7).

This notebook is executed in CI (via `nbmake`) so it cannot rot. It walks the four phases of
the local-anesthetic subsystem, each from the live, tested `hypnos.la` kernels:

- **LA0 — site-of-injection dominance.** Systemic absorption is *site-driven, not
  milligram-driven*: the same dose gives wildly different peaks by injection site, so the
  mg/kg ceiling is the wrong mental model.
- **LA1 — the double-uncertainty view.** Predicted plasma concentration against the
  toxicity-threshold *ranges* (never lines), whose honest punchline is that the threshold
  uncertainty *dwarfs* the PK spread.
- **LA2 — stereochemistry & cardiotoxicity.** Why agent choice changes the cardiovascular
  margin at a similar CNS threshold (racemic bupivacaine vs the S-enantiomers).
- **LA3 — binding saturation.** Total concentration under-predicts free-drug (toxic)
  concentration exactly when risk is highest, because plasma-protein binding saturates.

The design thesis: *the honest science and the safety posture are the same thing.*
""")

code(r"""
import numpy as np
import hypnos
from hypnos.la import (
    site_comparison, double_uncertainty, cardiotoxicity_comparison,
    free_concentration, saturable_free_concentration,
)

ds = hypnos.load()
assert hypnos.validate_dataset(ds) == [], "dataset must be valid"

la_models = sorted(m.id for m in ds if m.subsystem == "local_anesthetics")
print("local-anesthetic models:")
for mid in la_models:
    print(" ", mid, "(tier", ds[mid].tier + ")")
""")

md(r"""
## LA0 — site-of-injection dominance

The headline safety message *is* the science: the same dose of the same drug produces
materially different peak plasma concentrations depending on the vascularity of the
injection site (the robust rank intercostal > caudal/epidural > brachial plexus >
subcutaneous). `site_comparison` runs the forward Bateman kernel at every curated site.
""")

code(r"""
BUPI = "local_anesthetics.bupivacaine.systemic"
rows = site_comparison(ds, BUPI, dose_mg=100.0)
print(f"bupivacaine 100 mg by site (tier {rows[0].tier}):")
print(f"{'site':18s} {'rank':>4s} {'Cmax(ug/mL)':>12s} {'Tmax(min)':>10s}")
for r in rows:
    print(f"{r.site:18s} {r.rank:>4d} {r.cmax:>12.3f} {r.tmax_min:>10.1f}")
hi, lo = rows[0], rows[-1]
print(f"\nSame 100 mg -> Cmax {hi.cmax:.2f} at {hi.site} vs {lo.cmax:.2f} at {lo.site} "
      f"({hi.cmax/lo.cmax:.1f}x). The mg/kg ceiling is the wrong model.")
""")

code(r"""
# Optional plot (skipped if matplotlib is absent; the curated figure lives in docs/images/).
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4))
    for r in rows:
        ax.plot(r.t, r.cp, lw=2, label=f"{r.site} (Cmax {r.cmax:.2f})")
    ax.set_xlabel("time (min)"); ax.set_ylabel("total plasma conc. (ug/mL)")
    ax.set_title("Bupivacaine 100 mg — site-of-injection dominance"); ax.legend(fontsize=8)
    plt.show()
except ImportError:
    print("matplotlib not installed; skipping plot")
""")

md(r"""
## LA1 — the double-uncertainty view

The obvious LA plot — concentration vs a toxicity line — is a dosing calculator. Done in
Hypnos's idiom it is the opposite: a **double-uncertainty** instrument. The toxicity
thresholds are curated as **ranges** (the schema forbids a single-value threshold), and the
view names which uncertainty dominates. For LA the threshold band almost always dwarfs the
PK spread — so *no single safe-concentration line is defensible*, and that conclusion is the
safety message.
""")

code(r"""
du = double_uncertainty(ds, BUPI, site="lumbar_epidural", dose_mg=100.0)
print(f"{du.drug} {du.dose_mg:g} mg at {du.site} (tier {du.tier})")
print(f"predicted peak: total {du.peak_total:.3f} ug/mL", end="")
if du.peak_free is not None:
    print(f"  free {du.peak_free:.4f} ug/mL", end="")
print(f"  Tmax {du.tmax_min:.0f} min\n")
print(f"{'endpoint':18s} {'basis':12s} {'range(ug/mL)':>14s} {'peak':>8s}  position")
for e in du.endpoints:
    pk = f"{e.predicted_peak:.3f}" if e.predicted_peak is not None else "n/a"
    print(f"{e.endpoint:18s} {e.basis:12s} {f'[{e.low:g},{e.high:g}]':>14s} {pk:>8s}  {e.position}")
print("\n>>", du.dominant_uncertainty)
""")

md(r"""
The fold-range comparison is on a single multiplicative scale — the widest threshold band's
`high/low` vs the across-site `Cmax` spread — so the readout can never contradict itself.
""")

md(r"""
## LA2 — stereochemistry & the agent-choice cardiotoxicity margin

Racemic bupivacaine is the most cardiotoxic common amide ("fast-in/slow-out" avid cardiac
Na-channel binding); its S-enantiomers (levobupivacaine, ropivacaine) were developed to
widen the cardiovascular margin at similar potency. `cardiotoxicity_comparison` ranks the
agents with a **numeric** CNS-to-CVS fold-margin that is monotone with the qualitative class.
""")

code(r"""
print(f"{'drug':16s} {'class':>13s} {'stereochem':>13s} {'margin':>9s} {'fold':>6s}")
for a in cardiotoxicity_comparison(ds):
    fold = f"{a.cns_to_cvs_fold:.1f}x" if a.cns_to_cvs_fold is not None else "n/a"
    print(f"{a.drug:16s} {a.rank:>13s} {a.stereochemistry:>13s} {a.cns_to_cvs_margin:>9s} {fold:>6s}")
print("\nMost cardiotoxic first; the fold-margin widens with the (safer) qualitative class.")
""")

md(r"""
## LA3 — binding saturation: total under-predicts free when risk is highest

Amide LAs are highly, *saturably* bound (chiefly to a1-acid glycoprotein), and only the
**free** fraction is toxic. As total concentration rises the binding saturates, so the free
fraction climbs non-linearly. `saturable_free_concentration` is an exact capacity-limited
(Langmuir 1:1) kernel whose affinity is **pinned by the curated low-concentration
`fraction_bound`** — the only new curated quantity is the binding capacity (Tier-D
illustrative). The qualitative rise is the documented fact; the magnitude is not a fitted
prediction.
""")

code(r"""
pb = ds.drug("bupivacaine")["protein_binding"]
fb = pb["fraction_bound"]
cap = pb["free_fraction_model"]["binding_capacity_ug_ml"]
print(f"bupivacaine: fraction_bound {fb}, capacity {cap} ug/mL (Tier-{pb['free_fraction_model']['tier']})\n")
print(f"{'total':>7s} {'free_nl':>8s} {'ff%':>6s} {'free_lin':>9s} {'gap':>6s}")
for ct in [0.5, 2.0, 5.0, 10.0, 20.0]:
    fnl = float(saturable_free_concentration(np.array([ct]), fb, cap)[0])
    lin = ct * (1 - fb)
    print(f"{ct:>7.1f} {fnl:>8.3f} {100*fnl/ct:>5.1f}% {lin:>9.3f} {fnl/lin:>5.1f}x")
print("\nThe linear assumption (flat free fraction) under-predicts free-drug risk at high total.")
""")

code(r"""
# the curated model is picked up automatically by free_concentration; a non-saturable drug
# (lidocaine) correctly carries no model and stays linear.
fc_bupi = free_concentration(np.array([2.0, 20.0]), ds.drug("bupivacaine")["protein_binding"])
fc_lido = free_concentration(np.array([2.0, 20.0]), ds.drug("lidocaine")["protein_binding"])
print("bupivacaine free model:", fc_bupi.model, "(saturable)")
print("lidocaine   free model:", fc_lido.model, "(non-saturable -> linear, no fabricated curve)")
""")

md(r"""
## The thesis

Every LA view here is forward and comparative. The subsystem is designed so that the obvious
forbidden question — *"what is the maximum safe dose?"* — has no answer to read off:
thresholds are ranges, the dominant uncertainty is named, the free-fraction gap is shown, and
agent choice is compared, never ranked into a dose. **The honest science and the safety
posture are the same thing.**
""")

nb = {
    "cells": CELLS,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = Path(__file__).resolve().parent.parent / "notebooks" / "03_local_anesthetics.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
print(f"wrote {out} ({len(CELLS)} cells)")
