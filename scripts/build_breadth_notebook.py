#!/usr/bin/env python3
"""Generate notebooks/05_iv_breadth.ipynb, the Phase C breadth reference notebook.

Built deterministically from this script (run once; the committed artifact is the .ipynb).
Every code cell calls only the shipped, tested package API, and is executed in CI via
``nbmake`` so it cannot rot. Each claim the README makes about the new drugs is asserted.
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
# Hypnos reference notebook: IV breadth (Phase C)

**NOT FOR CLINICAL USE.** Research / education / simulation only.

Executed in CI (via `nbmake`), so it cannot rot. It reproduces four claims about the Phase C
breadth: the fentanyl family, midazolam, ketamine, etomidate, and the flumazenil reversal.
Every record here is `pending_human_review`: sourced and cross-checked by an automated process,
never human-verified.

1. **The same drug in a different patient** (midazolam): healthy-volunteer vs post-CABG ICU kinetics.
2. **Pedigree is not performance** (alfentanil): a covariate population model that prospectively
   under-predicted by about half, next to an older set that was nearly unbiased.
3. **Rational opioid selection, forward**: offset vs infusion duration across the fentanyl family.
4. **Resedation, computed**: flumazenil on a midazolam infusion, from both drugs' kinetics.
""")

code(r"""
import numpy as np
import hypnos
from hypnos.presets import default_schedule_for
from hypnos.simulate import simulate_reversal

ds = hypnos.load()
assert hypnos.validate_dataset(ds) == [], "dataset must be valid"
ADULT = dict(age=40, weight=70, height=175, sex="M")
print(len(ds), "models;", sorted({m.drug_name for m in ds if m.purpose == "pk" and m.kernel_implemented}))
""")

md(r"""
## 1. Midazolam: healthy volunteers vs the ICU

Bührer 1990 (healthy men) and Zomorodi 1998 (63 post-CABG ICU patients) were printed side by
side in Zomorodi's Table 3. For a 40-year-old both are in envelope, and they disagree by more
than 100% at the widest point. For a 62-year-old ICU patient only Zomorodi is in envelope, and
the compare says model-selection risk is unquantified there, not zero.
""")

code(r"""
t = np.linspace(0, 240, 241)
sched = default_schedule_for("midazolam")
cmp = hypnos.compare(ds, drug="midazolam", patient=ADULT, schedule=sched, t=t)
for r in cmp.included:
    print(f"{r.model_id:42s} tier {r.tier}  Cp at 4 h {r.cp[-1] * 1000:6.1f} ng/mL")
rel = cmp.divergence["cp"]["max_rel"]
print(f"peak relative plasma divergence: {100 * rel:.0f}%")
assert rel > 1.0

icu = hypnos.compare(ds, drug="midazolam", patient=dict(age=62, weight=85, height=175, sex="M"),
                     schedule=sched, t=t)
print("ICU patient in envelope:", [r.model_id for r in icu.included])
assert [r.model_id for r in icu.included] == ["hypnotics_iv.midazolam.zomorodi_1998"]
""")

md(r"""
## 2. Alfentanil: pedigree is not performance

Raemer 1990 ran both sets in one computer-controlled pump. Maitre 1987 (a covariate population
model) gave MDPE +53%: measured concentrations exceeded target in 94% of samples. The older
Scott set gave +1%. Positive MDPE means the model **under-predicts**, so Maitre's curve should sit
below Scott's for the same infusion.
""")

code(r"""
sched = default_schedule_for("alfentanil")
t = np.linspace(0, 60, 121)
for mid in ("opioids.alfentanil.maitre_1987", "opioids.alfentanil.scott_1985"):
    perf = {p["metric"]: p["value"] for p in ds[mid].raw["predictive_performance"]}
    r = hypnos.simulate(ds, mid, patient=ADULT, schedule=sched, t=t)
    print(f"{mid:34s} MDPE {perf['MDPE']:+5.1f}%  MDAPE {perf['MDAPE']:4.1f}%  Cp at 60 min {r.cp[-1] * 1000:6.1f} ng/mL")
m = hypnos.simulate(ds, "opioids.alfentanil.maitre_1987", patient=ADULT, schedule=sched, t=t)
s = hypnos.simulate(ds, "opioids.alfentanil.scott_1985", patient=ADULT, schedule=sched, t=t)
assert np.all(s.cp[t >= 10] > m.cp[t >= 10]), "Maitre should predict lower (it under-predicts)"
""")

md(r"""
## 3. Rational opioid selection, forward

The 50% plasma decrement after a fixed-rate infusion (a forward quantity; classic CSHT needs
inverse control and is out of scope). Remifentanil stays flat, sufentanil stays well below
fentanyl, and fentanyl keeps climbing.
""")

code(r"""
opioids = ["opioids.remifentanil.minto_1997", "opioids.sufentanil.gepts_1995",
           "opioids.alfentanil.scott_1985", "opioids.fentanyl.shafer_1990"]
durations = [30, 120, 480]
table = {}
for mid in opioids:
    table[mid] = [hypnos.decrement_time(ds, mid, patient=ADULT, infusion="1 mcg/kg/h",
                                        duration=float(d)).decrement_min for d in durations]
    print(f"{mid:34s}", "  ".join(f"{d:>4d} min -> {x:6.1f}" for d, x in zip(durations, table[mid])))
assert table["opioids.fentanyl.shafer_1990"][-1] > 3 * table["opioids.sufentanil.gepts_1995"][-1]
assert max(table["opioids.remifentanil.minto_1997"]) < 5
""")

md(r"""
## 4. Resedation, computed from both drugs' kinetics

Fiset 1995 held midazolam constant by infusion, reversed it with flumazenil, and then let the
volunteers resedate. The competitive shift is scaled by Fiset's published reversal EC50 (Tier C,
not a binding constant). The reversed midazolam is reported as a share of the un-reversed level.
Hypnos never computes a reversal or re-dosing schedule.
""")

code(r"""
t = np.linspace(0, 300, 601)
r = simulate_reversal(ds, "reversal.flumazenil.competitive", patient=dict(age=40, weight=80, height=177, sex="M"),
                      reversal_dose="0.5 mg", agonist_dose="0.1 mg/kg/h", t=t, reversal_time_min=60)
share = r.agonist_ce_effective / np.where(r.agonist_ce > 0, r.agonist_ce, 1)
for tt in (60, 70, 120, 300):
    i = int(np.argmin(np.abs(t - tt)))
    print(f"t={tt:>3d} min  reversed midazolam = {100 * share[i]:5.1f}% of un-reversed")
assert share.min() < 0.6 and share[-1] > 0.95 and r.renarcotization
""")

md(r"""
## The thesis

Breadth here means more than more drugs. Each addition carries the evidence about when to
distrust it: a population model that under-predicted prospectively, an ICU volume nearly three
times a volunteer's, an antidote that wears off before the drug it reverses. All of it is
computed forward from cited, `pending_human_review` records, and none of it is a dose.
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

out = Path(__file__).resolve().parent.parent / "notebooks" / "05_iv_breadth.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
print(f"wrote {out} ({len(CELLS)} cells)")
