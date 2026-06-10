#!/usr/bin/env python3
"""Generate notebooks/04_organ_function.ipynb — the v0.5 reference notebook.

Built deterministically from this script (run once; the committed artifact is the
.ipynb). Every code cell calls only the shipped, tested package API, so the
notebook is a deterministic projection of the dataset/kernels and is executed in
CI via ``nbmake`` so it cannot rot. Plotting is optional and guarded, matching
notebooks 01/02/03.
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
# Hypnos reference notebook: the organ-function envelope (v0.5)

**NOT FOR CLINICAL USE.** Research / education / simulation only.

This notebook is executed in CI (via `nbmake`) so it cannot rot. It reproduces the v0.5
headline: Hypnos's silence on organ-failure patients is replaced by an **explicit, cited
statement**. Almost no anesthetic PK model was fitted in hepatic / renal / cardiac failure,
yet those are the patients in whom getting the dose wrong is most dangerous — so a model with
no standing there is **greyed out and Tier-D'd with a named extrapolation**, exactly as a
too-high-BMI simulation already is. A model with *cited* mechanistic standing (remifentanil's
organ-independent esterase clearance) keeps it, with a note.

Two record kinds, never blurred (v0.5 §B2): a model *fitted in* a disease group is **evidence**;
a documented adjustment factor is an **annotation** (opt-in, Tier-D, always cited). This
notebook exercises the envelope half (S0) and the protein-binding failure mode (S1).
""")

code(r"""
import numpy as np
import hypnos
from hypnos.simulate import evaluate_safety, compare

ds = hypnos.load()
assert hypnos.validate_dataset(ds) == [], "dataset must be valid"

# A normal patient: no organ covariates declared -> the envelope is silent (unaffected).
NORMAL = dict(age=55, weight=80, height=175, sex="M")
# An organ-failure patient: end-stage liver disease.
CIRRHOTIC = dict(age=55, weight=80, height=175, sex="M", child_pugh="C")
print("normal patient :", NORMAL)
print("cirrhotic patient:", CIRRHOTIC)
""")

md(r"""
## S0 — the physiological envelope speaks

`evaluate_safety(model, patient, drug_record)` returns `(tier_floor, warnings, envelope_violated)`.
Declare `child_pugh="C"` and every model with no hepatic standing is greyed to **Tier-D** with a
named `HEPATIC EXTRAPOLATION`; the staging cut-points (Child-Pugh = chronic liver disease,
KDIGO CrCl<60, EF<40%, albumin<3.5) are *definitional*, named in code, not fitted PK.
""")

code(r"""
prop = ds["hypnotics_iv.propofol.schnider_1998"]
for label, patient in [("normal", NORMAL), ("Child-Pugh C", CIRRHOTIC)]:
    tier, warns, violated = evaluate_safety(prop, patient, ds.drug("propofol"))
    print(f"propofol Schnider | {label:13s} -> tier {tier}  (envelope_violated={violated})")
    for w in warns:
        print("    -", w[:100])
""")

md(r"""
### The remifentanil exception — cited standing, not a blanket pass

Remifentanil is cleared by non-specific plasma/tissue esterases, **organ-independently** — so
it is the one agent that keeps standing in hepatic (and renal) failure. That standing is curated
as a cited `organ_tolerance` entry (Dershwitz 1996 / Hoke 1997), so the envelope does *not*
wrongly grey it. The honesty cuts both ways: a model is penalized for silence **and** spared only
with a citation.
""")

code(r"""
remi = ds["opioids.remifentanil.minto_1997"]
tier, warns, violated = evaluate_safety(remi, CIRRHOTIC, ds.drug("remifentanil"))
print(f"remifentanil Minto | Child-Pugh C -> tier {tier}  (envelope_violated={violated})")
for w in warns:
    print("    -", w[:120])
""")

md(r"""
## The eligible-model set shrinks — the extrapolation made visible

`compare()` carries the same envelope: in an organ-failure patient the included set **shrinks**,
and every greyed model is named with its reason. The honest message is *"model X has some standing
and everything else is extrapolation"* — a research instrument for exactly the patients the
literature is quietest about.
""")

code(r"""
t = np.linspace(0, 30, 181)
sched = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
for label, patient in [("normal", NORMAL), ("Child-Pugh C", CIRRHOTIC)]:
    cmp = compare(ds, drug="propofol", patient=patient, schedule=sched, t=t)
    print(f"propofol | {label:13s}: {len(cmp.included)} in-envelope, {len(cmp.excluded)} greyed")
    for e in cmp.excluded:
        print("    grey", e["model_id"].split(".")[-1], "—", e["reasons"][0][:70])
""")

md(r"""
## S1 — the protein-binding / free-fraction failure mode

For a **binding-sensitive** drug (highly albumin-bound: propofol, fentanyl, dexmedetomidine),
hypoalbuminemia raises the *free* (active) fraction, so a **total**-concentration prediction
**under-estimates** effect. v0.5 surfaces this as a cited `BINDING-SENSITIVE` caveat — the shift
is *named, never silently modeled* (the never-invent rule, carried into the binding axis).
""")

code(r"""
HYPO = dict(age=55, weight=80, height=175, sex="M", albumin_g_dl=2.2)  # severe hypoalbuminemia
print("binding-sensitive drugs:",
      [d["name"] for d in ds.drugs.values() if (d.get("protein_binding") or {}).get("binding_sensitive")])
prop = ds["hypnotics_iv.propofol.eleveld_2018"]
_, warns, _ = evaluate_safety(prop, HYPO, ds.drug("propofol"))
print("\npropofol Eleveld | albumin 2.2 g/dL:")
for w in warns:
    if "BINDING-SENSITIVE" in w:
        print("    -", w)
""")

md(r"""
### Why local anesthetics are correctly OFF the albumin axis

LA protein binding is driven by **α1-acid glycoprotein**, not albumin — so the v0.5
hypoalbuminemia caveat must *not* fire for local anesthetics. Their free-fraction story is the
binding-*saturation* failure mode of v0.6 LA3, a different axis. Confirm the separation holds:
""")

code(r"""
bupi = ds["local_anesthetics.bupivacaine.systemic"]
_, warns, _ = evaluate_safety(bupi, dict(age=55, weight=80, albumin_g_dl=2.2), ds.drug("bupivacaine"))
fired = any("BINDING-SENSITIVE" in w for w in warns)
print("bupivacaine | albumin 2.2 g/dL -> hypoalbuminemia BINDING-SENSITIVE caveat fired:", fired)
print("(expected False — LA binding is α1-AG-driven, not albumin; binding_sensitive =",
      ds.drug("bupivacaine")["protein_binding"]["binding_sensitive"], ")")
""")

md(r"""
## The thesis

Silence reads as *"fine"*; v0.5 replaces it with an explicit, cited statement. The envelope
greys organ-failure extrapolations to Tier-D and names them; cited mechanistic standing
(remifentanil's esterase clearance) is the only thing that survives; and the binding failure
mode is surfaced with its citation, never folded invisibly into a clearance number. **No model is
silently "scaled for cirrhosis."** That is the never-invent rule applied to the most
extrapolation-prone patients in anesthesia.
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

out = Path(__file__).resolve().parent.parent / "notebooks" / "04_organ_function.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n", encoding="utf-8")
print(f"wrote {out} ({len(CELLS)} cells)")
