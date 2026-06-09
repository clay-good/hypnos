# Contributing to Hypnos

Hypnos is **infrastructure, not a simulator, and honest about uncertainty by
default.** Contributions are welcome, but the bar is *correctness and
provenance*, not coverage. A wrong number with a high tier is worse than no
number at all — anesthetic drugs are lethal at the wrong dose.

Please read [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) first. The
non-negotiable safety guardrails are in §10.

## The single highest-leverage contribution

**Promoting a model record from `unverified` to `verified`** by reading the
source PDF, field by field. This is worth more than adding a new model.

## Confidence tiers (see spec §5)

| Tier | Meaning |
| --- | --- |
| **A** | Externally validated across ≥2 independent populations with acceptable predictive performance; broad covariate range. |
| **B** | One robust population model with at least one external validation; performance reported and reasonable. |
| **C** | Single small/narrow derivation cohort; little/no external validation. |
| **D** | Used outside its derivation envelope, or an extrapolation, or a hypothesized link with no quantitative validation. **Not predictive.** |

**Worst input wins.** A record's tier must equal the worst tier among its
parameters. A composed simulation (PK + ke0 + PD) inherits the worst tier among
its components. `hypnos validate` enforces this.

**Envelope violations force Tier D.** If a request lands outside a model's
`applicability_envelope` or triggers a `known_failure_mode`, the result is
auto-tiered to D. You cannot get an A-looking number from extrapolation.

## PDF-verification checklist (the gate on `verified`)

A contributor opens the source PDF and confirms, field by field:

1. **Every structural parameter** (volumes, clearances, or micro-rate constants).
2. **Every covariate equation _and the exact form of any LBM/FFM term_** —
   this is where transcription errors hide (e.g. the James vs. Janmahasatian
   LBM choice, the sign of a coefficient).
3. **The derivation population and _n_.**
4. **The stated applicability range** (covariate min/max).

Only then does `extraction.review_status` move `unverified -> verified`, and
`verified_by` / `verified_date` are filled in. **LLMs may assist (drafting,
cross-checking) but never promote on their own authority.**

## Adding a model record

1. Add `dataset/models/<drug>_<author>_<year>.json` conforming to
   `dataset/schema/model.schema.json`.
2. Add its citation to `dataset/citations/` (Crossref/PubMed-verified).
3. If you implement an executable kernel, register it in
   `python/hypnos/export/registry.py` and set `kernel.implemented = true`.
   **If you cannot verify the covariate equations, leave `kernel.implemented`
   `false`** — `simulate()` will refuse the model rather than risk a
   mis-transcribed equation. This is the project's honesty stance in code.
4. Run the checks below; all must pass.

## Local checks (must be green before a PR)

```bash
pip install -e ".[dev]"
hypnos validate          # JSON Schema + integrity (citations, tiers, kernels, envelopes)
pytest                   # analytic-vs-numeric, round-trip, envelope, tier propagation
```

## What is out of scope (hard lines, not roadmap)

* Any per-patient dosing recommendation, bedside calculator, or pump-driver logic.
* Inverse control / target-achieving output (computing an infusion to *reach* a target).
* Local-anesthetic systemic-toxicity thresholds (deferred; needs its own safety pass).

If you find yourself wanting a "just enter the patient and get a dose" feature,
that is the signal the contribution has crossed the line the project is designed
to stay behind.
