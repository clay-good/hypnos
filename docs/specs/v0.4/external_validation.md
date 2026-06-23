# Hypnos — design spec (v0.4): the external-validation layer

**Recompute each model's predictive performance — Varvel's MDPE / MDAPE / wobble / divergence — by running Hypnos's own reference kernels against *open raw datasets* (the Open-TCI propofol concentration database; VitalDB on PhysioNet for BIS/effect), so that confidence-tier assignment becomes *numeric, reproducible, and auditable* rather than editorial; stratify every metric by whether the subject was *inside* the model's applicability envelope, turning v0.1's asserted failure modes into empirically confirmed predictions; and keep the Hypnos-computed numbers strictly separate from the publisher-reported ones, so the two can be reconciled rather than conflated.**

> **Status:** **VE0–VE3 implemented and the subsystem is feature-complete — the cross-model leaderboard runs on real VitalDB (BIS) and on Open-TCI/measured concentrations (cp); only the VPC *visual* overlay remains.** `subjects_from_opentci` adds the gold-standard plasma-concentration (`cp`) adapter (the sibling of `subjects_from_vitaldb`), and `hypnos validate` now raises the advisory **tier-falsification flag** (§5): a Hypnos-computed in-envelope MDAPE above a model's tier band flags a tier mismatch for human review (a computed metric *informs* but never *sets* a tier). The cross-model, **envelope-stratified leaderboard** (§7.1) ships as `cross_model_leaderboard` + `partition_by_envelope` (`hypnos.analysis`) and the `hypnos leaderboard` CLI; the VitalDB harness (`scripts/fetch_vitaldb.py`) is now end-to-end (fetch the open cohort → cache the *derived* per-subject cohort, gitignored → leaderboard → a committable, aggregate-only report `docs/validation/vitaldb_bis_leaderboard.{json,md}` with a pinned manifest). A first 30-case run put all three propofol PK→BIS stacks at **MDPE ≈ −33%** (measured BIS lower/deeper than a propofol-only model predicts — the remifentanil-synergy bias the adapter names), the PK-model choice spanning only ~2% MDAPE. The Varvel metric engine (§6: `performance_error`, `varvel_metrics`, the seeded population roll-up, and the adapter-agnostic `validate_against_cohort` harness), the `external_validation[]` / `validation_status` schema block (§4), and the `hypnos validate` consistency checks (§4) all ship in `hypnos.analysis` / the dataset schema, tested against hand-computed answers and end-to-end against the real kernel via a self-consistency fixture. **VE1** now adds the invocation surface the spec describes (§9): `hypnos validate-cohort` (over a generic `subjects_from_csv` cohort or a `--self-consistency` known-answer fixture) and the **VitalDB adapter** `subjects_from_vitaldb` (PD-BIS, measured BIS as the independent observation) + the local `scripts/fetch_vitaldb.py` that writes only derived metrics + a manifest (raw records never committed, §3) — all unit-tested against a synthetic VitalDB-shaped fixture. The live VitalDB fetch + committed metrics await the maintainer confirming data-use terms + domain review; the cross-model leaderboard (VE2) and envelope stratification (VE3) ship, the Open-TCI concentration adapter and tier-falsification flag ship, and only the VPC *visual* overlay remains. Additive over [v0.1](../v0.1/spec.md)–[v0.3](../v0.3/estimation_uncertainty.md): no existing record changes meaning; the computed-metrics block is a new, *generated* artifact distinct from the human-curated `predictive_performance`. Dataset/schema target `0.4.0`. (VE0 landed ahead of v0.3's E0 because it needs **no** newly-curated clinical numbers — the never-invent ethos applied to the roadmap.) This spec operationalizes a promise the project has made since v0.1: that tier assignment "can be partly **numeric** rather than purely editorial" (v0.1 §5) and that "where a paper publishes observed-vs-predicted data, the kernel output is compared and the MDPE/MDAPE recorded" (v0.1 §6, §9). Today those metrics are **curated from the source papers** — `Model.predictive_performance` / `predictive_mdape` (models.py §318/§322) read published values; `analysis.py` computes decrement times and time-to-peak but has **no** MDPE/MDAPE/wobble engine. v0.4 builds that engine and points it at open data.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default.* It adds the empirical floor under it. v0.1–v0.3 curate what models *claim* about themselves (point, envelope, BSV, estimation precision). v0.4 asks the orthogonal question only data can answer: **when run forward against real recorded cases, how well does each model actually predict — and exactly where does it break?**

---

## 1. The problem this layer solves — claims vs. evidence

Hypnos's tiers (v0.1 §5) are defined partly by quantitative thresholds: Tier A means *"acceptable prospective predictive performance (|MDPE| ≲ 10–20%, MDAPE ≲ 20–30%)."* But those numbers, today, come from whatever the source paper happened to report — measured in *that* paper's cohort, by *that* paper's methods, often only in-sample. Three problems follow:

1. **The numbers are not reproducible.** A published MDAPE is a single value with no re-runnable provenance. If two papers report different MDAPEs for the same model, Hypnos has no way to adjudicate; it can only transcribe.
2. **The numbers are not comparable across models.** Marsh's MDAPE was measured in one cohort, Eleveld's in another. A leaderboard built from heterogeneous published numbers compares methods and populations as much as models — the apples-to-oranges that the model-comparison literature spends whole papers trying to control for.
3. **The envelope claims are untested.** v0.1 asserts, e.g., that "Schnider's James-LBM term inverts above BMI ≈ 42" and auto-tiers it to D there (README, divergence view). That is a *prediction*. It has never been checked, inside Hypnos, against data showing Schnider's error actually exploding in the obese.

The field's own standard solves all three — and Hypnos already has every ingredient except the engine. **Varvel's framework** (the canonical anesthesia PK/PD validation methodology) computes, from a series of *observed* vs. *predicted* concentrations:

```
performance error:   PE_ij = 100 · (C_obs,ij − C_pred,ij) / C_pred,ij      (% , per sample j of subject i)
per-subject:
  MDPE_i      = median_j(PE_ij)                 — bias        (signed)
  MDAPE_i     = median_j(|PE_ij|)               — inaccuracy  (magnitude)
  wobble_i    = median_j(|PE_ij − MDPE_i|)      — intra-individual variability over time
  divergence_i= slope of |PE_ij| vs t_j         — drift of error with time (%/h)
population:  median (and spread) of each across subjects i
```

`C_pred` is exactly what Hypnos's reference kernels already produce (reference.py §147 `simulate`, the matrix-exponential PK solver), driven by the *recorded* dose history and the subject's covariates. The missing piece is the harness that aligns predictions to observations and computes the four metrics — and the open data to run it on.

**The deliverable** (§7): a reproducible, envelope-stratified, cross-model leaderboard plus a *published-vs-reproduced* reconciliation, computed by Hypnos's own kernels, with a pinned data manifest, so anyone can re-run it and get byte-identical numbers. That is the difference between a dataset that *asserts* its tiers and one that *earns* them in public.

---

## 2. What we validate against — open data, by validation mode

Validation needs *observations*. The two open-data modalities differ sharply in what they observe, so v0.4 defines **two validation modes**, and is explicit that a given dataset supports only one.

### 2.1 PK-concentration validation (the gold standard, scarce)

Requires **measured plasma/blood concentrations** with timestamps and the dosing history. Sources:

- **The Open-TCI / pooled propofol databases** — the assembled arterial-concentration datasets behind models like Eleveld 2018, several openly available. These are the canonical material for PK validation.
- **Other published observed-vs-predicted datasets** where the source releases the raw concentration-time data (rarer for opioids; curated per source).

Output: PE/MDPE/MDAPE/wobble/divergence on **plasma concentration** (`cp`) — the literal Varvel quantity.

### 2.2 PD-effect (BIS/TOF) validation (abundant, indirect)

Requires a **measured effect** (BIS for hypnotic depth, train-of-four for NMB) plus the drug infusion record; plasma assays are not needed. Source:

- **VitalDB** (the PhysioNet/open high-resolution perioperative database) — thousands of cases with synchronized BIS, vitals, and drug infusion-pump data.

Output: PE-style metrics on the **predicted effect** — Hypnos drives PK → effect-site (the `ke0` link) → the Hill PD model (reference.py §357 `sigmoid_emax`) from the recorded infusion, and compares the predicted BIS trajectory to the recorded BIS. This validates the *composed* PK+link+PD stack (v0.1 §5's "worst input wins" chain) end-to-end, which is closer to what users actually run.

> **Honesty about modality (inherited stance).** A BIS-validated model is not a concentration-validated model; the two are reported as distinct metrics with the modality named, never averaged into one "performance" number. VitalDB has limited plasma assays; Open-TCI has limited continuous effect — Hypnos validates each model on what is actually observable, and says which.

---

## 3. Data handling — reproducible, ethical, never redistributed

This layer touches human-subjects data, so its handling is a first-class part of the design, not an afterthought.

- **Hypnos stores derived metrics only, never raw patient records.** The repository commits the *computed* MDPE/MDAPE/wobble/divergence and a **data manifest** (dataset name, version/DOI, the subject-ID list, per-file checksums, the access date), not the concentrations or waveforms. Re-running requires the user to obtain the source data themselves under its own terms.
- **Access terms are respected.** VitalDB / PhysioNet credentialing and any data-use agreements are honored; the harness reads from a user-provided local path and refuses to bundle or upload the inputs. The manifest records *which* version was used so a reviewer with the same access reproduces the numbers exactly.
- **Determinism, same discipline as v0.2/v0.3.** Any Monte-Carlo step (e.g. propagating v0.2 bands into a VPC, §7.3) is seeded; a computed-metric record pins `(dataset_version, subject_set, seed, hypnos_version, git_commit)`. Identical inputs → byte-identical metrics. The runtime's no-wall-clock / no-unseeded-RNG rule applies here as everywhere.
- **No raw data in CI.** The metric *engine* (§6) is unit-tested in CI against tiny synthetic fixtures with known answers; the *full external validation* is a separately-invoked, locally-run command whose outputs are committed as curated artifacts (with their manifest), reviewed like any other dataset change. CI never needs credentialed data to be green.

---

## 4. Schema extension

Additive and backward-compatible (v0.2 §3 discipline). The key design decision: **computed metrics live in a new block, separate from the publisher-reported ones**, so the two are *comparable*, never *conflated*.

### 4.1 Model-level: `external_validation[]` (Hypnos-computed)

```jsonc
"external_validation": [
  {
    "dataset": "opentci_propofol",          // identifier of the open source
    "dataset_version": "2016-09",
    "dataset_doi": "10.xxxx/...",
    "mode": "pk_concentration",             // pk_concentration | pd_bis | pd_tof
    "target": "cp",                         // cp | ce | bis | tof
    "cohort": {                             // exactly which subjects entered this metric
      "n_subjects": 1031,
      "filter": "all",                      // or a covariate predicate, e.g. "bmi > 35"
      "in_envelope": true                   // §4.2 — the stratification axis
    },
    "metrics": [
      { "name": "MDPE",  "value": -3.2, "ci95": { "low": -5.1, "high": -1.4 }, "units": "%" },
      { "name": "MDAPE", "value": 18.9, "ci95": { "low": 17.0, "high": 20.8 }, "units": "%" },
      { "name": "wobble","value":  9.7, "units": "%" },
      { "name": "divergence", "value": -1.1, "units": "%/h" }
    ],
    "provenance": {                         // what makes it reproducible — NOT a human PDF verification
      "computed_by": "hypnos",
      "hypnos_version": "0.4.0",
      "git_commit": "…",
      "seed": 7,
      "manifest_checksum": "sha256:…"
    },
    "reproducible": true                    // distinguishes this from human-curated predictive_performance
  }
]
```

### 4.2 The envelope-stratification axis

Every `external_validation` entry carries `cohort.in_envelope`. The harness computes each model's metrics **twice**: over the subjects inside its `applicability_envelope` and over those outside (and at envelope-violating covariate slices, e.g. `bmi > 42` for Schnider). This makes v0.1's `known_failure_modes` *testable*: a confirmed failure mode is one where the out-of-envelope MDAPE materially exceeds the in-envelope MDAPE, by data. (This pairs with the README's existing discipline of only ever badging a model's *in-envelope* MDAPE.)

### 4.3 Schema-diff summary

| Location | New optional key | Holds |
| --- | --- | --- |
| model root | `external_validation` | array of Hypnos-computed metric sets, each with dataset, mode, cohort (+ `in_envelope`), metrics (+ CI), reproducible provenance |
| model root | `validation_status` | rollup: `"none" \| "internal_only" \| "external_pk" \| "external_pd" \| "external_both"` — UI badge + tier-falsification input (§5) |

`predictive_performance` (v0.1, human-curated, from-source) is **unchanged**. The new block sits beside it. `required` lists are unchanged. A new `$def` `validation_metric` is added; the existing `range` `$def` is reused for `ci95`.

---

## 5. Tier coupling — making tiers falsifiable (humans still gate promotion)

v0.1 §5 defines the tiers partly by MDPE/MDAPE thresholds but had no mechanical way to check a claim against them. v0.4 closes the loop *without* handing tier promotion to a machine:

- **Tier falsification (automatic, advisory).** If a model is labeled Tier A but its Hypnos-computed **in-envelope** MDAPE exceeds the Tier-A band (≳ 30%), `hypnos validate` raises a **tier-mismatch flag** — not an auto-demotion, a flag for human review. The tier becomes a *falsifiable claim* with reproducible evidence attached, exactly the scientific upgrade the project wanted.
- **Tier corroboration.** Conversely, a computed in-envelope MDAPE comfortably inside the threshold is recorded as supporting evidence in `tier_rationale`, so the editorial tier now cites a reproducible number, not just a paper.
- **Humans still promote.** Consistent with v0.1 §8/§9 ("Humans verify; LLMs do not promote"), a computed metric *informs* but never *sets* a tier. The harness produces evidence; a contributor reads it and decides. The difference from today: the evidence is now reproducible and envelope-stratified, not a single transcribed number.
- **Worst-input-wins is unchanged.** A composed simulation still inherits the worst tier among its components; external validation simply gives several of those component tiers an empirical basis.

> The distinction from [v0.3](../v0.3/estimation_uncertainty.md) matters and is preserved: estimation uncertainty is *internal* (how well the fitting data pins the model down); external validation is *external* (how well it predicts independent data). A model can be tightly estimated yet validate poorly, or loosely estimated yet generalize well. Hypnos keeps both and never lets one stand in for the other.

---

## 6. The metric engine — what gets built

A pure-NumPy addition to `analysis.py` (which today holds `time_to_peak_effect` §52 and `decrement_time` §92), plus a thin data-adapter layer. No inverse control (§10); everything runs the existing forward solver.

- **`analysis.performance_error(c_obs, c_pred) -> np.ndarray`** — the elementwise PE%, with the zero/near-zero-prediction guard the metric requires.
- **`analysis.varvel_metrics(pe, times) -> VarvelResult`** — MDPE, MDAPE, wobble, divergence for one subject; the population roll-up aggregates across subjects (median + nonparametric CI via the seeded bootstrap).
- **`analysis.validate_against_cohort(ds, model_id, adapter, *, stratify_by_envelope=True, seed) -> List[ExternalValidation]`** — for each subject: reconstruct the dosing schedule from the recorded infusion (reusing `build_dosing`, simulate.py §158), simulate forward with the subject's covariates, **interpolate the prediction to the observation timestamps**, compute PE, then the metrics; stratify by `in_envelope` via the existing `Envelope.check` (models.py §161) / `evaluate_safety` (simulate.py §211).
- **Data adapters** — one per source (`adapters/opentci.py`, `adapters/vitaldb.py`), each mapping the source's native records to a common `(covariates, dose_events, observations[(t, value, kind)])` shape. Adapters are the only source-specific code; the metric math is shared.
- **Time-alignment is explicit and tested.** Aligning a continuous prediction to discrete, irregular samples (and handling infusion-pump step records, dead time, and clock skew between pump and monitor) is where validation harnesses silently go wrong. The alignment policy is documented, unit-tested on synthetic fixtures with known answers, and recorded in the manifest.

---

## 7. The headline features — leaderboard, reconciliation, confirmed failure modes

### 7.1 The cross-model leaderboard (apples-to-apples, at last)

Run *every eligible model* for a drug against the *same* cohort of subjects, in-envelope only, and rank by computed MDAPE/MDPE. Because the cohort is identical, the comparison isolates the *model* — the controlled comparison the literature does by hand, now reproducible and open:

```text
$ hypnos validate-external --drug propofol --dataset opentci_propofol --mode pk
propofol — Open-TCI pooled (in-envelope subjects only), Hypnos 0.4.0, commit abc123, seed 7
  model              n     MDPE     MDAPE    wobble   divergence
  eleveld_2018      1031   -3.2%    18.9%     9.7%    -1.1%/h     tier A  (corroborated)
  schnider_1998      612   -7.4%    24.1%    12.0%    +0.3%/h     tier B
  marsh_1991         588  -14.8%    29.6%    13.4%    +2.1%/h     tier B
```

### 7.2 Published-vs-reproduced reconciliation (a kernel self-check)

For each model, compare Hypnos's recomputed metric to the publisher-reported `predictive_performance`. Agreement **validates the reference kernel** end-to-end on real data (a far stronger check than the synthetic round-trip). Disagreement is *itself* a finding: it flags either a transcription error in the curated parameters (the exact bug v0.1 §9 says hides in the covariate equations) or a methodology difference (sampling scheme, error definition) — both worth surfacing. The reconciliation report names which.

### 7.3 Confirmed failure modes + the VPC

- **Failure-mode confirmation.** For each `known_failure_mode`, report the in-envelope vs. out-of-envelope (or at-the-failing-slice) MDAPE side by side. Schnider's obesity claim stops being an assertion and becomes a row: *in-envelope MDAPE 24%, BMI>42 MDAPE 48% — failure mode confirmed by N subjects.* The dataset's safety claims become empirically backed.
- **Visual predictive check (VPC).** Where v0.2 prediction bands exist, overlay the curated 5–95% band on the binned observed quantiles from the cohort — the field-standard graphical validation, and what lets the v0.2 *band* tier (v0.2 §9) be partly numeric. The honesty note carries through: a VPC is only as good as the curated Ω; models with `variability_status: none` get a point-prediction goodness-of-fit instead, named as such.

---

## 8. Export & reporting

External validation is primarily a *reporting* layer, but it threads into the existing exports:

| Surface | v0.4 behavior |
| --- | --- |
| **TCI-sim JSON / dataset** | `external_validation[]` blocks pass through verbatim (lossless JSON), so a downstream consumer sees the reproducible metrics and their manifest. |
| **PharmML / SO** | The SO has a native goodness-of-fit / diagnostics section; computed MDPE/MDAPE attach there, alongside v0.3's estimation object, with provenance. |
| **NONMEM / nlmixr2 / Pumas / SBML** | The computed metrics + manifest attach as `hypnos:externalValidation` MIRIAM/RDF on the model annotation (these formats describe the model, not its post-hoc validation), with `hypnos:validationStatus`. |
| **A standalone validation report** | `hypnos validate-external --report` emits a human-readable + machine-readable report (the leaderboard, the reconciliation, the confirmed failure modes, the VPC figures) with the full manifest — a citable artifact. |

Every validation artifact inherits the universal `clinicalUse = "PROHIBITED"` annotation. Good external performance is **not** a clinical clearance and the report says so prominently (§10).

---

## 9. Validation of the validator

The metric engine is itself verified, with the project's usual layered discipline:

- **Engine unit tests (CI).** `performance_error` / `varvel_metrics` are checked against hand-computed answers on synthetic fixtures (including the textbook edge cases: a single sample, all-zero error, monotone drift for divergence), to algebraic tolerance.
- **Kernel-on-real-data check (§7.2).** The published-vs-reproduced reconciliation *is* a validation of the whole forward stack against the source's own reported metrics — the strongest available check that the curated parameters and the kernel agree with reality.
- **Adapter tests.** Each data adapter has fixtures (tiny, synthetic, committed) asserting correct dose reconstruction and time alignment, since that is where harnesses fail silently (§6).
- **Reproducibility test.** Re-running with the same manifest + seed reproduces committed metrics byte-for-byte; a drift is a CI failure (against a stored synthetic-cohort metric, not patient data).

**Human role.** A contributor confirms the *adapter mapping* against the source's data dictionary (units, what each column means, infusion-rate encoding) — the data-layer analog of the PDF parameter verification. As always, **LLMs assist but never promote**: a computed metric does not move a tier; a human reading the evidence does (§5).

---

## 10. Safety & scope guardrails (non-negotiable, inherited and extended)

External validation is forward-only and changes none of v0.1 §10 / v0.2 §10 / v0.3 §10. It adds guardrails specific to *performance numbers on real cases*, because a good MDAPE is exactly the kind of result that invites over-trust:

- **Validation uses *recorded* dose histories — it never tunes a dose.** The harness drives the forward solver from the case's actual infusion record. It does not search for the dose that would have improved the prediction, or recommend any dose. That search is inverse control, forbidden as everywhere.
- **Good external performance is not clinical validation.** A reproducible 18% MDAPE on an open cohort is a *research* statement about a *population*, not a fitness-for-use clearance for any patient or device. Every report repeats this; the `clinicalUse = "PROHIBITED"` flag is universal and prominent.
- **No patient-data redistribution; access terms honored** (§3). This is both an ethical and a legal hard line.
- **Metrics describe a cohort, not a person.** An MDAPE is a median over subjects; it says nothing about where the next individual will land — that is the v0.2 prediction band's job, and the two are never substituted for each other.

---

## 11. Worked example (illustrative — numbers pending real runs)

```text
$ hypnos validate-external --drug propofol --dataset opentci_propofol --mode pk \
                           --stratify-envelope --reconcile --seed 7
Open-TCI pooled propofol  ·  manifest sha256:9f3c…  ·  Hypnos 0.4.0  ·  commit abc123

leaderboard (in-envelope, identical cohort):
  eleveld_2018   MDPE -3.2%  MDAPE 18.9%  wobble 9.7%   divergence -1.1%/h
  schnider_1998  MDPE -7.4%  MDAPE 24.1%  wobble 12.0%  divergence +0.3%/h
  marsh_1991     MDPE -14.8% MDAPE 29.6%  wobble 13.4%  divergence +2.1%/h

reconciliation (Hypnos-computed vs publisher-reported MDAPE):
  eleveld_2018   computed 18.9%  vs published ~19–24%   AGREE (kernel corroborated)
  schnider_1998  computed 24.1%  vs published ~18%      DISAGREE — flag: check cohort/method or covariate eqs

failure-mode confirmation:
  schnider_1998  in-envelope MDAPE 24.1%  ·  BMI>42 slice MDAPE 47.8% (n=39)  -> failure mode CONFIRMED
  marsh_1991     in-envelope MDAPE 29.6%  ·  age>75 slice MDAPE 41.2% (n=58)  -> elderly degradation CONFIRMED

tier check:
  eleveld_2018  tier A  -> corroborated (in-envelope MDAPE within Tier-A band)
  no tier-mismatch flags raised
```

The reading: the leaderboard ranks the three propofol models on *one* cohort, reproducibly; the reconciliation cross-checks Hypnos's own kernel against the literature (and flags one disagreement worth a human look); and the failure-mode rows turn v0.1's *assertions* about Schnider-in-obesity and Marsh-in-the-elderly into *measurements*. Every number carries a manifest, so a reviewer re-runs and gets the same thing.

(All numbers above are illustrative placeholders. Per the project ethos, real metrics are produced only by running the engine against the actual open data under its access terms, and are committed with their manifest.)

---

## 12. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **VE0 — Metric engine + reconciliation** | ✅ engine shipped | `performance_error` / `varvel_metrics` / `pooled_performance` + the `validate_against_cohort` harness in `analysis.py`; synthetic-fixture CI tests + a real-kernel self-consistency test; schema `external_validation` / `validation_status` block + `validate` consistency checks. *Reconciliation against real published numbers needs observed data → demonstrated once an adapter (VE1/VE2) is wired.* | Hypnos can recompute Varvel metrics reproducibly and run them end-to-end against the forward solver; the schema carries the computed block, kept separate from the curated one. |
| **VE1 — PD/BIS validation (VitalDB)** | 🟡 adapter + CLI shipped; live run pending data-use terms | `hypnos validate-cohort` (the spec's "separately-invoked, locally-run command", §9) over a generic `subjects_from_csv` cohort or a `--self-consistency` known-answer fixture; the **VitalDB adapter** `subjects_from_vitaldb` (PD-BIS: measured BIS as the independent observation; `Orchestra/PPF20_RATE` reconstructed as a step-infusion schedule) + the local `scripts/fetch_vitaldb.py` writing derived metrics + a manifest to a gitignored `data/` (no raw records committed, §3); all unit-tested against a synthetic VitalDB-shaped fixture. **Pending:** the live fetch + committed metrics (maintainer confirms VitalDB data-use terms) + domain review of the track/unit choices and the propofol-only-vs-remifentanil-synergy caveat. | The depth-of-anesthesia stack is validated end-to-end on thousands of open cases, with no patient data committed. *(Engine + adapter + command shipped; the validated run is a local, terms-gated step.)* |
| **VE2 — Cross-model leaderboard (run on VitalDB)** | 🟢 shipped (BIS); Open-TCI `cp` adapter proposed | `cross_model_leaderboard` + `partition_by_envelope` + `hypnos leaderboard`: every eligible PK[→PD] stack scored on ONE cohort by the same kernels + seed (apples-to-apples), ranked by overall MDAPE; the VitalDB harness fetches → caches the derived cohort (gitignored) → writes the committable, aggregate-only `docs/validation/vitaldb_bis_leaderboard.{json,md}`. **First real run:** 30 VitalDB cases, all three propofol PK→BIS stacks at MDPE ≈ −33% (synergy-dominated), PK choice ~2% MDAPE. **Shipped:** the Open-TCI plasma-concentration (`cp`) adapter `subjects_from_opentci` for the gold-standard concentration leaderboard (`hypnos leaderboard --target cp`). | The open, reproducible propofol leaderboard exists and has run on real data. |
| **VE3 — Envelope stratification** | 🟢 shipped (stratification + report); tier-falsification flag + VPC proposed | each model re-scored on its in-/out-of-envelope sub-cohorts (`partition_by_envelope` via `Envelope.check`), carried in the leaderboard record + the `--report` markdown — turning v0.1's asserted failure modes into measured ones. **Shipped:** the `validate` tier-mismatch flag (§5, advisory). **Remaining:** VPC overlays on the v0.2 bands (a visual nicety). | The tiers are now empirically stratified and the report is a reproducible, citable artifact. |

VE0 alone is a useful, self-contained increment: even before any external dataset is wired in, recomputing each model's metrics from the curated parameters and reconciling them against the published numbers is a real, reproducible check on the dataset's own integrity.

---

## 13. Cheat sheet (target API)

```python
import hypnos
from hypnos.analysis import performance_error, varvel_metrics, validate_against_cohort
from hypnos.adapters import opentci
ds = hypnos.load()

# The metric primitives (pure NumPy; the same math the harness uses)
pe = performance_error(c_obs, c_pred)        # elementwise PE%
m  = varvel_metrics(pe, times)               # MDPE, MDAPE, wobble, divergence

# Validate a model against an open cohort (forward only; recorded doses; seeded)
cohort = opentci.load("/path/to/opentci")    # user-provided; never bundled
ev = validate_against_cohort(ds, "hypnotics_iv.propofol.eleveld_2018",
                             cohort, stratify_by_envelope=True, seed=7)
ev[0].metrics            # [{name: "MDAPE", value: 18.9, ci95: {...}}, ...]
ev[0].cohort.in_envelope # True / False stratum
ev[0].provenance         # {hypnos_version, git_commit, seed, manifest_checksum}

m = ds["hypnotics_iv.propofol.eleveld_2018"]
m.validation_status                          # "external_pk" | "external_both" | ...
m.predictive_performance                     # v0.1: publisher-reported (unchanged)
m.external_validation                        # v0.4: Hypnos-computed, reproducible (distinct)
```

```bash
hypnos validate-external --drug propofol --dataset opentci_propofol --mode pk \
                         --stratify-envelope --reconcile --seed 7 --report
hypnos validate                                  # now also raises tier-mismatch flags (advisory)
```

---

## 14. Open questions & explicit deferrals

- **Bayesian / individualized forecasting.** Maximum-a-posteriori "Bayesian TCI" individualization (updating predictions from a patient's measured samples) is a powerful validation refinement *and* edges toward per-patient adaptation — deferred, with the dosing-line reservation recorded, exactly as v0.2 §14 deferred Bayesian model averaging.
- **Effect-site validation without assays.** BIS validates the *composed* stack but cannot isolate the `ke0` link from the PD model; teasing those apart needs either concentration+effect in the same subjects or a published `t_peak` validation — tracked as a refinement, not assumed.
- **Sampling-design bias.** Sparse, opportunistic sampling (VitalDB monitor data) and rich arterial sampling (Open-TCI) weight the metrics differently; the harness records the sampling design in the manifest and does not pretend the two are interchangeable.
- **Additional drugs & datasets.** The engine is drug-agnostic; remifentanil/dexmedetomidine/volatile validation lands as suitable open datasets are identified and adapters written — coverage growth within the declared envelope, not new mechanism.
- **Relationship to the uncertainty layers.** v0.2 (BSV) and v0.3 (estimation) describe uncertainty the model *declares about itself*; v0.4 measures error *against the world*. A VPC (§7.3) is where they meet — observed scatter vs. declared band — and is the one place all three layers are checked at once. They remain conceptually distinct and separately tiered.
