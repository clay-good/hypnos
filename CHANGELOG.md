# Changelog

All notable changes to Hypnos are documented here. The project follows the
phased roadmap in [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) §11.
Dates are ISO-8601. Versioning is [SemVer](https://semver.org/); the dataset
version is pinned in `dataset/VERSION` and stamped into every export as
`hypnos:datasetVersion`.

## [Unreleased]

### v0.2 — the population-variability layer (V0–V2; V3 partial) (2026-06-09)
Curate the **random-effects** structure of population models — between-subject
variability (η / Ω) and residual error (ε / Σ) — alongside the typical-value
parameters Hypnos already curates, propagate it as seeded Monte-Carlo prediction
bands, and upgrade the divergence view from "the point estimates differ by *X*"
to the sharper question *do the models disagree beyond their own stated
uncertainty?* Full design in [`docs/specs/v0.2/variability.md`](docs/specs/v0.2/variability.md).
Dataset/schema bump **0.1.0 → 0.2.0**; additive and backward-compatible (every
v0.1 record stays valid; the variability block is optional).

- **Schema (V0).** New optional `parameters[].variability` (per-parameter `bsv.omega2`
  — the canonical η-scale variance, the NONMEM `$OMEGA` diagonal entry — plus a
  checked-consistent `cv_percent`, optional `shrinkage_percent`/`iov`, its own tier
  and extraction), model-root `residual_error` (Σ: `proportional`/`additive`/`log`,
  scale labeled in the key name), `omega_block` (off-diagonal Ω), and a
  `variability_status` rollup (`none|partial|diagonal|full`). `required` lists
  unchanged; `additionalProperties:false` preserved.
- **Curation (V0).** **Eleveld 2018 propofol** carries the full Ω diagonal (ω² on
  V1/V2/V3/CL/Q2/Q3/ke0) + the log-additive Σ (σ=0.191), cross-checked against the
  `tci` R package (the same source Hypnos already cross-checks its kernels against),
  curated **`unverified`** with the §4-trap checklist in each `tier_rationale`.
  The variability layer carries its **own tier** (B), so the band is honestly
  labeled below the Tier-A median line it surrounds.
- **`hypnos validate` (V0).** Three new consistency checks (spec §9): `cv_percent`
  recomputes from `omega2` within tolerance (catches the variance/SD/CV% confusion,
  Trap 1); every variability/residual/omega-block citation resolves; and
  `variability_status` matches the curated contents (no `full` without an
  `omega_block`, no `none` hiding a curated ω²).
- **Bands (V1).** New `reference.sample_individual` draws a virtual individual
  (`P_i = P_typ·exp(η)`, seeded) + `residual_std`/`apply_residual` for the Σ layer.
  `simulate(..., bands=True, percentile=(5,95), samples=2000, seed=…)` returns
  byte-reproducible `cp_quantiles`/`ce_quantiles` and a propagated `band_tier`.
  The **never-synthesize rule** is enforced: a model with no published BSV (Marsh,
  Schnider, the pediatric pair) draws **no band** and says why — a missing band is
  a true statement; a borrowed one is a lie with error bars.
- **Uncertainty-aware divergence (V2).** `compare(..., bands=True, seed=…)` adds two
  machine-readable readouts to `divergence["cp"]`/`["ce"]`: a **separation index**
  (at the instant of peak median spread, are the driver models' bands disjoint? —
  the share of the trajectory where they separate is model-selection risk you cannot
  variability-away) and a **variance decomposition** (structural vs BSV vs residual
  share of the total predictive variance). Models with no BSV are **named** in
  `excluded_from_bands`, never silently dropped. Surfaced in `hypnos compare --bands
  --percentile 5,95 --samples … --seed …`. New `variability.png` figure.
- **Exports (V3, partial).** **NONMEM** now emits the real `$OMEGA` diagonal (with
  `EXP(ETA(.))` wired into `$PK` and CV annotations) and `$SIGMA` from the residual
  model — the single most natural upgrade from v0.1's `$OMEGA 0 FIX` placeholder;
  models with no BSV keep `0 FIX` and name the missing component. **TCI-JSON** carries
  the `variability` block losslessly. Remaining: `$OMEGA BLOCK`, PharmML/nlmixr2/Pumas
  random effects, and BSV backfill for the other models.
- **Safety.** Unchanged and tightened in proportion: bands describe a *given* forward
  dose history and never invert one (no quantile-targeting), and every band is labeled
  a statement about the *model's stated uncertainty*, not a claim about a real patient.
  `clinicalUse = "PROHIBITED"` remains universal. 26 new tests; suite 215 green.

### Inhalational wash-out (FA/FA₀ emergence) — the offset mirror of wash-in (2026-06-09)
- **`hypnos.washout` / `washout_comparison` + `hypnos washout` CLI** — the
  emergence (offset) companion to wash-in, completing the volatile onset→offset
  symmetry exactly as `decrement_time` completed the IV onset/peak/offset triad.
  Same single-compartment alveolar mass balance, run with the inspired fraction
  set to zero and the symmetric early-phase idealization (mixed-venous ≈ the
  initial alveolar fraction, i.e. tissues still saturated), giving
  `FA/FA₀(t) = floor + (1 − floor)·e^(−t/τ)` with `floor = λ·Q̇/(V̇_A + λ·Q̇) = 1 −`
  the wash-in plateau and the *same* τ. New `alveolar_washout` reference kernel.
- **Reproduces the correct clinical ordering from the curated λ alone:**
  desflurane (λ 0.42, floor 0.34) and nitrous oxide (0.47, 0.37) wash out fast,
  isoflurane (1.4, 0.64) slowly — the discriminator is the *floor* (lower = more
  complete, faster emergence), the exact complement of the wash-in plateau. This
  is the physicochemical reason desflurane is preferred for long cases. New
  `washout.png` figure (mirror of `washin.png`); dashboard gains a wash-out chart.
- **Honesty boundary (identical to wash-in):** the math is exact but rests on the
  stated, standard 70-kg-adult ventilation constants and captures only the early,
  lung-dominated phase — the long tissue-release tail (full multi-compartment
  Mapleson) stays out of v0.1 scope. Tests assert FA/FA₀(0)=1, monotone decay, the
  floor bound, floor = 1 − plateau and the shared τ, the one-time-constant 63.2%
  point, and the solubility ordering.

### Dashboard covers the volatiles + export well-formedness hardening (2026-06-09)
- **Dashboard inhalational section:** the Streamlit dashboard showed zero of the 4
  volatile agents (it filtered to compartmental drugs only). It now has an
  "Inhalational agents" panel — age-corrected MAC at the patient's age, the MAC-vs-age
  curve, and the blood:gas-driven wash-in (FA/FI) curves — all computed by the same
  tested package functions (`hypnos.mac`, `alveolar_washin`, `mac_age_corrected`) the
  CLI uses, so the UI never drifts from the data.
- **Export hardening:** the XML well-formedness test now spans *every* PK model, not
  just the three kernel-backed adults — so the kernel-pending exporters (fentanyl,
  rocuronium take a distinct "no instantiated parameters" branch) are covered too.
  `pharmml.py` and `sbml.py` go to 100% line coverage; a malformed-XML or bad-escape
  regression in any model's export can no longer ship.

### Divergence quantification: name the driver + show plasma spread (2026-06-09)
- The divergence metric now names the **driver pair** — the two models furthest
  apart at the instant of peak disagreement — so the view reports not just *how
  much* the models disagree but *which* one is the outlier (e.g. propofol elderly:
  `driver: schnider_1998 vs marsh_1991`). `compare().divergence[...]` gains a
  `driver` field; surfaced in `hypnos compare` and the dashboard.
- **Fixed a real gap:** `hypnos compare` printed only the *effect-site* spread, so a
  comparison whose included models are PK-only (the pediatric Kataria/Paedfusor pair,
  or remifentanil with the PK-only Kim) showed *no* divergence line at all — hiding
  the plasma divergence that is the whole point there. It now reports plasma and
  effect-site spread separately, each with its driver.
- **Housekeeping:** refreshed the stale `hypnos status` example in the README
  (0/15 → 0/19, tier-A models now lead the verify-next list) and completed the
  architecture diagram (it omitted the analysis, inhalational, and verification
  modules and most CLI subcommands). Tests cover the driver and the PK-only plasma
  spread.

### Kataria 1994 pediatric propofol — a live Kataria-vs-Paedfusor divergence (2026-06-09)
- **`hypnotics_iv.propofol.kataria_1994`** executable kernel completes the spec's
  named pediatric pair (Kataria/Paedfusor, spec §1). Weight-proportional
  volumes/clearances with the distinctive age term on V2
  (`V2 = 0.78·WGT + 3.1·AGE − 16`); PK-only (the 1994 disposition study published no
  ke0). Transcribed from the standard published parameter set (same STANPUMP/Shafer
  lineage as Minto — Shafer is a Kataria co-author) and checked to reproduce the
  reference child (23 kg, 7 y: V1 9.43 L, V2 23.64 L). Tier C (single cohort,
  n=53, ages 3–11 y); `review_status` **unverified** pending human PDF confirmation
  of every covariate equation, especially the V2 age term — the field most prone to
  transcription error, and so flagged for the verifier.
- **Pediatric model-divergence is now live:** for a child the view compares Kataria
  vs Paedfusor vs Eleveld on plasma — Kataria and Paedfusor are both in-envelope yet
  peak at 4.88 vs 4.36 µg/mL, making the "Kataria vs Paedfusor" question (a standard
  pediatric-TCI research comparison) measurable. The `pediatric.png` figure
  regenerates to show all three.
- 19 models · 17 executable kernels. New citation `kataria-1994-propofol-pediatric`.
- Alfentanil/fentanyl stay deferred (contested parameterizations); Kataria differs —
  its parameters are consistent across the quality literature, like Paedfusor's.

### Accuracy in the divergence view — connect performance to model-selection (2026-06-09)
- The model-divergence view (`hypnos compare` + dashboard) now reports each
  included model's **published in-envelope MDAPE** next to its curve, so it answers
  both halves of model-selection risk: *how much do the models disagree?* **and**
  *how accurate is each where it's valid?* Previously the curated performance data
  was reachable only via `hypnos performance`, disconnected from the headline view.
- New `Model.predictive_mdape` returns only the in-envelope MDAPE entries: a model's
  out-of-envelope/failure-mode number (e.g. Minto's 53.4% in morbid obesity) applies
  exactly where that model would itself be greyed out, so it is never attached to an
  *included* (in-envelope) model. Multiple in-envelope studies render as a range
  (e.g. Eleveld propofol MDAPE 22–30%); a model with no published MDAPE reads `n/a`.
- Tests assert the in-vs-out-of-envelope split (Minto shows 24.6%, never 53.4%) and
  that both CLI views render the badge.

### Reproducible figures (fix a stale one) + wash-in visualization (2026-06-09)
- **`scripts/regenerate.py` now regenerates every README figure** from the live
  dataset/kernels, in place into `docs/images/`. Previously only `divergence.png`
  was reproducible (and even that was stale — two models, not the committed three);
  `synergy.png`, `pediatric.png`, and `mac_age.png` were orphaned artifacts the
  script could not reproduce, contradicting the spec's "everything is a
  deterministic projection of the dataset" claim.
- **Fixed a factually stale figure + README claim:** `pediatric.png` labelled
  `eleveld_2018 (kernel pending)` and the README said "only Paedfusor is
  in-envelope" — both predate the Eleveld kernel. The broad-envelope Eleveld (Tier
  A) in fact *covers a 6-year-old by design*; the figure and text now show both
  Eleveld and Paedfusor in-envelope, with the adult-only Marsh/Schnider as Tier-D
  pediatric extrapolations. A figure can no longer silently rot out of sync: it is
  redrawn from `compare()`/`simulate()`, so implementing a kernel updates every
  figure that depends on it.
- **New `washin.png`** — the FA/FI wash-in curves + the λ→plateau mechanism, the
  visual for the wash-in feature.
- Smoke test (`test_figures.py`, skips without matplotlib) regenerates all five
  figures to a temp dir and asserts each renders, enforcing the no-rot guarantee.

### Inhalational wash-in (FA/FI uptake) (2026-06-09)
- **`hypnos.washin` / `washin_comparison` + `hypnos washin` CLI** — the volatile
  *uptake* relation the spec names (§3/§6, Phase D) but that was the one Phase-D item
  not yet implemented. A single-compartment alveolar mass balance turns the curated
  blood:gas partition coefficient into the FA/FI wash-in: the early plateau
  `V̇_A/(V̇_A+λ·Q̇)` (the wash-in "knee") and time constant `τ = FRC/(V̇_A+λ·Q̇)`.
- Reproduces the canonical solubility ordering straight from the data — desflurane
  (λ 0.42) and nitrous oxide (0.47) wash in fast, isoflurane (1.4) slowly; the less
  soluble agent has the higher early plateau. New `alveolar_washin` reference kernel.
- **Honesty boundary (same stance as the Greco surface):** the math is exact but
  rests on *stated, standard 70-kg-adult* ventilation constants (V̇_A 4 L/min, FRC
  2.5 L, Q̇ 5 L/min, all overridable). Comparative/education-grade, **not** a
  per-patient FA/FI predictor; full multi-compartment tissue uptake (Mapleson) stays
  out of v0.1 scope. Tests assert FA/FI(0)=0, monotonicity, the plateau bound, the
  one-time-constant 63.2% point, and the solubility ordering.

### Predictive-performance: extend coverage + citation guardrail (2026-06-09)
- **Dexmedetomidine (Hannivoort)** gains its first metrics — best of nine published
  models in a spinal-anesthesia external validation (Obara 2018, *J Anesth*: MDPE
  5.6%, MDAPE 18.1%, wobble 6.2%).
- **Minto remifentanil** now carries an *in-envelope* number alongside its
  out-of-envelope failure-mode number: MDPE −17.3% / MDAPE 24.6% in cardiac surgery
  (Scherrer 2022, *BJA*), against MDAPE 53.4% in morbid obesity (La Colla 2010). The
  record now shows both sides of the envelope at once. Two new citation records.
- **Integrity guardrail:** `hypnos validate` now checks that every
  `predictive_performance` citation resolves to a real citation record (it already
  checked record/parameter/failure-mode citations — performance was the one
  un-guarded path; a typo'd source would previously have passed). Backfilled the
  symmetry across both headline drugs and the α₂-agonist class.

### Predictive-performance backfill + surfacing (2026-06-09)
- **External-validation MDPE/MDAPE backfill** (the named-incomplete Phase E item;
  spec §5 "tier assignment can be partly numeric"). Two new citation records back
  the additions: an independent head-to-head external validation of the propofol
  trio (Hüppe 2020, BJA — Eleveld MDAPE 22%, Marsh 25%, Schnider 26% over 50
  surgical adults; Marsh gains its first metric, Eleveld/Schnider gain an external
  datapoint) and a number that **quantifies a documented failure mode**: the Minto
  remifentanil James-LBM-in-obesity failure now carries the published MDPE −53.4% /
  MDAPE 53.4% (La Colla 2010, *Clin Pharmacokinet*) measured in morbidly obese
  patients — the authors' "not clinically acceptable" made machine-readable.
- **Surfacing:** `hypnos performance [--drug X]` + `hypnos.performance_table(ds)`
  list every metric with its population and resolved DOI (the data was previously
  read only by the BibTeX exporter — invisible to users). `hypnos info` now reports
  `models_with_predictive_performance`. Every row is asserted to resolve to a real
  citation with a DOI — a performance number is never shown bare.
- **Housekeeping:** removed a duplicate `default_schedule_for` import in `cli.py`.

### Offset analysis (decrement time) + solver memoization (2026-06-09)
- **Performance:** `reference.simulate` now memoizes the augmented matrix
  exponential by `(dt, rate)`. On a uniform grid the same propagator recurs every
  step, so this turns O(n) `expm` calls into ~1 — n=8001 dropped 230→41 ms,
  n=361 dropped 68→1 ms; results are bit-identical (it is the same computation,
  cached). The full test suite went from ~26 s back to ~2 s.
- **`hypnos.analysis.decrement_time`** + `hypnos decrement` CLI — the *offset*
  companion to `tpeak` (onset): plasma decline after a constant-*rate* infusion.
  Forward-only; it lengthens with infusion duration for propofol (accumulation,
  1.0→4.2 min over 10→600 min) yet stays near-flat for remifentanil (~2.2 min,
  its celebrated context-insensitivity). Explicitly **not** the classic
  constant-concentration context-sensitive half-time (which needs inverse
  control, out of scope). Hypnos now offers an onset/peak/offset analysis triad.

### Drug-aware dashboard + shared dosing presets + `hypnos models` (2026-06-09)
- **`hypnos.presets`** — drug-appropriate default dose schedules, now a single
  source of truth shared by the CLI and the dashboard (moved out of `cli.py`).
- **Dashboard modernized** (`dashboard/app.py`): a **drug selector** (every
  simulatable drug, not just propofol), drug-appropriate default doses,
  conventional concentration units (ng/mL for opioids), plasma + effect-site
  charts, and an onset (time-to-peak-effect) table. Reuses the tested package
  logic so it never drifts from the CLI.
- **Bug fix:** the dashboard had the same propofol-dose-for-every-drug overdose
  bug previously fixed in the CLI (it hardcoded `drug="propofol"` and a 2 mg/kg
  schedule); both are now fixed at the shared-preset source.
- **`hypnos models [--drug X]`** + `filter.pk_drugs` — list models
  (id/purpose/tier/kernel/review) and discover which drugs are simulatable.

### Time-to-peak-effect (onset) analysis (2026-06-09)
- **`hypnos.analysis.time_to_peak_effect`** + `hypnos tpeak` CLI — the spec's
  `effect_link` "time-to-peak-effect parameterization" (§3). Pure forward
  simulation of a bolus; the effect-site peak is where dCe/dt=0, i.e. Ce=Cp
  (used as an internal sanity check). Dose-independent; envelope-enforced;
  requires a ke0 link (raises for PK-only models like Kim/Paedfusor). Validated:
  Schnider 1.54 min and Minto 1.61 min match the published ~1.6 min; lower ke0
  gives a later peak (Marsh 3.92 > Schnider 1.54 min).
- Documented why **context-sensitive half-time is out of scope**: it requires a
  target-controlled infusion (constant-plasma), i.e. inverse control, which
  Hypnos forbids (§10). The safety boundary shapes which derived metrics exist.

### CSV flat-parameter exporter (2026-06-09)
- **`csv` exporter** (`hypnos.export.csv_flat`) — completes the spec's export
  matrix (§7 "CSV / BibTeX"). One row per parameter across all models, with the
  record tier/review status and the resolved DOI/PMID joined in; proper CSV
  quoting (covariate equations contain commas). Per-model and whole-dataset
  (`hypnos export --format csv` → `parameters.csv`). Wired into the regenerate
  script and CI export matrix.
- Fixed a latent signature bug: `bibtex.build_for_model` (in BUILDERS) did not
  accept the `patient` arg that `export_model` passes; both bibtex and csv
  per-model builders now take an ignored `patient` for a uniform signature.

### Conventional concentration units + drug-appropriate CLI dosing (2026-06-09)
- **Conventional concentration units.** Each drug declares a `concentration_unit`
  (opioids and dexmedetomidine: ng/mL; propofol/rocuronium: µg/mL). The kernels
  still compute one internal unit (µg/mL = mg/L); `SimulationResult`/`Comparison`
  expose `cp_peak_display`/`ce_peak_display` and the unit, and the CLI prints in
  it. Removes the awkward "remifentanil Ce 0.008 µg/mL" readout.
- **Bug fix:** the CLI `simulate`/`compare` previously applied a fixed propofol
  schedule (2 mg/kg) to every drug — a ~1000× overdose for remifentanil
  (26,874 ng/mL artifact). Added drug-appropriate default schedules
  (`_default_schedule_for`): remifentanil mcg/kg + mcg/kg/min, dexmedetomidine
  mcg/kg/h, rocuronium 0.6 mg/kg, etc.

### Kim remifentanil + effect-site divergence fix (2026-06-09)
- **`remifentanil_kim_2017` executable kernel** — completes the spec's named
  remifentanil trio (Minto, Eleveld, Kim). Derived in obesity; Janmahasatian FFM;
  PK-only (the Kim disposition study published no ke0). Transcribed from the
  published equations, cross-checked against the `tci` `pkmod_kim`, validated to
  reproduce the reference individual (37 y, 74.5 kg, FFM 52.3): V1=4.76, V2=8.4,
  V3=4, CL=2.77, Q2=1.94, Q3=0.197. `review_status` unverified. Kim's BMI-to-~70
  envelope covers morbidly-obese patients where Eleveld (BMI 52) and Minto (James
  LBM > 40) are greyed.
- **Bug fix:** `compare()` effect-site (`ce`) divergence now excludes models
  without an effect compartment (`ke0 = 0`). PK-only models (Kim, Paedfusor) have
  `ce = 0` and were spuriously inflating the effect-site spread (e.g. the
  pediatric Paedfusor-vs-Eleveld comparison). Plasma (`cp`) divergence still spans
  every PK model. 18 models, 16 kernels.

### Eleveld two-slope BIS PD model (2026-06-09)
- **`pd_effect.propofol.eleveld_bis`** + `sigmoid_emax_twoslope` kernel — the
  validated PD companion to the Eleveld PK kernel, giving a fully-Eleveld PK-PD
  BIS trajectory. Asymmetric Hill: γ=1.47 below Ce50, γ=1.89 above (steeper);
  Ce50 age-corrected (3.08·exp(-0.00635·(age-35))). Transcribed/validated against
  the `tci` `emax_eleveld`; continuous at Ce50; `review_status` unverified.
- `simulate()` now dispatches the PD kernel by `kernel.function` (single-slope vs
  two-slope with age-adjusted Ce50). Composed Eleveld PK (A) + BIS (B) → Tier B;
  older patients reach deeper BIS at the same dose (lower Ce50). 17 models, 15
  kernels.

### Remifentanil Eleveld 2017 kernel (2026-06-09)
- **`remifentanil_eleveld_2017` executable kernel** — completes the spec's named
  remifentanil pair (Minto + Eleveld). Allometric (fat-free-mass) general-purpose
  model; transcribed from the published equations, cross-checked against the `tci`
  R package, and validated to reproduce the reference individual exactly: V1=5.81,
  V2=8.82, V3=5.03, CL=2.58, Q2=1.72, Q3=0.124, ke0=1.09. `review_status` stays
  `unverified`.
- **Bug found + fixed:** the `tci` source computes V3 from `V2ref` (8.82), which
  does not reproduce the published V3 reference (5.03); Hypnos uses `V3ref`,
  matching the source paper. Documented in the kernel and the record.
- The model-divergence view now works for remifentanil: Minto and Eleveld
  cross-validate for adults (~5% effect-site spread), and Eleveld's broad envelope
  stays valid for the obese/pediatric where Minto is greyed. 16 models, 14 kernels.

### Eleveld propofol kernel (2026-06-08)
- **`propofol_eleveld_2018` executable kernel** — the general-purpose model is no
  longer kernel-pending. Transcribed from the published equations and
  cross-checked against the open `tci` R package (`pkmod_eleveld_ppf`); validated
  to reproduce the reference individual (35 y, 70 kg, 170 cm, male) exactly:
  V1=6.28, V2=25.5, V3=273, CL=1.79, Q3=1.11, ke0=0.146. Arterial-sampling PK
  arm; optional covariates `opiate_coadministration`, `arterial`,
  `postmenstrual_age`; FFM via Al-Sallami 2015.
- `review_status` stays **`unverified`**: an LLM transcription validated to one
  reference point is not a human PDF verification. Run
  `hypnos verify hypnotics_iv.propofol.eleveld_2018` for the checklist.
- The model-divergence view now overlays all three adult propofol models
  (Marsh/Schnider/Eleveld); for the obese patient Eleveld stays in-envelope at
  Tier A while Schnider is greyed out. 13 executable kernels now.

### Verification workflow & docs (2026-06-08)
- **`hypnos.verification`** module + `hypnos status` / `hypnos verify <id>` CLI:
  per-model field-by-field checklists (parameters, covariate equations, envelope,
  derivation population, source DOI), a verification-coverage report, and a
  prioritized "what to verify next" list (implemented kernel + best tier first).
  Guides human verification; never promotes (`LLMs do not promote`).
- **Essay** `docs/about/essay.md` — *Why model-selection risk is the load-bearing
  idea* (the spec's declared conceptual cornerstone).
- **Reference notebook** `notebooks/01_model_divergence.ipynb`, executed in CI via
  `nbmake` so it cannot rot; reproduces the divergence comparison + verification
  workflow. New `notebooks` optional-dependency extra.

### Phase E — Hardening (2026-06-08)
- **COMBINE `.omex` exporter** (`hypnos.export.combine`): bundles SBML (master) +
  PharmML + TCI-JSON + provenance `metadata.rdf` + `citations.bib` with an
  `omex-manifest`. Per-model and whole-dataset archives. Archives are written
  **deterministically** (fixed entry timestamps) → byte-identical for a given
  dataset version.
- **BibTeX citation exporter** (`hypnos.export.bibtex`): per-model and
  whole-dataset citation export.
- CLI: `hypnos export --format omex` (single archive if `--output` ends in
  `.omex`, else per-model) and `--format bibtex`.
- **Predictive-performance backfill** (start): Eleveld MDPE/MDAPE entries
  (older-adult bias −27%; pooled MDAPE < 30%) cited to the source.
- Release/reproducibility: `.zenodo.json`, this `CHANGELOG.md`, and
  `scripts/regenerate.py` (deterministically regenerates all exports + figures).
- CI: builds and validates the `.omex` archive as an artifact.

### Phase D — Inhalational + NMB (2026-06-08)
- Volatiles (`physicochemical` model type): sevoflurane, desflurane, isoflurane,
  nitrous oxide — MAC40 + Mapleson/Nickalls age-correction + partition
  coefficients. `hypnos.mac` API + `hypnos mac` CLI; nitrous-oxide additivity.
- Neuromuscular blockers seeded: rocuronium PK (Wierda, kernel pending) +
  train-of-four PD sigmoid. Sugammadex deferred.

### Phase C — Breadth + pediatrics (2026-06-08)
- Dexmedetomidine (Hannivoort 2015, new `alpha2_agonists` subsystem); pediatric
  propofol (Paedfusor); fentanyl (Shafer, curated, kernel pending).
- Explicit pediatric/geriatric Tier-D extrapolation labeling.
- Loader fix: in-package dataset copy no longer shadows the repo-root
  source-of-truth in source checkouts.

### Phase B — Opioids + interaction (2026-06-08)
- Remifentanil (Minto 1997); propofol–remifentanil Greco response surface
  (`simulate_interaction`); nlmixr2/rxode2 (R) + Pumas (Julia) exporters.

### Phase A — Propofol spine (2026-06-08)
- Dataset (JSON Schema + JSON-LD), propofol PK (Marsh, Schnider; Eleveld curated,
  kernel pending) + ke0 + BIS PD; pure-NumPy/SciPy reference kernels; NONMEM /
  PharmML / SBML / TCI-JSON exporters with round-trip validation; the
  model-divergence view; CLI; CI.
