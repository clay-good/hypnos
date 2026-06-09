# Changelog

All notable changes to Hypnos are documented here. The project follows the
phased roadmap in [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) §11.
Dates are ISO-8601. Versioning is [SemVer](https://semver.org/); the dataset
version is pinned in `dataset/VERSION` and stamped into every export as
`hypnos:datasetVersion`.

## [Unreleased]

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
