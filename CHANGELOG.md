# Changelog

All notable changes to Hypnos are documented here. The project follows the
phased roadmap in [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) §11.
Dates are ISO-8601. Versioning is [SemVer](https://semver.org/); the dataset
version is pinned in `dataset/VERSION` and stamped into every export as
`hypnos:datasetVersion`.

## [Unreleased]

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
