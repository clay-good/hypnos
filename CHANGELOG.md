# Changelog

All notable changes to Hypnos are documented here. The project follows the
phased roadmap in [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) §11.
Dates are ISO-8601. Versioning is [SemVer](https://semver.org/); the dataset
version is pinned in `dataset/VERSION` and stamped into every export as
`hypnos:datasetVersion`.

## [Unreleased]

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
