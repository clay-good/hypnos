# Changelog

All notable changes to Hypnos are documented here. The project follows the
phased roadmap in [`docs/specs/v0.1/spec.md`](docs/specs/v0.1/spec.md) §11.
Dates are ISO-8601. Versioning is [SemVer](https://semver.org/); the dataset
version is pinned in `dataset/VERSION` and stamped into every export as
`hypnos:datasetVersion`.

## [Unreleased]

### Spec completion — every proposed roadmap phase across v0.3–v0.7 is now implemented (2026-06-23)
- **v0.7 C3 — covariate-aware exports.** The derived-covariate equation now rides inside every export:
  NONMEM `$PK` computes the named LBM/FFM equation the model's way (round-trip-equal to the library),
  SBML/PharmML carry it as `hypnos:covariateEquation` RDF, TCI-JSON passes the structured block through
  verbatim. **v0.7 is complete (C0–C3).**
- **v0.4 VE2/VE3 finish — Open-TCI adapter, tier-falsification, VPC.** `subjects_from_opentci` adds the
  gold-standard plasma-concentration (`cp`) adapter; `hypnos validate` raises the advisory tier-mismatch
  flag when an in-envelope MDAPE exceeds a model's tier band; `visual_predictive_check` bins observed vs
  predicted percentiles. **v0.4 is complete (VE0–VE3).**
- **v0.3 E1/E2/E3 — confidence bands.** `simulate(bands=["confidence"])` draws a seeded estimation
  confidence band (Eleveld's is width ~8% vs the prediction band's ~huge — reducible vs irreducible);
  the estimation component folds into `compare`'s variance decomposition; the NONMEM export carries a
  per-θ estimation block. **Eleveld 2018 propofol carries sourced per-θ estimation uncertainty** (an agent
  read the open-access Table 2's 99% profile-likelihood CIs directly; `pending_human_review`). **v0.3 → E0–E3.**
- **v0.5 Part A — the reversal subsystem (B0+B1) + S1 disease modifiers.** Three forward-only kernels
  (`hypnos.reversal`): sugammadex **encapsulation** (rocuronium→TOF recovers ~1%→98% by 30 min), naloxone
  **competitive** antagonism with the **renarcotization** relapse rendered forward, neostigmine **indirect**
  inhibition with the hard deep-block **ceiling** — via `simulate_reversal` + `hypnos reverse`, never a dose.
  The opt-in, Tier-D, cited **`disease_state_modifier`** ships (dexmedetomidine clearance reduced in hepatic
  impairment, De Wolf 2001). The **rocuronium PK kernel** is now implemented from sourced micro-constants
  (`pending_human_review`), unlocking the reversal pairing. **v0.5 → S0/S1/B0/B1 complete.**
- 4 new subsystems' worth of code, 6 new citations (CrossRef-verified), 3 reversal drugs + 3 reversal
  records, the `reversal_mechanism` / `disease_state_modifier` / `source_review` schema blocks. 27 models,
  17 drugs, 38 citations. 499 tests pass; ruff clean; validate clean. The remaining residue is **never-invent
  data backfill** (per-model Ω/Σ and RSE for models other than Eleveld) — the perpetual curation the
  `pending_human_review` mechanism exists to support, not a code gap.

### Curation — the `pending_human_review` state + an automated source cross-check (2026-06-23)
- **A new review state that respects the governance rule while paying down verification debt.**
  Every record was `unverified`, which conflated *illustrative placeholders* with *transcribed-but-
  unconfirmed real values*. A new **`pending_human_review`** state distinguishes them: an automated
  process fetched a real, cited source and numerically compared the curated values — but **no human has
  confirmed against the primary PDF**. The governance line "Humans verify; LLMs do not promote" is
  preserved by construction: an automated check can reach `pending_human_review` but **never** `verified`
  (`source_review.human_verified` is always `false`; `hypnos validate` enforces it).
- **`extraction.source_review`** — structured, auditable provenance for the cross-check: who/what reviewed
  it, the date, the method, the **exact source URLs actually fetched and compared**, the outcome
  (`match` / `partial_match` / `discrepancy`), and notes. A human verifier can audit the machine's work,
  then promote or contest it. Surfaced in `hypnos status` (a new count) and the `hypnos verify` checklist
  (the blocking message names the automated check and what remains), and exported as `hypnos:reviewStatus`.
- **The cross-check, run for real — all 24 models triaged.** Parallel research agents fetched open-access
  reproductions, primary-paper PubMed abstracts (verbatim reference values for Eleveld-remi and Hannivoort),
  the `tci` R package source, opentiva/PyTCI/simtiva/OpenTCI/PAS implementations, FDA/EMA labels (DailyMed,
  SmPC), and CrossRef, then compared each curated parameter. **18 now `pending_human_review`** (Eleveld
  propofol+remi, Minto, Kim, Hannivoort, Schnider, Marsh, Kataria, Paedfusor, Eleveld-BIS, iso/N₂O/sevo MAC,
  bupivacaine + lidocaine disposition, and the source-**filled** fentanyl Shafer-1990 set V1 6.09/V2 28.1/V3
  228 L · Cl 0.504/2.87/1.37 L/min and the source-**filled** rocuronium Wierda set V1 3.08 L / ke0 0.1).
- **Four real problems surfaced, never silently fixed — now `contested`:** (1) the propofol→BIS sigmoid is
  cited to Schnider 1999, which does not model BIS and lacks the curated Ce50 4.7 / γ 1.43 (a Bouillon/Vuyk
  sigmoid); (2) **ropivacaine** disposition is cited to Tucker 1979, but ropivacaine had no human PK until
  Lee 1989 (whose Vss 59 L the curated value matches exactly — GT Tucker is a Lee co-author, explaining the
  drift); (3) **levobupivacaine** is cited to Bardsley 1998, which *explicitly could not compute* Vss/CL —
  the values match the Chirocaine SmPC; (4) the **propofol–remifentanil BIS interaction** is curated as a
  Greco α=4.0 supra-additive surface, but Bouillon 2004's BIS surface is *additive* (Minto-style β=0; synergy
  there applies to laryngoscopy/no-response endpoints, not BIS). Each carries a `source_review` with the
  correct source and a recommended human action.
- **Side-findings + a fixed bug:** the agents caught a parameter-mislabel bug in the upstream `tci` package
  (its Schnider Cl1 mislabels weight as LBM — *our* data is correct), and an **objectively wrong DOI** on the
  Wierda rocuronium citation (`10.1007/BF03007580` resolves to a *scopolamine* paper; **fixed** to
  `…578`, CrossRef-verified).
- **Honest residue:** desflurane MAC40 (curated 6.0 vs sources 6.4–6.6 at age 40) flagged as a `discrepancy`,
  value unchanged; only **2 records stay `unverified`** — succinylcholine PK (the only published set, Hoshi
  1993, is a 1-compartment fit with ~100% CV and a seconds-scale t½, not a usable population model) and the
  rocuronium TOF Ce50/γ (paywalled, unsourceable this session). 12 new tests. Verification coverage went
  from **0 sourced** to **0 verified · 18 pending_human_review · 4 contested · 2 unverified** — the queue a
  human verifier now works top-down.

### v0.4 VE2/VE3 — the cross-model, envelope-stratified external-validation leaderboard, run on real VitalDB (2026-06-23)
- **Tiers that are *earned*, not asserted.** v0.4 VE0/VE1 built the Varvel engine and the VitalDB
  adapter but had only ever run on self-consistency fixtures. VE2/VE3 close the loop: a reproducible,
  **envelope-stratified, cross-model leaderboard** that scores every eligible PK[→PD] stack on **one**
  cohort by the **same** kernels and seed, so the ranking compares *models* — not the cohorts-and-methods
  that the published-number leaderboards conflate. ([Spec §7.1.](docs/specs/v0.4/external_validation.md))
- **`hypnos.cross_model_leaderboard` + `partition_by_envelope` + `hypnos leaderboard`** — the
  adapter-agnostic engine ranks candidates by overall MDAPE, and (VE3) re-scores each model on its
  **in-** and **out-of-envelope** sub-cohorts (`Envelope.check`, the same test that greys a model), turning
  v0.1's *asserted* failure modes into *measured* ones. Deterministic; a kernel-pending or unscorable
  candidate is skipped, never fabricated.
- **The VitalDB harness is now end-to-end** (`scripts/fetch_vitaldb.py`): fetch a cohort from the open
  VitalDB API → **cache the derived per-subject cohort** (gitignored `data/`, raw waveforms never kept) →
  run the cross-model leaderboard → write a **committable, aggregate-only report** (`docs/validation/
  vitaldb_bis_leaderboard.{json,md}`) with a pinned manifest. `--use-cache` re-scores offline, byte-identical.
- **First real run (30 VitalDB cases, measured BIS, seed 7):** all three propofol PK→BIS stacks land at
  **MDPE ≈ −33%** — measured BIS is systematically *lower* (patients *deeper*) than a propofol-only model
  predicts, exactly the remifentanil-synergy bias the adapter caveat names; the PK-model choice barely moves
  the needle (Marsh 32.3% / Schnider 33.7% / Eleveld 34.4% MDAPE, all within ~2%). An honest, interpretable,
  reproducible result — and a textbook demonstration that the *missing synergy term* dominates the *choice of
  propofol PK model* for depth prediction. Every number carries its manifest and its domain-review caveats.
- 5 new tests (per-model envelope partition, self-consistency ranking, determinism, kernel-pending skip,
  record shape). No schema change. **Pending (data, not code):** the Open-TCI plasma-concentration adapter
  (VE for the `cp` target) and a propofol-mono sub-cohort filter to separate the synergy bias from PK error.

### v0.9 — the pharmacogenomic-envelope subsystem (2026-06-23)
- **The genetic special-population axis — and the one distinction it exists to enforce: a susceptibility is
  not a dose change.** v0.9 adds genotype as a *declared* envelope axis with two rigorously-separated record
  kinds. A **`pharmacogenomic_modifier`** (kinetic) scales a forward prediction — opt-in, Tier-D by
  construction, substrate-scoped, always cited. A **`pharmacogenomic_safety_flag`** (susceptibility) forbids a
  trigger (malignant hyperthermia) or warns of a prolonged effect (atypical BCHE) and carries **no kinetic
  effect** — `affects_kinetics` is `False` by construction, the machine-checkable form of the category-error
  guard. The schema `$defs` are disjoint (a modifier cannot carry `trigger_agents`, a flag cannot carry an
  `adjustment`), so MH can never be smuggled in as a clearance change. ([Full design spec.](docs/specs/v0.9/pharmacogenomics.md))
- **`hypnos.pgx_overlay(ds, genotype=...)` + `hypnos pgx` — avoidance first.** For a *declared* genotype the
  overlay leads with the loudest, highest-stakes facts (MH → AVOID volatiles + succinylcholine, rendered as a
  contraindication with no number), then the prolonged-block awareness flags, then the opt-in Tier-D kinetic
  modifiers, then the **substrate guardrails** (the explicitly NOT-affected drugs) and the honest **no-evidence
  majority**. The ordering is the message: most actionable anesthetic pharmacogenetics is avoidance and
  prolonged-effect awareness, not titration.
- **The never-infer rule (new, an ethics line).** A modifier or flag fires **only** on a genotype/phenotype the
  caller explicitly declares — never from ancestry, frequency, or likelihood. Absent a declaration there is no
  genetic effect, shown as "genotype not declared," never "wild-type assumed." Hypnos curates *what a declared
  genotype implies for a model*, never a person's genotype.
- **The substrate guardrail, enforced in code and by the validator.** BCHE acts on succinylcholine, mivacurium,
  and the *ester* local anesthetics only — **not** remifentanil (nonspecific esterases), **not** the amide LAs
  (hepatic), **not** rocuronium (not a cholinesterase substrate). A genetic effect never leaks across a class;
  `hypnos validate` rejects a modifier whose `substrate_scope` excludes its own model's drug.
- **Opt-in forward simulation** (`simulate(..., pharmacogenomics=True)`) applies a declared, substrate-scoped
  kinetic modifier (forcing Tier-D); `evaluate_safety` surfaces a declared safety flag as a contraindication
  warning with no numeric effect. **No genotype-guided dose, ever** (spec §9): a genotype shifts a *forward*
  prediction or raises an avoidance flag; it never yields a dose.
- **Curated** (all `unverified`, cited): the RYR1/CACNA1S→MH avoidance flag (CPIC, Gonsalves 2019) on the
  volatile agents and succinylcholine; the atypical-BCHE prolonged-block awareness flag and the opt-in Tier-D
  hydrolysis-clearance modifier (Lockridge 2015) on a new succinylcholine record (the BCHE-substrate anchor,
  kernel-pending). Three new citations; the `succinylcholine` drug. Exports carry the blocks as `hypnos:`
  RDF (safety flags **never** as a parameter change); every safety-flag-bearing model is `hypnos:safetyCritical`.
- 18 new tests (the flag/modifier separation, never-infer, the substrate guardrail, avoidance-first overlay,
  the no-evidence majority, validator rejections, the non-kinetic RDF). Dataset/schema → `0.9.0`.

### v0.8 — the developmental-pharmacology subsystem (2026-06-23)
- **Making the pediatric/neonatal envelope *speak* mechanistically.** v0.1 could only grey a model below its
  derivation age; v0.8 adds the two mechanisms developmental PK actually turns on: theory-based **allometric
  size** (clearance ∝ weight^¾, volume ∝ weight¹) and the **PMA clearance-maturation** sigmoid (Anderson–Holford,
  `MF(PMA) = PMA^Hill / (TM50^Hill + PMA^Hill)`). A model *fitted in children* (Kataria, Paedfusor, Eleveld's
  broad-application model) is **evidence** — it keeps its tier; an adult model carried into a neonate is an
  explicit, bounded, **Tier-D extrapolation** — never a silent linear-per-kg rescaling. ([Full design spec.](docs/specs/v0.8/developmental.md))
- **`hypnos.developmental` kernels** — `allometric_scale` (¾ on clearances, 1.0 on volumes, relative to a
  reference individual; rate constants untouched), `maturation_factor` (the verbatim Anderson–Holford sigmoid,
  driven by **PMA, never chronological age**), and `apply_developmental` (compose size then maturation onto an
  adult model's volumes/clearances). `linear_per_kg_scale` exposes the naive shortcut **only** as a labeled,
  visibly-wrong reference line.
- **`hypnos.developmental_overlay` + `hypnos developmental` — the headline.** For a neonate/infant the
  eligible-model set shrinks to those with *actual* pediatric standing (greyed + reasoned where out of band);
  every greyed adult model with a curated block fans out as a Tier-D `allometry_plus_maturation` or
  `allometry_only` band (the latter carrying the over-dose caveat: un-modeled maturation OVER-states neonatal
  clearance); the linear-per-kg line is overlaid as the visibly-wrong reference. The honest output is the
  **spread** and the labels — **never a pediatric dose** (spec §9, the most acute inverse-control temptation).
- **Opt-in forward simulation** (`simulate(..., developmental=True)`) carries an adult model into a child by its
  curated block, forcing Tier-D and warning; it refuses (never invents) when no block is curated. New envelope
  dimensions `pma_weeks` / `postnatal_age_days` / `gestational_age_weeks` grey a below-derivation-age patient.
- **The five developmental traps** become validator checks (maturation driven by PMA not chronological age;
  exponents fixed-vs-fitted and on the right parameter; reference weight stated; maturation pathway labeled;
  pediatric-valid size descriptor) and a verification-checklist group. SBML exports carry
  `hypnos:developmentalModel` / `hypnos:maturationFunction` RDF with a comment forbidding chronological-age
  substitution of the PMA driver.
- **Curated** (all `unverified`, Tier-D, cited to Anderson–Holford 2008): Marsh (`allometry_plus_maturation`,
  illustrative TM50 = 47.7 wk / Hill = 3.4 pending source confirmation) and Schnider (`allometry_only`); Eleveld
  is labeled fitted-pediatric evidence (it already covers neonate→elderly, so it earns standing, not an
  annotation). 14 new tests. Dataset/schema → `0.8.0`.

### v0.7 C2 — the covariate-value band + the fifth variance component (2026-06-10)
- **The other half of the covariate uncertainty: how well do we know the covariate?** C1 answered
  *which equation* (a structural choice); C2 answers *how well do we even know the input* — a weight is
  often estimated, a height rounded, and v0.1–v0.3 treat the covariate vector as exact. Supply a covariate
  as a distribution (`weight={"mean":130,"sd":13}`) and `simulate(..., bands=["covariate"], seed=7)`
  propagates it as a **seeded, reproducible covariate band**, framed always as *"how much does
  not-knowing-the-weight move the prediction?"* — never *"what weight should I enter?"*. ([Spec §6/§7.2/§7.3.](docs/specs/v0.7/covariate_uncertainty.md))
- **`covariates.sample_covariate_vector` + `point_patient`** — the seeded primitive draws a perturbed
  covariate vector from `{mean, sd, dist}` marginals (independent, no covariance assumed absent a caller's —
  the honest default of v0.2's `omega_block.complete=false`); scalars are held fixed. A distribution-valued
  covariate is **collapsed to its mean for every deterministic path** (kernel, envelope, PD link), so
  scalar-covariate calls remain byte-identical — a tested backward-compatibility invariant. The
  never-invent rule is load-bearing: **no distribution → no value band** (Hypnos curates no patient's weight).
- **`bands` now accepts band kinds** (`True` still means the v0.2 prediction band; `["prediction","covariate"]`
  opts into both) across `simulate`/`compare`. The covariate band carries its own tier and an honest caveat
  (the dose schedule is held fixed at the point weight, isolating the covariate→PK effect).
- **The fifth variance component.** `compare(..., bands=["prediction","covariate"])` extends the
  reducible/irreducible decomposition: `variance_share` gains `covariate`, and `covariate_split` separates
  its two halves — **equation-choice** (which body-size equation, reducible by *agreeing on the model*) and
  **value-uncertainty** (the input distribution, reducible by *measuring better*) — with a refined
  `reducibility` rollup (`reducible = structural + covariate`). **Gated** on the covariate band being
  requested, so the legacy `bands=True` path keeps its exact 3-way decomposition (tested).
- Reaches the **CLI** (`hypnos simulate --covariate-band --weight-sd 12`; `hypnos compare --bands
  --covariate-band --weight-sd 13`) and the **dashboard** (a weight-uncertainty slider + the covariate-split
  readout). New `docs/images/covariate_band.png` (CI-regenerated): the covariate-value band nested inside the
  wider BSV band, and the decomposition showing — honestly — that for Eleveld the between-subject variability
  (~87%) dwarfs the weight-measurement uncertainty (~2%): *the patient is the uncertainty, not the weight*.
- 13 new tests (covariate band reproducibility, never-invent, widening with SD, point-patient collapse, the
  5-way decomposition + split, legacy-path invariance, the dashboard covariate readout). No dataset/schema
  change (still `0.7.0`). C3 (covariate-aware exports — NONMEM `$PK` computing LBM the model's way) remains.

### v0.7 C1 — the equation-divergence view: model-selection risk inside one model (2026-06-10)
- **Divergence *within* a model, over its covariate equations.** v0.1's `compare` overlays *different
  models* for one patient; C1 adds the orthogonal view — `covariate_divergence(ds, model_id, patient=...)`
  overlays the **same model** under each admissible body-size equation (the LBM↔FFM substitution some
  commercial TCI pumps make silently when they implement "the Schnider model" with Janmahasatian FFM
  instead of James LBM), greying any equation evaluated outside its own validity envelope. It reuses the
  v0.1 `_divergence` machinery, parameterized over equations instead of models, and reports the per-equation
  derived value + curve + envelope status, the peak spread, and the driver pair. ([Spec §7.1.](docs/specs/v0.7/covariate_uncertainty.md))
- **A backward-compatible kernel override.** The five covariate-scaled kernels (Schnider/Minto/Kim/Eleveld
  ×2) gain a one-line `_maybe_override` hook that substitutes a caller-supplied derived covariate when
  `patient["_<quantity>_override"]` is present; absent it, the kernel computes its own verbatim equation, so
  **default behavior is byte-identical** — a tested invariant asserts the verbatim divergence curve
  reproduces a plain `simulate()` exactly.
- **The Schnider-in-obesity failure mode, rendered at its source.** For an obese patient the model's own
  James LBM is greyed (it has inverted), while the FFM substitutions stay in-envelope — a ~9% divergence
  in a single model, from a covariate equation no one was looking at. The substitution is always a labeled
  `verbatim=false` alternative, never presented as the model's own prediction; the inversion is surfaced,
  not silently "fixed" (auto-substituting would erase the documented failure mode — v0.7 §10).
- `hypnos covariate-divergence --model propofol.schnider_1998 --weight 130 ...` (full or shorthand model id),
  the `covariate_divergence`/`CovariateDivergence`/`EquationCurve` API, and `docs/images/covariate_divergence.png`
  (CI-regenerated: the obese-patient curves + a BMI sweep showing the divergence open exactly where James
  leaves its BMI≤37 envelope). 8 new tests. No dataset/schema change (still `0.7.0`). C2 (covariate-value
  bands + the fifth variance component) and C3 (covariate-aware exports) remain.

### v0.7 C0 — the covariate-equation library + bindings (2026-06-10)
- **The fourth hidden uncertainty becomes a first-class object.** v0.1 made model-selection
  uncertainty first-class; v0.2 between-subject; v0.3 estimation. v0.7 C0 curates the uncertainty
  that sits *underneath* all of them — the **covariate transform** every covariate-scaled model
  silently depends on. The choice of, e.g., the James (1976) vs. Janmahasatian (2005) lean-mass
  equation reparameterizes a model with no change to a single published θ, and is the field's
  best-documented "same model, different pump, different concentration" divergence. C0 turns v0.1's
  lone `covariates.lbm_equation` string into a validated, cited, validity-bounded library. ([Full
  design spec.](docs/specs/v0.7/covariate_uncertainty.md))
- **A shared equation library** (`dataset/covariate_equations/`): `james_1976` (LBM), `janmahasatian_2005`
  and `al_sallami_2015` (FFM) — the three equations Hypnos's curated models are actually derived with —
  each carrying its `form`, sex branches, units, **its own `validity_envelope`**, `known_failure_modes`,
  tier, and citation (new: James, Janmahasatian, Al-Sallami). Curated once and shared, exactly like
  `dataset/citations/`. The Du-Bois/Mosteller/Devine equations the spec names are deferred until a
  curated model binds one (no untested, unbound equations — the never-invent rule).
- **`covariate_model` bindings on every covariate-scaled model** (Schnider, Minto → James; Kim →
  Janmahasatian; Eleveld propofol + remifentanil → Al-Sallami): each names the exact published equation
  per derived covariate, what θ it scales (`used_for`), whether it is the authors' `verbatim` choice, with
  its own tier and a `covariate_sensitivity_status` rollup. Schema-additive (`covariate_model` $def +
  `covariate_equation.schema.json`); no v0.1–v0.6 record changes meaning.
- **`hypnos.covariates.evaluate`** generalizes `reference.lbm_james` into a small registered, individually-
  validated library (James/Janmahasatian/Al-Sallami), each a verbatim transcription that checks its own
  validity envelope and **surfaces** the James inversion (never silently "fixes" it). The inversion is a
  *tested property*: James LBM is verified to peak-then-decline above its envelope, and `evaluate` flags it
  (`inverted`, Tier-D, a cited warning). `hypnos covariates` exposes the library on the CLI (list / evaluate
  for a patient / show a model's bindings, surfacing the inversion at its source).
- **`hypnos validate`** gains the covariate layer: every `derived_inputs[].equation` resolves to a library
  record, every `used_for` symbol resolves to a real parameter, the equation's inputs are available to the
  model, `covariate_sensitivity_status` matches contents, every covariate-layer citation resolves, and each
  library record is schema-checked + verified to be implemented in code (data/code stay in sync). The
  verification checklist gains a `covariate_equation` group with the five §4 transcription traps (named
  equation, units, sex coding, equation envelope, fixed-vs-fitted exponents) as explicit human line items.
- New `docs/images/covariate_equations.png` (CI-regenerated from the live library): James LBM peaking then
  inverting vs the monotone FFM equations, and the consequence — Schnider's clearance inflating in obesity
  because its LBM term rests on the inverted equation. 17 new tests. Dataset/schema → `0.7.0`. C0 is the
  self-contained increment the spec blesses; C1 (equation-divergence view), C2 (covariate-value bands + the
  fifth variance component) and C3 (exports compute LBM the model's way) remain.

### v0.4 VE1 — the first external-data adapter (VitalDB) (2026-06-10)
- **The first source-specific adapter for the Varvel engine.** [VitalDB](https://vitaldb.net) is an
  open intra-operative database; its cases carry the TCI pump's propofol delivery and the *measured*
  bispectral index. `analysis.subjects_from_vitaldb` (exported top-level) maps fetched cases to
  `SubjectRecord`s for **PD-BIS validation**: the independent observation is *measured* BIS (not the
  pump's own predicted Ce), and the propofol delivery (`Orchestra/PPF20_RATE`, mL/h, PPF20 = 20 mg/mL)
  is reconstructed as a faithful **step-infusion schedule** (one event per rate change) — exercising
  the multi-event infusion path the engine already supports. A measured-BIS artifact (>100) is dropped;
  a case missing BIS or infusion is dropped (never invents data).
- **`scripts/fetch_vitaldb.py`** — a LOCAL fetch (via the `vitaldb` package + the open clinical table)
  that pulls a cohort, runs `validate_against_cohort`, and writes derived metrics + a manifest to a
  **gitignored `data/`**. Raw records never enter the repo (spec §3); only derived metrics + manifest
  are committable, and only after the maintainer confirms VitalDB's data-use terms. The adapter's
  clinical choices ship with explicit caveats as domain-review items, not verified facts: a
  propofol-only PK→BIS stack does NOT model remifentanil's synergistic BIS deepening (expect systematic
  BIS over-prediction in balanced anaesthesia — a real, interpretable finding), and the track
  names + vial concentration must be confirmed per cohort.
- Tests cover the schedule reconstruction (rate changes → infusion events, mL/h→mg/h), the
  artifact/empty-case dropping, and an end-to-end run through the real Eleveld PK→BIS stack + Varvel
  engine on a synthetic VitalDB-shaped fixture (no network, CI-safe).
- No dataset/schema change (still `0.6.0`). This is the unblocked engineering half of VE1; the
  paper-gated curations (rocuronium PK for reversal; per-θ SEs for E1) remain human-verification work
  by the project's own rule (*LLMs assist, never promote*).

### v0.4 external validation reaches the CLI — `hypnos validate-cohort` (2026-06-10)
- **The Varvel external-validation engine is now a command.** VE0 shipped the metric engine
  (`performance_error` / `varvel_metrics` / `validate_against_cohort`) as a tested *library*, but the
  v0.4 spec frames external validation as a "separately-invoked, locally-run command" (§9) — and that
  surface was missing. `hypnos validate-cohort` exposes the engine: it simulates each subject's
  RECORDED dose history and scores predicted-vs-observed (MDPE / MDAPE / wobble / divergence with a
  seeded-bootstrap CI), forward-only, never tuning a dose.
- **Two modes, both honest.** `--observations <csv>` runs a researcher's own cohort (a long-format
  CSV: `subject,time_min,observed,kind` + covariates + dose specs), mapped to `SubjectRecord`s by a new
  generic `subjects_from_csv` adapter — the source-neutral entry point the Open-TCI / VitalDB adapters
  (VE1) will sit beside. `--self-consistency [--offset X]` needs **no external data**: it builds a
  known-answer fixture from the model's own predictions biased by a known offset, so a +X% offset must
  recover `MDPE ≈ +X%` — a CI-runnable correctness check on the engine (the fixture the spec validates
  the engine with). `--json` emits the schema `external_validation[]` record.
- New `analysis.subjects_from_csv` + `analysis.subjects_from_cohort_self_consistency`, exported from the
  top-level package; tests cover the CSV grouping/parse/skip-malformed, the known-answer recovery
  (offset → MDPE/MDAPE, wobble = 0), and the CLI (both modes + the source-required guard).
- README v0.4 section + cheat sheet + the CLI architecture-diagram node document the command. No
  dataset/schema change (still `0.6.0`); this exposes existing tested machinery, never touches
  credentialed data, and keeps the computed metrics strictly separate from publisher-reported ones.

### v0.5 organ-function reference notebook + architecture-diagram accuracy (2026-06-10) — housekeeping/docs
- **A CI-executed v0.5 reference notebook** (`notebooks/04_organ_function.ipynb`, built
  deterministically by `scripts/build_organ_notebook.py`) reproduces the organ-function envelope —
  the flagship "make the silence speak" feature — from the live kernels: propofol greyed to Tier-D
  with a *named* `HEPATIC EXTRAPOLATION` in cirrhosis, the eligible-model set shrinking under
  `compare()`, remifentanil keeping cited standing (esterase clearance, Dershwitz/Hoke), the
  hypoalbuminemia `BINDING-SENSITIVE` caveat for the highly-albumin-bound drugs, and the proof that
  local anesthetics stay *off* the albumin axis (their binding is α1-AG-driven). Reference notebooks
  now cover v0.1 (divergence), v0.2 (variability), v0.5 (organ envelope), and v0.6 (local anesthetics).
- **Architecture-diagram accuracy.** The README data-flow mermaid diagram omitted the v0.6
  local-anesthetic surfaces; the package node now lists `la` (absorption · double-uncertainty ·
  cardiotoxicity · free-fraction) and the Varvel metrics, the CLI node lists `la`, and the dashboard
  node lists the local-anesthetic panel — so the diagram matches the shipped code.
- No dataset/schema change (still `0.6.0`); documentation + accuracy for the shipped subsystems.

### v0.6 reaches the dashboard + a reference notebook (2026-06-10) — housekeeping/reach
- **The local-anesthetic subsystem now reaches the dashboard.** Every prior subsystem shipped into
  `compare` + the CLI + the Streamlit dashboard (volatiles, the organ-failure overlay, the prediction
  bands all do); v0.6 had reached only the CLI and library. LA models are `purpose: pk` but use the
  Bateman absorption kernel, not an IV-disposition kernel, so `compare()` reports them *unavailable*
  and the divergence view was empty when a local anesthetic was selected. A dedicated **"Local
  anesthetics" dashboard section** (mirroring the self-contained volatiles section) now renders the
  subsystem properly: site-of-injection dominance, the double-uncertainty view (total + non-linear
  free traces against the threshold RANGES, with the dominant-uncertainty readout), and the
  agent-choice cardiotoxicity comparison — all from the same tested `hypnos.la` functions the CLI
  uses, so no surface drifts. A new `test_dashboard.py` case exercises it.
- **A CI-executed v0.6 reference notebook** (`notebooks/03_local_anesthetics.ipynb`, built
  deterministically by `scripts/build_la_notebook.py`) walks all four LA phases from the live kernels —
  site dominance (LA0), the double-uncertainty view (LA1), the stereochemistry/cardiotoxicity margin
  (LA2), and the binding-saturation free-fraction model (LA3) — closing the notebook gap that had stood
  since v0.2 (the prior reference notebooks cover only v0.1 divergence and v0.2 variability). Executed
  in the CI `notebooks` job via `nbmake`, so it cannot rot.
- No dataset/schema change (still `0.6.0`); this is reach + documentation for the completed subsystem.

### Local-anesthetic binding-saturation free-fraction model (v0.6 LA3) (2026-06-10) — v0.6 subsystem COMPLETE
- **The failure mode total concentration hides, made visible.** Local anesthetics are highly,
  *saturably* protein-bound (chiefly to α1-acid glycoprotein); only the FREE fraction is toxic, and as
  total concentration rises the binding saturates so the free fraction climbs non-linearly — total
  concentration under-predicts free-drug risk exactly when risk is highest. LA1 surfaced this as a
  caveat on the linear trace; LA3 now *models* it.
- **`la.saturable_free_concentration` — an exact, reusable capacity-limited (Langmuir 1:1) kernel.**
  A single saturable high-affinity site of capacity `C` and affinity `Ka`, inverted total→free as the
  physical root of a quadratic. The affinity is **pinned by the already-curated low-concentration
  `fraction_bound`** (`C·Ka = fb0/(1−fb0)`), so no second number is invented; the only new curated
  quantity is the binding capacity (≈ AAG molar concentration ~17 µM × drug MW ≈ 4.7–4.9 µg/mL,
  cited). The low-concentration limit reduces exactly to the linear `free = total·(1−fb0)`; as total
  approaches `C` the free fraction climbs from ~5% toward 1.
- **`free_fraction_model` on bupivacaine/levobupivacaine/ropivacaine — Tier-D, illustrative.** The
  QUALITATIVE rise of the free fraction is the documented, load-bearing fact (Tucker & Mather 1979
  concentration-dependent binding); the magnitude is a representative mechanistic illustration, not a
  fitted model, and is labeled Tier-D throughout. Lidocaine (non-saturable) correctly carries no model
  and stays linear; `hypnos validate` rejects a `free_fraction_model` on a non-saturable drug.
- **`la.free_concentration` upgraded.** Where a model is curated it returns the non-linear free trace
  as the primary `c_free` with the linear baseline beside it (`c_free_linear`); where a saturable drug
  has no curated model it stays linear + the LA1 caveat (never a fabricated curve). The
  double-uncertainty view (`hypnos la --site`) now reports the non-linear free peak, the linear
  under-estimate, and the **under-prediction gap** (e.g. bupivacaine 200 mg intercostal: free peak
  1.76× the linear estimate).
- **Surfaced + exported.** `hypnos verify` gains a free-fraction-model checklist item; `hypnos:freeFractionModel`
  rides in SBML/PharmML RDF; TCI-JSON passes the block verbatim (inside `protein_binding`).
- **README figure `la_saturation.png`** (regenerated from the live kernel): the free fraction rising
  with total for the three saturable amides (vs the flat linear assumption), and bupivacaine's free
  concentration with the under-prediction gap shaded.
- With LA3 the **v0.6 local-anesthetic subsystem is complete** (LA0 absorption · LA1 thresholds +
  double-uncertainty · LA2 cardiotoxicity · LA3 saturation). 23 models · 13 drugs · 28 citations.

### Local-anesthetic stereochemistry & cardiotoxicity (v0.6 LA2) (2026-06-10)
- **The agent-choice safety argument, made quantitative and honest.** Racemic bupivacaine is the
  most cardiotoxic common amide ("fast-in/slow-out" avid cardiac Na-channel binding); its
  S-enantiomers (levobupivacaine, ropivacaine) were developed to widen the cardiovascular margin at
  similar local-anesthetic potency. LA2 curates that as a drug-level `cardiotoxicity_class` (`rank`
  ∈ high/intermediate/low, `stereochemistry`, and the qualitative `cns_to_cvs_margin`) for
  lidocaine, bupivacaine, and ropivacaine, and adds **levobupivacaine** (a new drug + systemic
  model, Bardsley 1998) to anchor the intermediate point of the three-way comparison.
- **`la.cardiotoxicity_comparison` / `hypnos la --agents` — the agent-choice axis (v0.6 §3.4/§6).**
  Lists every LA carrying a cardiotoxicity class, most-cardiotoxic first, with its **numeric**
  CNS-to-CVS fold-margin (cardiovascular midpoint / CNS-first-symptoms midpoint, total basis). The
  numeric margin is monotone with the qualitative ranking — racemic bupivacaine **1.4×** <
  levobupivacaine **1.7×** < ropivacaine **2.0×** < lidocaine **5.0×** — so the view shows *why* a
  similar CNS threshold can hide a very different cardiovascular margin, without ranking or
  recommending a dose (forward, comparative, research-only; §7).
- **Consistency fix:** widened ropivacaine's curated cardiovascular range (4–8 → 5–10 µg/mL) so the
  numeric fold-margin order matches the documented qualitative order (ropivacaine's margin is the
  widest of the long-acting amides). The load-bearing claim is the margin *order*, not the precise
  Tier-C numbers.
- **Surfaced everywhere.** The double-uncertainty view (`la.double_uncertainty` / `hypnos la --site`)
  now reports the agent's cardiotoxicity class and raises an explicit **NARROW-margin cardiotoxicity
  warning** for bupivacaine (cardiovascular toxicity can occur with little/no CNS prodrome). `hypnos
  verify` gains a cardiotoxicity checklist item; `hypnos validate` resolves the class citation and
  enforces the rank/margin vocabulary. Export (§9): `hypnos:cardiotoxicityClass` rides in SBML/PharmML
  RDF and TCI-JSON passes the block verbatim.
- **README figure `la_cardiotoxicity.png`** (regenerated from the live kernels): the CNS and
  cardiovascular threshold bands per agent on a shared log axis, the gap widening
  bupivacaine → levobupivacaine → ropivacaine → lidocaine.
- **Design note.** `cardiotoxicity_class` lives at the **drug** level beside `protein_binding` (both
  intrinsic chemistry shared across a drug's models), rather than the model-schema `$def` the v0.6
  spec §4 sketched — matching the established `protein_binding` precedent. 23 models · 13 drugs · 28 citations.

### Local-anesthetic toxicity thresholds + the double-uncertainty view (v0.6 LA1) (2026-06-10)
- **Thresholds are RANGES, never lines — enforced by construction.** A new `toxicity_threshold`
  block (`endpoint` ∈ cns_first_symptoms / cns_seizure / cardiovascular; a `concentration_range`
  with required, non-null `low`/`high`; a mandatory `basis` of total vs free plasma) carries the
  documented systemic-toxicity ranges with their honest, large uncertainty. The schema **forbids a
  single-value threshold** (both bounds required); `hypnos validate` adds the safety-critical
  integrity checks — `low < high`, basis present, only the `local_anesthetics` subsystem may carry
  thresholds, and — the load-bearing guard — a `total_plasma` threshold on a **saturable**-binding
  drug *must* carry a free-fraction `saturation_caveat` (total under-predicts free-drug risk exactly
  when risk is highest). The view's whole design exists to *prevent* the single-line reading (v0.6 §3.3/§4/§7).
- **The headline double-uncertainty view (`la.double_uncertainty`, `hypnos la --site …`).** Overlays
  the predicted total-plasma peak (and the free-concentration peak) against each curated threshold
  *band*, on the matching basis, and names which uncertainty dominates — compared on a consistent
  multiplicative fold-range scale (threshold high/low vs the across-site Cmax spread). For LA the
  threshold band almost always wins, and the readout says so: *"no single safe-concentration line is
  defensible — this is the answer, not a gap to be closed."* It computes **no** dose, ceiling, margin,
  or safe/unsafe verdict (v0.6 §6/§7). The CLI lists the ranges, then renders the full view for a
  chosen `--site`, every line framed RESEARCH/EDUCATION — NOT a dosing tool.
- **The free-concentration trace (`la.free_concentration`).** The bound→free transform applies the
  *linear* binding fraction `c_free = c_total·(1 − fraction_bound)` and, for a saturable drug,
  attaches the saturation failure-mode caveat — the linear free trace **under-predicts** the free
  (toxic) concentration at high total, and that gap is surfaced, never hidden (the nonlinear model is
  LA3). Toxicity tracks *free* drug, so the view shows free beside total wherever it matters.
- **Curated threshold ranges (`unverified`, Tier-C).** Lidocaine CNS-first-symptoms / seizure /
  cardiovascular on a total basis with its comparatively *wide* CNS-to-CVS margin (Tucker & Mather
  1979); bupivacaine CNS (total **and** free) + cardiovascular with the explicit *narrow / absent*
  CNS-to-CVS margin that is its cardiotoxicity story (Knudsen 1997); ropivacaine CNS + cardiovascular
  with the *wider* margin that sets up the LA2 agent-choice comparison (Knudsen 1997, Scott 1989).
  Threshold magnitudes are old, wide, method-dependent and ethically un-doseable in volunteers — Tier-C
  by design; the load-bearing facts are the *margins* and the *uncertainty*, not precise numbers.
- **Export (v0.6 §9): thresholds export AS ranges, and LA records are doubly warned.** Provenance gains
  `hypnos:safetyCritical = "true"` for every LA model (flowing to every text banner too); SBML/PharmML
  carry `hypnos:siteAbsorption` + `hypnos:toxicityThresholdRange` RDF (low/high/basis/units/tier — no
  projection collapses a range to a single value); TCI-JSON passes the `absorption`, `toxicity_thresholds`,
  and drug `protein_binding` blocks verbatim (lossless).
- **Verification: the LA checklist group (v0.6 §8).** `hypnos verify` surfaces the five safety-critical
  human line items — threshold basis (total vs free, Trap 1), salt-vs-base & units (Trap 2),
  speed-of-rise/method context (Trap 3), the saturation caveat for binding-sensitive agents, and that
  every threshold is a *range*, never a line. Given the stakes, threshold ranges require human source
  confirmation before `verified`.
- **README figure `la_double_uncertainty.png`** (regenerated from the live kernels like every Hypnos
  figure): bupivacaine by site against the CNS/cardiovascular threshold *bands*, plus the total-vs-free
  trace showing the binding-saturation gap. Dataset/schema/package version → **0.6.0** (the v0.6 target;
  previously pinned at 0.2.0 across the v0.3–v0.6 work).

### Local-anesthetic systemic-absorption subsystem (v0.6 LA0) (2026-06-11)
- **A new `local_anesthetics` subsystem — and its safety message is its science.** The one
  LA-specific fact that *is* the safety message: **systemic absorption is site-driven, not
  milligram-driven.** The same dose produces wildly different peak plasma concentrations by
  injection site (the documented rank intercostal > caudal/epidural > brachial plexus >
  subcutaneous), so a single mg/kg ceiling is unreliable on its face — and saying so is the
  point. LA0 curates **disposition + site absorption + binding ONLY — no toxicity thresholds**
  (deferred to v0.6 LA1, which needs its own safety framing): the first release teaches that
  the milligram ceiling is the wrong mental model without drawing a single threshold that
  could be misread (v0.6 §10, §7).
- **`hypnos.la` — a forward-only Bateman kernel** (`absorption_pk`: site-set first-order
  absorption into one-compartment disposition, the closed form `C(t) = F·D·ka/(V(ka−k10))·
  (e^−k10·t − e^−ka·t)`), plus `site_comparison` — the headline: the *same* dose of the *same*
  drug at every curated site, sorted by peak. `hypnos la --drug bupivacaine --dose 100` shows
  Cmax falling from intercostal to subcutaneous and peaking progressively later. No dose is
  ever computed (v0.6 §7) — systemic concentration only, **not** block efficacy, **not** a
  toxicity margin.
- **Three LA models** (lidocaine, bupivacaine, ropivacaine) — reference-adult one-compartment
  disposition + a site-absorption block whose `rank` is the robust curated direction even
  where the absolute `ka` is Tier-C (the rank order is established; the magnitudes are old,
  wide, method-dependent — Tier-C by design). All cited to **Tucker & Mather 1979** (the
  canonical LA clinical-PK review), `unverified`. Schema gains an additive `absorption` block.
- **Binding curated, and correctly separated from the v0.5 albumin axis.** LA protein binding
  (lidocaine ~65%, bupivacaine ~95% *saturable*, ropivacaine ~94%) is α1-acid-glycoprotein-
  driven, **not** albumin-driven — so `binding_sensitive` is `false` and the v0.5
  hypoalbuminemia caveat correctly does *not* fire for LAs (their free-fraction story is the
  saturation failure mode of v0.6 LA1/LA3, not the albumin one). `hypnos validate` resolves
  every `absorption`/binding citation and rejects duplicate site ranks.
- **Clean separation from the IV-disposition path.** LA models are `purpose: pk` but their
  kernel is `la_absorption` (not an IV-disposition kernel), so `simulate()`/`compare()` now
  give a clean *"use hypnos.la"* message instead of a `KeyError` (a latent crash this also
  guards for the volatiles' kernel class). 15 new tests.

### The estimation-uncertainty layer + the reducible/irreducible split (v0.3 E0) (2026-06-11)
- **Estimation uncertainty is now a first-class, *separate* layer from between-subject
  variability** — the conflation the field's parameter tables routinely commit (an RSE
  printed in the column beside a BSV CV, indistinguishable but physically different). v0.3
  gives the dataset the vocabulary to record the two apart, **closing the conflation by
  construction**: a new per-θ `estimation_uncertainty` block (`se`, `scale`, `rse_percent`,
  `ci95`, `method`, its own tier) lives **beside** `variability`, never inside it, so an
  estimation RSE can never be silently read as a between-subject CV. Plus a model-root
  `estimate_covariance` (the NONMEM `$COV` output, with `covariance_step_succeeded` recorded
  honestly), `uncertainty_status` (`none | marginal | correlated`), and `omega2_se` (the
  second-order SE on a BSV variance). Schema is additive; no existing record changes meaning.
- **`hypnos validate` enforces the numeric traps a machine *can* catch** (v0.3 §4): `scale`
  is mandatory when an SE is present (Trap 2 — a log-scale SE applied as natural is silently
  wrong); `rse_percent` recomputes from `se` and `value.central` on the declared scale; `ci95`
  is consistent with `se` for a symmetric asymptotic interval (Trap 3); every estimation
  citation resolves; and `uncertainty_status` matches the curated contents (no `correlated`
  without an `estimate_covariance`). The cardinal RSE-vs-CV disambiguation (Trap 1) stays a
  **human** checklist line item — both magnitudes are plausible, so the machine cannot guess;
  the schema separation is the structural guard. `hypnos verify` gains an **estimation** group.
- **The reducible/irreducible decomposition — the v0.3 headline (§7) — lands now, number-free.**
  `compare(..., bands=True)` already split the predictive variance into structural / BSV /
  residual; it now adds a `reducibility` rollup that reframes the same variance along the axis
  that decides *what to do about it*: **reducible** (between-model — curate/validate more models;
  and estimation — more data per model) vs **irreducible** (BSV + residual — the population is the
  limit, the assay is noisy). The punch-line a single curve can never give: *how much of this
  uncertainty could research buy down, and how much is the patient?* The estimation (more-data)
  component contributes 0 until per-θ SEs are curated, and the readout says so honestly
  (`estimation_curated: false`) rather than overstating what is reducible.
- **Why the Eleveld estimation values are not curated here.** E0's headline curation — each
  model's per-parameter RSE table — requires transcribing it from the source PDF, exactly the
  human-verification step the project refuses to let an LLM perform (and the Eleveld 2018 full
  text is paywalled to automated fetch). So this ships the **machinery + enforcement** (the
  spec's stated E0 value: *"closing, by construction, the conflation"*) and leaves the RSE
  *curation* — and the confidence bands (E1) it would unlock — as the explicit human-PDF step.
- 16 new tests (rse↔se on natural/log scales, the four traps on constructed records,
  parsing/tiers, estimate-covariance, and the reducibility rollup on real `compare` output).
- **Housekeeping:** removed two pre-existing unused `typing` imports (`csv_flat.py`,
  `inhalational.py`); `ruff check python/hypnos/` is now clean.

### The protein-binding / free-fraction failure mode (v0.5 S1) (2026-06-11)
- **The albumin axis of the v0.5 organ envelope now says *why* it matters, with a citation.**
  S0 greyed every model under hypoalbuminemia with a generic message. For a **binding-sensitive**
  drug, the view now adds the specific, safety-relevant failure mode: a highly protein-bound drug's
  **free (active) fraction rises** when albumin falls, so a model fit in normal-albumin patients
  **under-estimates** effect from a given *total* concentration (v0.5 §B3/§B7). The free-fraction
  shift is **surfaced, never silently modeled** — the never-invent rule: Hypnos flags the documented
  failure mode rather than fabricating a free-fraction correction it cannot source.
- **Drug-level `protein_binding` curated for the three binding-sensitive anesthetics** — propofol
  (~98% bound; **Zamacona 1997**, *Acta Anaesthesiol Scand* 41:1267-72, the critically-ill
  protein-binding study), fentanyl (~84%; **Meuldermans 1982**, *Arch Int Pharmacodyn Ther* 257:4-19),
  and dexmedetomidine (~94%; **Weerink 2017**, *Clin Pharmacokinet* 56:893-913) — each
  `{fraction_bound, binding_sensitive, citation, note}`, three new citation records. Remifentanil
  (~70% bound) is deliberately **not** flagged — no claim without standing.
- **`evaluate_safety(model, patient, drug_meta=None)`** gains the optional drug record; when a
  hypoalbuminemia organ-finding fires and the drug is binding-sensitive, it appends the cited
  `BINDING-SENSITIVE:` caveat. `simulate`/`compare` thread the drug record automatically, so the
  caveat appears in the CLI and the dashboard's organ-failure overlay with no UI change. Backward
  compatible — callers that pass no `drug_meta` (and normal-albumin patients) are unaffected.
- **`hypnos validate`** now resolves every drug-level `protein_binding` citation (the same
  never-assert-bare rule the rest of the dataset follows).
- 12 new tests (curated chemistry, the caveat firing only for a binding-sensitive drug + low
  albumin, backward-compatible defaults, the end-to-end `simulate` path, citation resolution).
- **Fix:** the per-model BibTeX export now collects the **organ-tolerance** (model-level) and
  **protein-binding** (drug-level) citations too — `hypnos export --format bibtex` was silently
  dropping the 5 v0.5 references (Dershwitz/Hoke/Zamacona/Meuldermans/Weerink), emitting 19 of 24.
  The dataset-level `.omex`/combine export already carried all of them; only the per-model
  ref-collection lagged. Regression test added.
- **Housekeeping:** removed a pre-existing unused `load` import in `simulate.py` (flagged by ruff;
  it lived only in a docstring word).

### The organ-function envelope — making Hypnos speak on organ failure (v0.5 S0) (2026-06-10)
- **The *physiological* envelope now speaks.** v0.1's envelope greys a model when a
  patient's *demographic* covariates (age/weight/BMI) fall outside the derivation range;
  it was silent on the highest-stakes extrapolations — **hepatic, renal, cardiac, and
  protein-binding (albumin) impairment** — and silence reads as "fine." A patient who
  declares organ impairment now greys every model with no cited standing in that state
  (Tier-D + a **named** extrapolation: `HEPATIC EXTRAPOLATION: Child-Pugh C … -> Tier D`),
  exactly as a too-high BMI already does. The eligible-model set visibly **shrinks** under
  organ failure (v0.5 §B6) instead of pretending nothing changed.
- **Remifentanil keeps its standing — with a citation.** The one well-established exception
  is encoded, not hand-waved: remifentanil is cleared by nonspecific blood/tissue esterases,
  so its clearance is essentially **independent of hepatic and renal function**. The three
  remifentanil models (Minto/Eleveld/Kim) carry a cited `organ_tolerance` on the hepatic and
  renal axes (**Dershwitz 1996**, *Anesthesiology* 84:812-20; **Hoke 1997**, *Anesthesiology*
  87:533-41 — two new citation records). A cirrhotic/renal-failure patient therefore greys
  every propofol model but leaves remifentanil standing, with the honest **caveat** that the
  active acid metabolite (GR90291) accumulates in renal failure (the parent PK these models
  predict is unaffected). Cardiac and albumin axes still grey remifentanil — no model was
  fit there, and the never-invent rule holds: no standing is claimed without a source.
- **Schema (additive, backward-compatible):** `applicability_envelope` gains optional
  numeric organ-function ranges (`crcl_ml_min`, `albumin_g_dl`, `ejection_fraction_pct` —
  the forward-compatible "fitted-in-disease" slot) and an `organ_tolerance[]` block (cited
  mechanistic standing). No existing record changes meaning. The impairment cut-points
  (CrCl<60 = KDIGO CKD stage ≥3; albumin<3.5; EF<40; any Child-Pugh class = chronic liver
  disease) are **definitional clinical staging**, not fitted PK, so they live named in code —
  not curated per model.
- **Surfaced everywhere the demographic envelope is:** `evaluate_safety`/`Envelope.organ_check`
  drive `simulate` and `compare`; the **CLI** gains `--child-pugh / --crcl / --albumin /
  --ejection-fraction`; the **dashboard** gains an "Organ function" panel that renders the
  envelope-shrinkage overlay live. `hypnos validate` resolves every `organ_tolerance` citation
  and well-orders the new ranges. A normal simulation (no organ covariates) is **unaffected**.
- 23 new tests (organ-envelope behavior, exact staging cut-points, standing-via-fitted-range,
  the compare overlay incl. the graceful all-greyed case, citation-resolution, dashboard overlay).

### External-validation metric engine — Varvel's framework (v0.4 VE0) (2026-06-10)
- **New `hypnos.analysis` metric engine** implementing the field-standard anesthesia
  PK/PD validation methodology (Varvel 1992): `performance_error(c_obs, c_pred)` (the
  signed PE%, with the non-positive-prediction guard the metric requires) and
  `varvel_metrics(pe, times)` → **MDPE** (bias), **MDAPE** (inaccuracy), **wobble**
  (intra-individual variability), and **divergence** (the %/h drift of `|PE|` with time).
  `pooled_performance(...)` rolls subject-level metrics up to a population median with a
  **seeded** nonparametric bootstrap 95% CI — identical `(subjects, seed)` → identical CIs
  (v0.4 §3 determinism), with all-nan metric columns handled quietly (never imputed).
- **`validate_against_cohort(ds, model_id, subjects, target=...)`** — the adapter-agnostic
  harness: it drives the *existing* forward solver from each subject's **recorded** dose
  history + covariates, interpolates the prediction to the observation timestamps, and
  computes the metrics. Forward-only — it runs the recorded doses and never searches for a
  dose (v0.4 §10). Source-specific adapters (Open-TCI, VitalDB) sit above this and are
  deferred to VE1/VE2 with their data-manifest + ethics handling (v0.4 §3); the engine
  itself needs no credentialed data, so CI stays green on synthetic fixtures.
- **Schema (additive, backward-compatible):** new model-root `external_validation[]`
  (Hypnos-**computed** Varvel metric sets — reproducible, kept strictly separate from the
  human-curated publisher-reported `predictive_performance`, per v0.4 §4.1) and
  `validation_status` rollup (`none | internal_only | external_pk | external_pd |
  external_both`), plus `external_validation_entry` / `validation_metric` `$defs`. No
  existing record changes meaning; no model carries the block yet (the metrics are a
  *generated* artifact, produced only by running the engine against real data under its
  access terms). `CohortValidation.to_record()` serializes a run into a schema-valid entry.
- **`hypnos validate` gains dormant-until-populated consistency checks** (v0.4 §4): a
  computed block can never be mislabeled — `validation_status` must match the curated
  entries, each entry's `target` must match its `mode` (a BIS validation can't be filed as
  a concentration validation), and metric CIs must be well-ordered.
- **Tested against hand-computed answers** on synthetic fixtures (the textbook edge cases:
  a single sample, all-zero error, monotone drift for divergence, a non-positive
  prediction) **and** end-to-end against the real Eleveld kernel via a self-consistency
  fixture (a model's own prediction fed back as the observations scores ~0 error; a uniform
  +20% offset scores MDPE = MDAPE = 20% exactly) — proving the alignment + solver wiring
  without curating or inventing a single clinical concentration. 23 new tests.
- *Why VE0 (v0.4) landed ahead of the v0.3 estimation-uncertainty curation:* VE0 is pure
  algorithm and needs **no** new curated numbers, whereas v0.3 E0's headline requires
  transcribing each model's per-parameter RSE table from the source PDF — exactly the
  human-verification step the project refuses to let an LLM perform. Building the
  fabrication-free engine first, and leaving the RSE *curation* as the explicit
  human-PDF gap, is the never-invent ethos applied to the roadmap itself.

### Verification checklist surfaces (and flags missing) per-parameter source locators (2026-06-10)
- **`hypnos verify <id>` now shows each structural parameter's curated `source_locator`**
  (`@ Schnider 1998, Table 2`) — provenance that was already in the dataset but hidden
  from the one workflow that needs it most, the human PDF verification that is the
  project's single highest-leverage contribution. A verifier now goes straight to the
  right table instead of hunting the whole PDF (it even surfaces provenance nuances like
  Schnider's ke0 coming from the 1999 PD-companion paper).
- **Where a structural parameter has no curated locator, the checklist flags the gap**
  (plain text `[! no source locator curated — add one]` + a summary nudge; markdown
  `⚠️ no source locator curated`), so verification *adds* provenance over time. Audit
  found 35 of 99 parameters lack a locator today (concentrated in Marsh, Paedfusor, the
  Greco surface, the PD sigmoids, and the volatiles) — now visible instead of silent.
  The dataset stays honest about its own provenance, not just its parameter values.
- New `ChecklistItem.locator` field; rendered by both the CLI and the markdown PR
  checklist. Tests cover the present-and-surfaced and absent-and-flagged cases.

### Dashboard surfaces the PD effect (BIS) prediction band (2026-06-10)
- **The dashboard now renders the v0.2 effect band**, closing the same "feature in the
  package but not in the UI" gap a prior commit closed for the volatiles. With the
  prediction-band toggle on, it composes the band-eligible PK model with its matching
  PD model (Eleveld PK + Eleveld BIS for propofol, discovered generically by drug +
  shared author token) and draws the BIS effect ribbon — the curated PK between-subject
  variability pushed through the (fixed) Hill link, quantiles taken on the effect draws
  directly. Labeled an honest **lower bound** on true effect spread, because PD-parameter
  BSV (Ce50, γ) is not curated (never-invent, carried into effect space). All compute is
  the tested `simulate(..., pd_model=…, bands=True)` path, so the view never drifts from
  the CLI/figure. The `AppTest` smoke test now asserts the panel renders.

### Fix the `compare --bands` band readout + add the v0.2 reference notebook (2026-06-09)
- **Fixed a band-display bug in `hypnos compare --bands`.** The per-model band was
  printed as `[max(q_lo), max(q_hi)]` — the *independent temporal maxima* of the 5th-
  and 95th-percentile curves, which peak at different instants, producing an incoherent
  interval (e.g. Eleveld's elderly-patient Ce band read `[1.20, 8.20]`). It now reports
  the band at the **median-peak instant** `i = argmax(q[50])` → `[q_lo[i], q_hi[i]]`,
  matching `hypnos simulate` and the dashboard ribbon exactly (`[0.87, 6.90]`). The
  underlying quantile arrays were always correct; only this one CLI summary line was
  wrong. New regression test (`test_compare_band_matches_simulate_at_peak`) asserts the
  two views agree. Stale README band number refreshed.
- **New reference notebook** `notebooks/02_population_variability.ipynb`, executed in CI
  via `nbmake` (the v0.2 sibling of the divergence notebook). It reproduces the seeded
  prediction bands, the never-synthesize rule (no band for a no-BSV model), seed
  determinism (and the refusal of unseeded bands), the variance decomposition, and the
  PD effect band — entirely from curated data, with assertions that lock the behavior.

### v0.2 — the population-variability layer (V0–V3 export code complete) (2026-06-09)
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
- **PD effect bands (§14).** `simulate(..., bands=True, pd_model=…)` now propagates the
  curated **PK** between-subject variability through the (deterministic) PD link to an
  **effect band**: each virtual individual's true effect-site curve is pushed through
  the Hill model and quantiles are taken on the effect draws directly (correct under the
  monotone, non-linear transform — the population-median BIS differs from the typical-
  individual line, as it should). New `result.effect_quantiles`; surfaced in the CLI
  (`hypnos simulate … --pd … --bands`, which also gains `--seed/--percentile/--samples`)
  and drawn in `docs/images/effect_band.png`. PD-parameter BSV (Ce50, γ) is **not**
  curated, so the effect band is labeled an honest **lower bound** on true effect spread
  — the never-invent rule, carried into effect space. No band is drawn for a no-BSV PK
  model (never-synthesize) even with a PD model attached.
- **Dashboard ribbons (V2).** The Streamlit dashboard now renders the prediction
  bands as shaded **Altair ribbons** — a *Seeded 5–95% prediction bands* toggle (with
  seed + Monte-Carlo-samples controls) overlays each band-eligible model's percentile
  ribbon on its median line; no-BSV models stay bare dashed lines and are named, never
  given a fabricated band. Beside the chart: the **separation index** (are the bands
  distinguishable?), the **structural/BSV/residual variance decomposition**, and the
  `excluded_from_bands` models. All from the same seeded `compare(..., bands=True)`
  call, so the view never drifts from the CLI. Smoke-tested end to end via Streamlit
  `AppTest` (skipped where the dashboard extra is absent, e.g. the CI `test` job).
- **Exports (V3 — the random-effects layer now round-trips through *every* population
  format).** A shared `export/_variability.py` projects the curated Ω/Σ once; each
  exporter renders it in its own idiom:
  - **NONMEM** emits the real `$OMEGA` diagonal (with `EXP(ETA(.))` wired into `$PK`
    and CV annotations) and `$SIGMA`, and now a real **`$OMEGA BLOCK`** when a
    `complete` `omega_block` spans a contiguous, front-anchored η run (covariance
    `r·√(ωᵢ²ωⱼ²)`); non-block/incomplete correlations fall back to a diagonal `$OMEGA`
    plus a named caveat; no-BSV models keep `0 FIX` and name the missing component.
  - **PharmML** gains a first-class `VariabilityModel`: each η → `RandomEffect`/
    `VariabilityLevel` (log transform, variance + CV%), correlations → `Correlation`,
    ε → `ResidualError`. The interop anchor now carries the whole NLME object.
  - **nlmixr2/rxode2** gains a runnable `<id>_pop` companion (V/Cl re-parameterized
    with `exp(eta.X)`, micro-constants derived, `cp ~ lnorm/prop/add` residual) plus
    the Ω as a `lotri` matrix — round-trips to the same typical micro-constants at η=0.
  - **Pumas** gains a `<id>_pop` `@model` with `@param` ω² + `@random` η, the V/Cl
    `@pre` with `exp(η)`, and the residual `dv ~` distribution in `@derived`.
  - **TCI-JSON** carries the `variability` block losslessly (unchanged).
  - **SBML** — whose L3v2 core is deterministic and cannot express population random
    effects — now carries the curated Ω diagonal, Σ, and any off-diagonal correlations
    as `hypnos:` RDF predicates inside the model annotation (`hypnos:betweenSubjectVariability`
    / `hypnos:residualError` / `hypnos:omegaCorrelation`), with an inline SBML comment that
    a downstream deterministic consumer (COPASI/Tellurium) sees only the typical patient.
    The Ω/Σ ride in the annotation, not the `<parameter>` block, so the v0.1 micro-constant
    round-trip is untouched. With this, **every** Hypnos export carries the random-effects
    layer — closing the §8 matrix (SBML was the last "deferred" row).
  Remaining for V3 is **data, not code**: BSV backfill for Schnider, the remifentanil
  trio, and dexmedetomidine awaits source-table confirmation (the common-toolchain
  Schnider CV%s look like RSEs, not BSV), per the never-invent rule — the exporters
  above will carry it the moment it is curated.
- **Safety.** Unchanged and tightened in proportion: bands describe a *given* forward
  dose history and never invert one (no quantile-targeting), and every band is labeled
  a statement about the *model's stated uncertainty*, not a claim about a real patient.
  `clinicalUse = "PROHIBITED"` remains universal. 47 variability tests; CI suite 236
  green (+2 dashboard `AppTest` smoke tests that run wherever the dashboard extra is
  installed, skipped in the dev-only CI `test` job).

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
