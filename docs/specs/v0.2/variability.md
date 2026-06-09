# Hypnos — design spec (v0.2): the population-variability layer

**Curate the *random-effects* structure of each population PK/PD model — between-subject variability (η / Ω), residual error (ε / Σ), and where published, inter-occasion variability — alongside the typical-value parameters Hypnos already curates; propagate it through forward simulation as seeded Monte-Carlo prediction bands; and upgrade the model-divergence view from "the models' point estimates differ by *X*" to the sharper, more honest question: *do the models disagree beyond their own stated uncertainty?***

> **Status:** **V0–V3 implemented (export code complete); only the never-invent BSV data backfill remains** (dataset/schema `0.2.0`). Extends [v0.1](../v0.1/spec.md) (which curates the fixed-effect / typical-value layer). Additive only — every v0.1 record stays valid; the variability block is optional and defaults to "not curated," which renders as an explicit gap, never as a fabricated band. Eleveld propofol now carries a verified-pathway (`unverified`) Ω-diagonal + Σ; `hypnos validate` enforces the §4 traps; `simulate`/`compare` draw seeded bands and the uncertainty-aware divergence; and **all five population exports** (NONMEM, PharmML, nlmixr2/rxode2, Pumas, TCI-JSON) carry the random-effects layer, including off-diagonal `$OMEGA BLOCK` where a complete correlated block is curated. See the roadmap (§12) for per-phase status and what remains (Schnider/opioid/dexmedetomidine BSV backfill — blocked on source-table confirmation, not on code).

> This spec inherits v0.1's stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input.* The variability layer is that stance taken one level deeper. v0.1 made **model-selection** uncertainty a first-class, machine-readable field. v0.2 makes **parametric / between-individual** uncertainty a first-class field too — and, crucially, makes the *relationship between the two* measurable.

---

## 1. The problem this layer solves — the second hidden uncertainty

A v0.1 forward simulation answers: *for the typical patient matching these covariates, model M predicts this concentration-time curve.* That single curve hides **two** distinct uncertainties, not one:

1. **Structural / model-selection uncertainty** — *which model did you pick?* Marsh, Schnider, and Eleveld disagree for the same patient and dose. v0.1 made this visible: the divergence view overlays every eligible model and reports the spread.
2. **Parametric / between-subject uncertainty** — *even within one model, real individuals scatter around the typical-value prediction.* This is the entire reason population models are fitted as **non-linear mixed-effects (NLME)** models: a fixed-effect (typical) layer *plus* a random-effect layer (Ω for between-subject variability, Σ for residual error). v0.1 curates the first layer and silently drops the second.

Dropping the second layer is not a rounding error — it is the difference between a line and a distribution. Two failure modes follow, and they are opposite:

- **False alarm.** Two models' typical-value curves differ by 40% at the peak, so the divergence view flags large model-selection risk. But each model's published between-subject variability is ±50%, and their prediction bands overlap almost completely. The "disagreement" is *inside the noise both models already declare*: choosing between them changes little that their own uncertainty doesn't already swamp.
- **Missed alarm (the dangerous one).** Two models' typical-value curves differ by only 15%, so the divergence view reads as "models broadly agree." But both models report *tight* between-subject variability, their bands are narrow, and the bands **do not overlap**. That 15% is a genuine, irreducible, structural disagreement — neither model's own stated uncertainty can explain it away. This is exactly the case where model selection matters most, and a point-estimate-only view under-weights it.

> The field already separates these ideas in its vocabulary — Varvel's **inaccuracy** (MDAPE) and **divergence** are between-model/between-prediction, while **wobble** (intra-individual variability of the performance error) is within-subject. Hypnos v0.1 curates MDAPE and divergence; **wobble and the Ω/Σ that generate it are the missing half.** Curating them lets the dataset say, per patient and per instant, *how much of the total predictive uncertainty is "which model" versus "which patient" versus "measurement noise."*

The headline deliverable (§7) is the machine-readable answer to: **are these models actually distinguishable for this patient, or is their disagreement within the variability each one already admits?** No open resource answers this today.

---

## 2. What we curate — the NLME variability structure

Hypnos already treats a model record as a structured object, not a scalar (v0.1 §4). The variability layer adds three components of the standard NLME description. We curate them in the **canonical, unambiguous form** and let exporters project to each target's convention.

### 2.1 Between-subject variability (BSV) — the η / Ω layer

The near-universal convention for disposition parameters is **log-normal (exponential) BSV**:

```
P_i = P_typical(θ, covariates_i) · exp(η_i),   η_i ~ N(0, ω²)
```

so each structural parameter (CL, V1, Q2, …) may carry a between-subject variance **ω²** on the η (log) scale. Collected across parameters these are the diagonal of the variance-covariance matrix **Ω**. Where a model fits **off-diagonal** terms — correlated random effects, e.g. a CL–V1 correlation — we curate the correlation/covariance too (§3.3); where it does not, we record that the off-diagonals were **not published** (an honest gap, not an assumed zero — see §5).

### 2.2 Residual unexplained variability (RUV) — the ε / Σ layer

The observation model maps a prediction to an observation with residual error, most commonly **combined proportional + additive**:

```
C_obs = C_pred · (1 + ε_prop) + ε_add,   ε_prop ~ N(0, σ²_prop),  ε_add ~ N(0, σ²_add)
```

with proportional-only and additive-only as special cases, and the log-error form for some assays. Σ is the residual layer. It is conceptually distinct from BSV: BSV is *real between-person differences in physiology*; RUV is *assay noise + model misspecification + sampling-time error*. A prediction band that includes Σ answers "where will a measured sample land?"; one that includes only Ω answers "where will a person's true concentration be?" Hypnos curates both and keeps them separable, because the right one depends on the question.

### 2.3 Inter-occasion variability (IOV) — curated when published, never invented

Some models (notably newer general-purpose ones) fit **inter-occasion variability** — within-subject, between-occasion random effects. Where a source publishes IOV we curate it as a labeled, separable component; where it does not, the field is absent. IOV is **never** folded silently into BSV (a frequent transcription error that inflates apparent between-subject spread).

---

## 3. Schema extension

Additive and backward-compatible. Every addition is an *optional* property; `additionalProperties: false` is preserved by enumerating the new keys. No v0.1 record changes meaning.

### 3.1 Parameter-level: `variability` (BSV per structural parameter)

Added to each entry of `parameters[]`:

```jsonc
{
  "symbol": "Cl",
  "value": { "central": 1.79, "units": "L/min" },
  "covariate_model": "...",
  "tier": "A",
  "extraction": { "review_status": "unverified" },
  "primary_citation": "eleveld-2018-propofol",

  "variability": {
    "bsv": {
      "omega2": 0.0625,            // CANONICAL: variance on the η (log) scale; the NONMEM $OMEGA diagonal entry
      "cv_percent": 25.0,          // OPTIONAL human-readable convenience; relation made explicit (§4)
      "shrinkage_percent": null    // η-shrinkage if the source reports it (interpretation caveat, §4)
    },
    "iov": { "omega2": null },     // inter-occasion variance, only if separately published
    "tier": "B",                   // variability has its OWN tier — usually <= the point-estimate tier (§5)
    "extraction": {
      "review_status": "unverified",
      "tier_rationale": "BSV transcription unverified; confirm scale = log-variance, not SD or CV%, against Table N.",
      "source_locator": "Eleveld 2018, Table 3 (Omega)"
    },
    "primary_citation": "eleveld-2018-propofol"
  }
}
```

`omega2` is the single source of truth (an η-scale variance, exactly the value a NONMEM `$OMEGA` diagonal carries). `cv_percent` is a derived convenience whose exact relation to `omega2` is fixed in §4 so the two can never silently disagree — `hypnos validate` checks their consistency.

### 3.2 Model-level: `residual_error` (the Σ layer)

```jsonc
"residual_error": {
  "model": "combined",                         // proportional | additive | combined | log
  "proportional": { "variance": 0.04 },        // σ²_prop on the multiplicative scale
  "additive": { "sd": 0.05, "units": "ug/mL" },// σ_add in concentration units (note: SD, labeled)
  "tier": "B",
  "extraction": { "review_status": "unverified", "source_locator": "Eleveld 2018, Table 3 (Sigma)" },
  "primary_citation": "eleveld-2018-propofol"
}
```

Each term explicitly labels its scale (`variance` vs `sd`) in the key name, because variance-vs-SD is the most common Σ transcription error (§4).

### 3.3 Model-level: `omega_block` (off-diagonal Ω — correlated random effects)

Optional; present only when the source publishes covariances/correlations between η's:

```jsonc
"omega_block": {
  "correlations": [
    { "between": ["Cl", "V1"], "correlation": 0.60, "citation": "eleveld-2018-propofol" }
  ],
  "complete": false,   // false => only the listed pairs are published; the rest are UNKNOWN, not zero
  "tier": "C",
  "extraction": { "review_status": "unverified" }
}
```

`complete: false` is the honest default: unlisted pairs are treated as **unpublished** (the sampler warns and assumes independence *with a recorded caveat*), never as a confident zero.

### 3.4 Schema-diff summary

| Location | New optional key | Holds |
| --- | --- | --- |
| `parameters[].` | `variability` | per-parameter `omega2` (+ optional `cv_percent`, `shrinkage_percent`, `iov`), its own tier & extraction |
| model root | `residual_error` | Σ structure (`model`, `proportional`/`additive`/`log` terms, tier, extraction) |
| model root | `omega_block` | published off-diagonal Ω correlations + a `complete` flag |
| model root | `variability_status` | rollup: `"none" \| "partial" \| "diagonal" \| "full"` — drives §5 band eligibility and the UI badge |

The JSON-Schema additions live in `dataset/schema/model.schema.json` under new `$defs` (`variability`, `residual_error`, `omega_block`) and four new optional `properties`. `required` lists are unchanged.

---

## 4. The canonical representation and the transcription traps

Random effects are *more* error-prone to transcribe than fixed effects — they are where published-vs-implemented divergence most often hides — so the representation is pinned down to one canonical form and the traps are enumerated for the verifier. This mirrors v0.1 §9's insight that "the covariate equations are the part most worth double-checking," extended to the part that is even harder.

**Trap 1 — variance vs SD vs CV%.** A reported "25%" might be a CV, while the NONMEM `$OMEGA` it came from is a *variance* of 0.0625 — or 0.25, if the author tabulated the SD. Hypnos stores `omega2` (the η-scale variance) as canonical and fixes the log-normal relations:

```
CV ≈ sqrt(exp(omega2) − 1)          (exact, log-normal)
CV ≈ sqrt(omega2)  = ω              (small-ω approximation often used in papers)
```

`hypnos validate` recomputes `cv_percent` from `omega2` (exact form) and flags any record whose stored `cv_percent` disagrees beyond tolerance — catching the variance/SD/CV confusion automatically.

**Trap 2 — scale.** ω² lives on the **log** scale for exponential BSV but on the **natural** scale for additive BSV. The `bsv` block assumes exponential (the field default); any additive-BSV parameter must say so explicitly, because a log-scale variance applied as if natural-scale (or vice-versa) is silently wrong, not loudly broken.

**Trap 3 — Σ variance vs SD.** Identical to Trap 1 for residual error; handled by the scale-in-the-key-name convention (§3.2).

**Trap 4 — shrinkage.** A reported ω² with high η-shrinkage describes the *model's fitted prior*, not well-identified individual variability; we store `shrinkage_percent` when available and surface it as an interpretation caveat (a band drawn from a high-shrinkage ω² is a statement about the model, not a confident statement about people).

**Trap 5 — IOV folded into BSV.** Covered in §2.3: kept separable, never merged.

These five become explicit line items in the verification checklist (§9), each one a yes/no a human confirms against the source table.

---

## 5. Confidence tiers for variability — and the never-synthesize rule

Variability gets its **own** tier, independent of the point-estimate tier, because a model can have a rigorously externally-validated typical-value layer (Tier A) while its random-effects layer is single-cohort, sparsely reported, or only diagonal (Tier C). The general rule, consistent with v0.1's "worst input wins":

> The tier of a **prediction band** is `worst_tier([structural_tier, variability_tier, residual_tier, envelope_floor])`. Because the variability layer is typically the least externally validated component, the band tier is usually **≤** the point (median-line) tier. The median line keeps its v0.1 tier; the *band around it* may be labeled lower. Honest, and visible in the UI.

**Band eligibility, driven by `variability_status`:**

| `variability_status` | Meaning | Band behavior |
| --- | --- | --- |
| `none` | no published random effects (e.g. the classic Marsh TCI parameter set) | **No band drawn.** The model renders as a v0.1 point/line with a "no published variability" annotation. |
| `partial` | some parameters carry BSV; others do not | Band drawn only from the published η's; a warning records which parameters were treated as fixed (band is a *lower bound* on true BSV). |
| `diagonal` | full diagonal Ω, no off-diagonals published | Band drawn; off-diagonals assumed independent **with a recorded caveat** (`omega_block.complete = false`). |
| `full` | diagonal + published off-diagonal block | Band drawn from the full curated Ω. |

**The never-synthesize rule (non-negotiable).** If a model publishes no BSV, Hypnos draws **no band** for it. It does **not** borrow a sibling model's Ω, impute a "typical" CV, or otherwise fabricate uncertainty. A missing band is a true statement ("this model's authors did not publish between-subject variability"); a borrowed band is a lie with error bars. This is the variability-layer analog of v0.1's refusal to simulate a kernel-pending model.

---

## 6. Reference kernels — forward-only, seeded, reproducible

The v0.1 kernels solve the forward problem from typical-value parameters. The variability layer adds **one** new primitive and threads a sample count through the existing API. Nothing here is inverse control (§10).

- **`reference.sample_individual(typical: MicroParams, omega, omega_block, rng) -> MicroParams`** — draw η ~ N(0, Ω) (respecting any curated off-diagonal block; independent-with-caveat otherwise) and return a perturbed `MicroParams`. Each draw is a *virtual individual*; the existing solver simulates it unchanged. This composes cleanly with `MicroParams` (reference.py §43) — it perturbs the structural parameters before the matrix-exponential solve, so the solver, exporters, and validation all keep their single canonical representation.
- **Residual draw** — apply Σ to the predicted trajectory to produce *observation*-level samples when the question is "where will a measured sample land" rather than "where is the individual's true curve."
- **Determinism is mandatory.** Every band-producing call takes an explicit integer `seed`; identical `(seed, dataset_version, request)` → byte-identical quantiles. (The workflow/runtime forbids wall-clock and unseeded RNG precisely so results stay reproducible — the same discipline applies here: no band is ever drawn from an unseeded generator.)

Quantiles (default 5/50/95, configurable) are summarized from `n_samples` draws rather than stored per-draw, so memory stays O(grid), not O(grid × samples).

---

## 7. Uncertainty-aware divergence — the headline feature

This is where v0.2 earns its keep. v0.1's `compare()` (simulate.py §335) overlays every eligible model and `_divergence()` (simulate.py §308) reports the spread of the **point** curves plus the **driver pair**. v0.2 extends the same machinery so the view reports not just *how far apart the lines are* but *whether the bands around them are actually distinguishable.*

### 7.1 The separation index — are the models distinguishable?

For two models A and B at time *t*, let `μ(t)` be each median and `[q_lo(t), q_hi(t)]` each curated-percentile band. Define a band-aware separation:

```
separation(t) = ( q_lo_high(t) − q_hi_low(t) ) / pooled_scale(t)
```

where `high`/`low` are the models ranked by median at *t*. Interpretation:

- **separation > 0** — the bands are **disjoint**: a genuine, irreducible structural disagreement that neither model's own stated variability explains. *Model selection matters here.*
- **separation ≤ 0** — the bands **overlap**: the point-estimate disagreement is within the uncertainty both models already declare. *Model selection matters less than the line-spread suggests.*

`compare()`'s `divergence["cp"]` / `divergence["ce"]` dicts gain:

```jsonc
"separation": {
  "value": -0.31,                 // at the instant of peak median spread (t*)
  "bands_disjoint_at_tstar": false,
  "fraction_trajectory_disjoint": 0.07,   // share of the grid where the chosen-percentile bands separate
  "percentile": [5, 95],
  "band_tier": "C"                // worst variability tier among the two driver models (§5)
}
```

So the view answers, in one machine-readable object: the **driver pair** (which models — v0.1), the **magnitude** (how far apart the lines are — v0.1), and now **whether that gap survives the models' own uncertainty** (v0.2). The dashboard renders bands as shaded ribbons; where they separate, the gap is highlighted — *that* shaded sliver is model-selection risk you cannot variability-away.

### 7.2 Variance decomposition — what dominates the uncertainty here?

A second, complementary readout. Consider the generative process: pick an eligible model, draw a virtual individual (η), add residual (ε). The total predictive variance of the concentration at a fixed (patient, dose, time) decomposes:

```
Var_total  =  Var_structural   +   E[Var_BSV]   +   E[Var_residual]
              (between-model)      (within-model, η)   (Σ)
```

- **Var_structural** — variance of the per-model medians across the eligible set (essentially v0.1's divergence, recast as a variance component).
- **Var_BSV** — average within-model spread from Ω.
- **Var_residual** — from Σ.

Reported as three fractions summing to 1, time-resolved. The punchline: **for some patients model-choice dominates the uncertainty; for others between-subject variability does.** That tells a researcher *when curating more models helps* (structural-dominated regimes) versus *when it doesn't* (BSV-dominated regimes — the honest answer there is "no model will pin this down; the patient is the uncertainty"). This is a genuinely novel, publishable decomposition that the curated dataset makes possible for free.

> **Honesty note.** The decomposition is only as complete as the curated layers. A model with `variability_status: none` cannot contribute a BSV component; a comparison containing such a model reports the decomposition over the band-eligible subset and **names the excluded models**, exactly as v0.1 names kernel-pending and envelope-greyed models rather than silently dropping them.

---

## 8. Export — projecting Ω/Σ to each target's native convention

The variability layer turns several v0.1 exports from "typical value only" into "full population model," and is explicit where a target format cannot represent random effects.

| Format | v0.1 behavior | v0.2 behavior | Status |
| --- | --- | --- | --- |
| **NONMEM** | `$OMEGA 0 FIX` / `$SIGMA 0 FIX` ("variability out of scope") | Emit the curated **`$OMEGA`** diagonal (and **`$OMEGA BLOCK`** when a `complete` `omega_block` spans a contiguous, front-anchored η run — covariance `r·√(ωᵢ²ωⱼ²)`) and **`$SIGMA`** from `residual_error`; non-block or incomplete correlations fall back to a diagonal `$OMEGA` plus a named caveat; no-BSV models keep `... 0 FIX` and name the missing component. | ✅ |
| **PharmML + SO** | typical-value `IndividualParameter`s | Random effects are **first-class**: a `VariabilityModel` maps each η to a `RandomEffect`/`VariabilityLevel` (log transform, variance + CV%), correlations to `Correlation`, and ε to a `ResidualError`. The durable interop anchor now carries the whole NLME object. | ✅ |
| **nlmixr2 / rxode2** | structural model | A runnable `<id>_pop` companion re-expresses the structural params in V/Cl form with log-normal η (`P·exp(eta.X)`), derives the micro-constants, adds the `cp ~ lnorm/prop/add` residual endpoint, and emits the Ω as a `lotri` matrix — `rxSolve(<id>_pop, ev, omega=<id>_omega, nSub=…)`. | ✅ |
| **Pumas** | structural model | A `<id>_pop` `@model` adds `@param` ω² + `@random` η blocks, the V/Cl `@pre` with `exp(η)`, and the residual `dv ~` distribution in `@derived` — Pumas expresses NLME natively. | ✅ |
| **SBML (L3v2)** | deterministic ODE | SBML core cannot express population random effects. Emit the typical-value model as today **plus** the Ω/Σ as the `distrib` package's uncertainty annotations where supported, and as MIRIAM/`hypnos:` RDF otherwise — with an explicit comment that a downstream deterministic SBML consumer sees the typical patient only. Honest about the format's limit, same stance as v0.1's SBML notes. | deferred |
| **TCI-sim JSON** | params + envelope + clinicalUse flag | Add the curated `variability` / `residual_error` blocks verbatim (it is JSON; lossless). | ✅ |

Every variability-carrying export inherits the universal `clinicalUse = "PROHIBITED"` annotation and, additionally, a `hypnos:bandTier` and a `hypnos:variabilityStatus` so a consumer can pin reproducibility *and* know how much of the NLME object is curated.

---

## 9. Validation & verification

**Round-trip (automated, CI).** The exported `$OMEGA`/`$SIGMA` must equal the curated `omega2`/Σ (algebraic, exact). Where the toolchain permits a `$SIMULATION` run, the moments of its simulated population must match the reference Monte-Carlo band within a distributional tolerance (e.g. relative error on the 5th/50th/95th percentiles), the variability-layer analog of v0.1's 1e-4 ODE round-trip.

**Literature validation.** Where a source publishes a **visual predictive check (VPC)** or prediction intervals, the kernel's band is compared to the published interval and the agreement recorded — the variability analog of v0.1's "check the kernel against published example simulations," and what lets the *band* tier (§5) be partly numeric rather than purely editorial.

**Consistency checks (`hypnos validate`).** (1) `cv_percent` recomputes from `omega2` within tolerance (Trap 1). (2) Every `variability`/`residual_error`/`omega_block` `primary_citation` resolves to a real citation record (the same guardrail v0.1 added for `predictive_performance`). (3) `variability_status` matches the actual contents (no `full` without an `omega_block`).

**Human verification — the BSV checklist.** `review_status` for a variability block moves `unverified → verified` only when a human confirms, against the source table, the five traps of §4 as explicit yes/no items:

1. Is the reported number a **variance**, an **SD**, or a **CV%**? On what **scale** (log / natural)?
2. Does the stored `omega2` match the source's `$OMEGA` diagonal (not its square root, not its CV)?
3. Is Σ a **variance** or an **SD**, proportional / additive / combined?
4. Is any reported value shrinkage-inflated (record `shrinkage_percent`)?
5. Is IOV reported separately, and kept separate from BSV?

As in v0.1, **LLMs assist but never promote**: a transcription validated against one VPC figure is not a human PDF verification.

---

## 10. Safety & scope guardrails (non-negotiable, inherited and extended)

The variability layer is forward-only and changes none of v0.1 §10. It adds two guardrails specific to prediction bands, because a percentile band is exactly the kind of artifact that *invites* a dosing question:

- **No quantile-targeting. No inverse control through the back door.** Hypnos will not compute "the dose that keeps the 95th-percentile concentration below *X*," or "the dose for which *P*(BIS < 60) > 0.9," or any other percentile/probability **target**. That is inverse control wearing a statistical hat — precisely the regulated-device function v0.1 forbids. Bands describe a *given* forward dose history; they never invert one.
- **A band is a statement about the model's stated uncertainty, not a claim about a real individual.** The UI and every export label it as such. A 5–95% ribbon is "the published population model says individuals scatter this much," **not** "your patient will be in this range." This distinction is load-bearing: it is the line between an honest research artifact and an implied per-patient guarantee. The dashboard frames the whole feature as *"how distinguishable are these models, accounting for the variability each one declares?"* — a research question — never *"what range will this patient be in?"*

The `clinicalUse = "PROHIBITED"` annotation remains universal. Adding uncertainty bands makes the output *look* more like a clinical tool; the guardrails tighten in proportion, by design.

---

## 11. Worked example (illustrative — numbers pending verification)

The intended end-to-end behavior, on v0.1's own elderly-patient case (72 y, 60 kg, F; the README's left panel, where Marsh/Schnider/Eleveld point estimates already span 169% at the effect-site peak):

```text
$ hypnos compare --drug propofol --age 72 --weight 60 --height 162 --sex F \
                 --bands --percentile 5,95 --samples 2000 --seed 7
included (3, with bands):
  - propofol.eleveld_2018   tier A  band-tier B  Ce peak 3.7 [2.4, 5.6] ug/mL
  - propofol.marsh_1991     tier B  band-tier — (no published BSV; line only)
  - propofol.schnider_1998  tier B  band-tier C  Ce peak 8.0 [5.1, 12.0] ug/mL

divergence (effect-site):
  point spread        5.6 ug/mL (169%)        driver: schnider_1998 vs marsh_1991
  separation@t*       +0.4 (bands DISJOINT)   <- the 169% survives both models' stated BSV
  variance share      structural 0.71 | BSV 0.24 | residual 0.05
  note: marsh_1991 excluded from band metrics (variability_status: none)
```

The reading: Schnider's effect-site band sits clearly above Eleveld's even at the 5–95% level, so this is **real** model-selection risk, not noise — and the variance decomposition confirms *which-model* (0.71) dominates *which-patient* (0.24) here. Marsh is honestly excluded from the band math because its classic parameter set publishes no BSV — shown as a bare line with the reason named, never as a fabricated ribbon.

(Every bracketed number above is illustrative. Per the project's ethos, the actual ω²/Σ values are curated `unverified` and await human PDF confirmation before any of these bands is presented as authoritative.)

---

## 12. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **V0 — Schema + Eleveld** | ✅ done | schema extension (§3); curate the Ω diagonal + Σ for **Eleveld propofol** (the model with the most completely published random-effects structure and an existing Tier-A kernel); `validate` consistency checks. | One model carries a verified-pathway variability block end to end; `hypnos validate` enforces the traps. |
| **V1 — Bands** | ✅ done | `sample_individual` kernel; seeded Monte-Carlo bands in `simulate`; the never-synthesize rule; band tier propagation. | `simulate(..., bands=True, seed=...)` returns reproducible quantiles; no-BSV models correctly draw no band. |
| **V2 — Uncertainty-aware divergence** | ✅ done | separation index + variance decomposition in `compare()`; the headline figure (`docs/images/variability.png`). Dashboard ribbons pending. | The divergence view answers "are the models distinguishable?" and "what dominates the uncertainty?" |
| **V3 — Exports + breadth** | 🟡 export code done; data breadth pending | **done:** the random-effects layer now projects to **every** population format — NONMEM `$OMEGA` diagonal **and** `$OMEGA BLOCK` (off-diagonal, emitted when a `complete` block spans a contiguous, front-anchored η run; honest diagonal-plus-caveat otherwise), PharmML first-class `VariabilityModel` (η → `RandomEffect`/`VariabilityLevel`, ε → `ResidualError`), nlmixr2/rxode2 `_pop` companion (`lotri` Ω + `cp ~ lnorm/prop/add` Σ), Pumas `@random`/`@param` block, and TCI-JSON passthrough. The four shared projections live in `export/_variability.py`. **remaining (data, not code):** backfill BSV for Schnider, the remifentanil trio, dexmedetomidine (the published Schnider CV%s in common toolchains look like RSEs, not BSV — backfill awaits source-table confirmation, per the never-invent rule; the exporters above will carry it the moment it is curated). | The population model — not just its typical patient — round-trips through the pharmacometric formats. |

V0 alone is a useful, self-contained increment: it makes Hypnos the first open resource that curates the random-effects layer with the *same* tiering, citation, and verification discipline it already applies to fixed effects.

---

## 13. Cheat sheet (target API)

```python
import numpy as np, hypnos
ds = hypnos.load()

m = ds["hypnotics_iv.propofol.eleveld_2018"]
m.variability_status                      # "diagonal"
m.param("Cl").variability.omega2          # canonical η-scale variance
m.param("Cl").variability.cv_percent      # derived, checked-consistent
m.residual_error.model                    # "combined"
m.band_tier                               # worst of structural/variability/residual tiers

# Forward simulation WITH a prediction band (seeded, reproducible; still no inverse control)
patient = dict(age=72, weight=60, height=162, sex="F")
schedule = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
traj = hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018",
                       patient=patient, schedule=schedule, t=np.linspace(0, 60, 600),
                       bands=True, percentile=(5, 95), samples=2000, seed=7)
traj.ce_quantiles      # {5: array, 50: array, 95: array}; None if variability_status == "none"
traj.band_tier, traj.warnings

# Uncertainty-aware divergence — the v0.2 headline
cmp = hypnos.compare(ds, drug="propofol", purpose="effect_site",
                     patient=patient, schedule=schedule, bands=True, seed=7)
cmp.divergence["ce"]["separation"]        # {value, bands_disjoint_at_tstar, fraction_trajectory_disjoint, ...}
cmp.divergence["ce"]["variance_share"]    # {structural, bsv, residual}
cmp.excluded_from_bands                   # models with variability_status == "none", named with reasons
```

```bash
hypnos compare --drug propofol --age 72 --weight 60 --height 162 --sex F \
               --bands --percentile 5,95 --samples 2000 --seed 7
hypnos export  --format nonmem --output exports/nonmem/   # now emits real $OMEGA/$SIGMA where curated
hypnos validate                                            # adds the omega2<->cv consistency + scale checks
```

---

## 14. Open questions & explicit deferrals

- **Parameter-uncertainty (SE on θ/Ω) vs. between-subject variability.** This spec curates **BSV** (how individuals scatter). The *estimation* uncertainty of the parameters themselves (standard errors, the covariance step) is a different, also-valuable layer — deferred to a possible v0.3, kept conceptually distinct so the two are never conflated.
- **Correlated residual error / autocorrelation.** Most sources report independent ε; serially-correlated residual models are out of v0.2 scope.
- **Covariate uncertainty.** Bands here propagate Ω/Σ at a *fixed* covariate vector; uncertainty in the covariates themselves (e.g. an estimated weight) is deferred.
- **Bayesian model averaging.** Tempting — combine the eligible models weighted by their performance — but it edges toward producing a single "best" curve, which softens the very model-selection-risk signal v0.1 exists to make visible, and flirts with the dosing-tool line. Deferred, with that reservation recorded.
- **PD variability.** This spec is written around PK BSV; the same structure applies to PD parameters (Ce50, γ), and the schema additions are purpose-agnostic. Backfilling PD random effects (e.g. the Eleveld BIS model) is a natural V3+ extension.
