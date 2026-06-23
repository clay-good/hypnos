# Hypnos — design spec (v0.3): the parameter-estimation-uncertainty layer

**Curate the *estimation* uncertainty of each model's parameters — the standard errors / relative standard errors on the typical-value θ (and, where published, on Ω and Σ themselves), the estimate correlation matrix, and how it was obtained (asymptotic covariance step, bootstrap, SIR) — alongside the typical values (v0.1) and the between-subject random effects (v0.2); propagate it through forward simulation as a seeded *confidence* band (distinct from v0.2's *prediction* band); and split the total predictive uncertainty into its *reducible* part (estimation — shrinks with more data) and its *irreducible* part (between-subject variability — does not).**

> **Status:** **E0–E3 implemented — vocabulary + traps + the reducible/irreducible split (E0), seeded confidence bands (E1), the four-way variance decomposition (E2), and the estimation object in the NONMEM export (E3); Eleveld 2018 propofol now carries SOURCED per-θ estimation uncertainty (pending_human_review).** The schema carries the per-θ `estimation_uncertainty` block (beside `variability`, never inside), the model-root `estimate_covariance`, `uncertainty_status`, and `omega2_se`; `models.py` parses them and exposes `uncertainty_status` / `estimation_band_tier`; `hypnos validate` enforces the numeric traps (§4: scale, RSE↔SE, CI↔SE, status↔contents); `hypnos verify` gains an estimation checklist group; and `compare(..., bands=True)` reports the reducible/irreducible decomposition (§7) over the curated layers. **Not yet curated:** the per-model RSE/covariance values themselves — that is the human-PDF transcription E0 is designed around (the cardinal RSE-vs-CV disambiguation is a human line item, not an LLM call), and it gates the confidence bands of E1. Additive over [v0.2](../v0.2/variability.md) exactly as v0.2 was additive over [v0.1](../v0.1/spec.md): every existing record stays valid, the new blocks are optional, and "not curated" renders as an explicit gap, never as a fabricated interval. Dataset/schema target `0.3.0`. This spec deliberately closes the first open question v0.2 §14 left open: *"Parameter-uncertainty (SE on θ/Ω) vs. between-subject variability … a different, also-valuable layer — deferred to a possible v0.3, kept conceptually distinct so the two are never conflated."* Keeping them distinct is the whole point of this document.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input.* v0.1 made **model-selection** uncertainty first-class. v0.2 made **between-individual** uncertainty first-class. v0.3 makes **estimation** uncertainty first-class — and, crucially, makes the line between *uncertainty you could buy down with more data* and *uncertainty that is the population itself* a measurable, machine-readable thing.

---

## 1. The problem this layer solves — the third hidden uncertainty

A v0.2 forward simulation answers two questions at once: *which model did you pick* (the divergence view), and *where do individuals scatter around that model's typical patient* (the prediction band, drawn from Ω/Σ). A single band still hides a **third** uncertainty, conceptually distinct from the other two:

1. **Structural / model-selection** — *which model?* (v0.1).
2. **Parametric / between-subject (BSV)** — *which patient, within a model?* This is **Ω**: real, physiological, between-person scatter. It does **not** shrink as the derivation cohort grows; a bigger study measures it more precisely but does not make people more alike (v0.2).
3. **Estimation / parameter uncertainty** — *how well is the model itself pinned down?* Every θ in a fitted model is an **estimate** with a standard error. The typical-value clearance is not 1.79 L/min; it is 1.79 ± its SE. This uncertainty **does** shrink as the cohort grows — it is *reducible*. A small single-cohort model (Kataria, n≈40) and a large pooled model (Eleveld, n in the thousands) can publish the *same* typical CL with wildly different confidence in it.

Conflating (2) and (3) is the single most common uncertainty error in applied pharmacometrics, and Hypnos has already tripped over it in print: v0.2 §12 records that *"the published Schnider CV%s in common toolchains look like RSEs, not BSV — backfill awaits source-table confirmation."* That sentence **is** this spec's motivation. A "25%" attached to a parameter might be:

- a **CV of the between-subject distribution** (BSV — a v0.2 ω, irreducible), or
- a **relative standard error** of the estimate (RSE — a v0.3 confidence, reducible).

They are different numbers with different physical meaning, they propagate to different bands, and they are routinely tabulated in adjacent columns of the same table with nothing but a header to tell them apart. v0.2 refused to curate the ambiguous Schnider number until a human resolves the ambiguity. v0.3 gives the dataset the *vocabulary* to record the resolution — a slot for BSV **and** a separate slot for estimation uncertainty — so the confusion can never recur silently.

> The clinical stakes are concrete. Two models with identical typical curves and identical BSV can still differ in trustworthiness: if model A's CL has a 6% RSE (large pooled fit) and model B's has a 60% RSE (tiny cohort, near-unidentifiable), they are not equally good even though their *point* and *prediction-band* views look identical. v0.1 and v0.2 cannot distinguish them. v0.3 can.

The headline deliverable (§7) is the machine-readable answer to: **how much of this band would more data remove, and how much is the patient?** No open anesthesia PK/PD resource answers this today.

---

## 2. What we curate — the estimation-uncertainty structure

A fitted NLME model reports, for its estimated quantities, an *estimation* uncertainty — typically the inverse of the Fisher information (the "$COV step" in NONMEM), or an empirical bootstrap / SIR distribution. We curate it in a canonical, unambiguous form and let exporters project.

### 2.1 Per-parameter estimation uncertainty (SE / RSE / CI on θ)

For each structural θ, where published:

```
θ̂  with standard error SE(θ̂);   RSE% = 100 · SE(θ̂) / |θ̂|;   95% CI ≈ θ̂ ± 1.96·SE (asymptotic)
```

RSE% is the field's most-tabulated form. We store **SE on the estimated scale** as canonical (mirroring v0.2's choice to store the η-scale `omega2`), record the **scale** the SE lives on (natural vs log — many modern fits report θ on a log scale), and derive RSE%/CI with the relation fixed in §4 so they can never silently disagree.

### 2.2 Uncertainty *on the random effects* (SE on Ω and Σ)

The variance components Ω and Σ are themselves estimates with their own SEs. A reported ω² = 0.0625 with a 40% RSE on that variance is a *poorly identified* BSV — its v0.2 band is real but the band's *width* is itself uncertain. We curate `omega2_se` / `sigma_se` where published, as a second-order uncertainty on the v0.2 layer. This is rarely reported in full; the field is optional and usually absent (an honest gap).

### 2.3 The estimate correlation/covariance matrix (the $COV output)

Estimates are correlated — CL and V1 are typically negatively correlated in the fit. Drawing each θ independently from its marginal SE (ignoring these correlations) produces wrong confidence bands. Where a source publishes the **correlation matrix of the estimates** (or the full covariance), we curate it — the direct analog of v0.2's `omega_block`, but for *estimation* covariance rather than *random-effect* covariance. The two are different matrices and are never merged.

### 2.4 The method — because it changes the interpretation

Asymptotic SEs (Fisher/sandwich) assume the log-likelihood is locally quadratic; for skewed or near-boundary parameters (variances especially) they understate uncertainty. **Bootstrap** and **SIR** (sampling-importance-resampling) give empirical, often-asymmetric intervals. We record `method` because a 20% bootstrap RSE and a 20% asymptotic RSE are not the same claim.

---

## 3. Schema extension

Additive and backward-compatible, same discipline as v0.2 §3: every addition is an *optional* property; `additionalProperties: false` is preserved by enumerating new keys; no v0.1/v0.2 record changes meaning. The JSON-Schema additions live in `dataset/schema/model.schema.json` under new `$defs` (`estimation_uncertainty`, `estimate_covariance`) joining the existing `parameter_variability` / `residual_error` / `omega_block` / `range`.

### 3.1 Parameter-level: `estimation_uncertainty`

Added to each entry of `parameters[]`, *beside* (never inside) the v0.2 `variability` block — the physical separation in the schema mirrors the conceptual separation:

```jsonc
{
  "symbol": "Cl",
  "value": { "central": 1.79, "units": "L/min" },
  "variability": { "bsv": { "omega2": 0.0625, "cv_percent": 25.0 }, "tier": "B" },  // v0.2: BSV (irreducible)

  "estimation_uncertainty": {                  // v0.3: estimation (reducible) — a DIFFERENT number
    "se": 0.107,                               // CANONICAL: SE on the estimated scale
    "scale": "natural",                        // "natural" | "log"  (Trap 2)
    "rse_percent": 6.0,                         // OPTIONAL convenience; relation fixed in §4
    "ci95": { "low": 1.58, "high": 2.00 },     // OPTIONAL; asymptotic unless method says otherwise
    "method": "asymptotic_covariance",         // asymptotic_covariance | bootstrap | sir | profile_likelihood | not_reported
    "tier": "B",
    "extraction": {
      "review_status": "unverified",
      "tier_rationale": "RSE column confirmed as estimation RSE (header 'RSE%'), NOT the BSV CV in the adjacent column.",
      "source_locator": "Eleveld 2018, Table 2 (Estimate, RSE%)"
    },
    "primary_citation": "eleveld-2018-propofol"
  }
}
```

### 3.2 Model-level: `estimate_covariance` (the $COV output)

Optional; present only when the source publishes the estimate correlation/covariance:

```jsonc
"estimate_covariance": {
  "correlations": [
    { "between": ["Cl", "V1"], "correlation": -0.42, "citation": "eleveld-2018-propofol" }
  ],
  "complete": false,        // false => only listed pairs published; the rest are UNKNOWN, not zero
  "method": "asymptotic_covariance",
  "covariance_step_succeeded": true,   // NONMEM $COV often fails; record it honestly
  "tier": "C",
  "extraction": { "review_status": "unverified" }
}
```

`complete: false` is the honest default, identical in spirit to `omega_block.complete`: unlisted estimate pairs are treated as **unpublished** (the sampler warns and assumes independence *with a recorded caveat*), never as a confident zero.

### 3.3 Schema-diff summary

| Location | New optional key | Holds |
| --- | --- | --- |
| `parameters[].` | `estimation_uncertainty` | per-θ `se` (+ `scale`, optional `rse_percent`, `ci95`), `method`, its own tier & extraction |
| `parameters[].variability.bsv.` | `omega2_se` | SE on a curated BSV variance (second-order; §2.2) |
| `residual_error.` | `*_se` per term | SE on a curated Σ component (§2.2) |
| model root | `estimate_covariance` | published estimate correlations + `complete` + `covariance_step_succeeded` |
| model root | `uncertainty_status` | rollup: `"none" \| "marginal" \| "correlated"` — drives §5 band eligibility & the UI badge |

`required` lists are unchanged.

---

## 4. The canonical representation and the transcription traps

Estimation uncertainty is *even more* error-prone to transcribe than BSV, because it sits in the same tables, often the same rows, as the very BSV it must not be confused with. The representation is pinned to one canonical form and the traps enumerated for the verifier, extending v0.2 §4.

**Trap 1 — RSE vs BSV-CV (the cardinal sin).** A tabulated "25%" is *either* an estimation RSE *or* a between-subject CV. They are different quantities. The verifier must read the **column header** (and footnotes) and place the number in `estimation_uncertainty` or in `variability`, never guess. This is the trap that froze the Schnider backfill (v0.2 §12); `hypnos validate` cannot auto-resolve it (both are plausible magnitudes), so it is a hard human-verification line item (§9).

**Trap 2 — scale.** An SE on a θ fitted on the **log** scale is not an SE on the natural-scale value. RSE relations differ:

```
natural scale:  RSE% = 100 · SE / |θ̂|
log scale:      RSE% ≈ 100 · SE_log         (SE of log θ̂ ≈ relative SE, for small SE)
```

The `scale` field is mandatory whenever an SE is present; a log-scale SE applied as natural-scale is silently wrong, not loudly broken.

**Trap 3 — SE vs CI vs the multiplier.** Some sources publish 95% CIs, not SEs; the back-conversion `SE = (high − low) / (2·1.96)` holds only for **symmetric, asymptotic** intervals. For bootstrap/SIR CIs (often asymmetric) we store the CI verbatim and leave `se` null rather than fabricate a symmetric SE. `hypnos validate` checks SE↔CI↔RSE consistency only when `method = asymptotic_covariance`.

**Trap 4 — covariance-step failure mistaken for "no uncertainty."** A NONMEM run whose $COV step failed reports *no* SEs — which is *missing*, not *zero*. `covariance_step_succeeded: false` records this so an absent SE is never read as a confident point. (This is the estimation-layer analog of v0.2 Trap 5's "never fold IOV into BSV": never read "not reported" as "negligible.")

**Trap 5 — RSE on a variance vs RSE on an SD.** Mirrors v0.2 Trap 1/3 but for the second-order `omega2_se`/`sigma_se`: an RSE on ω² is not an RSE on ω. The scale-in-the-key-name convention from v0.2 §3.2 is reused.

These become explicit yes/no items in the verification checklist (§9), each confirmed against the source table.

---

## 5. Confidence tiers for estimation uncertainty — and the never-synthesize rule

Estimation uncertainty gets its **own** tier, independent of the point-estimate tier and the variability tier, ranked by *method strength*:

| Tier | Estimation-uncertainty evidence |
| --- | --- |
| **A** | Bootstrap or SIR intervals across the full parameter set; covariance step succeeded; correlations published. |
| **B** | Asymptotic SEs (successful $COV) for the structural parameters, marginal only. |
| **C** | Partial / RSE-only with no covariance; or a single small cohort where SEs are wide and near-boundary. |
| **D** | SE inferred indirectly, or applied to a parameter outside the fit's identifiability — not a defensible confidence statement. |

**Band-tier propagation extends v0.2 §5.** A v0.2 *prediction* band already inherits `worst_tier([structural, variability, residual, envelope_floor])`. A v0.3 *confidence* band (§6) additionally inherits the estimation-uncertainty tier. A combined band (confidence + prediction, the full picture) inherits the worst of all of them. The median line keeps its v0.1 tier; each band kind around it may be labeled lower, and the UI shows which.

**Band eligibility, driven by `uncertainty_status`:**

| `uncertainty_status` | Meaning | Confidence-band behavior |
| --- | --- | --- |
| `none` | no published estimation uncertainty | **No confidence band.** The v0.2 prediction band (if any) is unchanged; the model renders with a "no published estimation uncertainty" annotation. |
| `marginal` | per-θ SEs/RSEs, no correlation matrix | Confidence band drawn from marginal SEs; correlations assumed independent **with a recorded caveat** (`estimate_covariance.complete = false` or absent). |
| `correlated` | SEs + a published estimate covariance | Confidence band drawn from the full curated estimate covariance. |

**The never-synthesize rule (non-negotiable), restated for this layer.** If a model publishes no estimation uncertainty, Hypnos draws **no confidence band**. It does not impute a "typical RSE," borrow a sibling model's covariance, or assume a symmetric SE from a bootstrap CI. A missing confidence band is a true statement ("the authors did not report estimation uncertainty for these parameters"); a fabricated one is a lie with error bars. This is the same rule v0.2 applied to BSV, applied now to estimation uncertainty.

---

## 6. Reference kernels — nested, seeded, reproducible

The v0.2 kernels added `reference.sample_individual(typical, omegas, rng)` (reference.py §221) to draw a virtual individual from Ω. v0.3 adds **one** new primitive and a **nested** sampling loop. Nothing here is inverse control (§10).

- **`reference.sample_parameter_vector(theta, theta_cov, rng) -> ThetaVector`** — draw a perturbed typical-value vector θ\* from the estimate covariance (respecting any curated off-diagonal correlations; independent-with-caveat otherwise), on the curated `scale` (log-normal for log-scale θ, normal for natural-scale). Each draw is a *plausible alternative fit* of the same model. This composes with `MicroParams.from_volumes_clearances` (reference.py §62) exactly as `sample_individual` does — it perturbs the typical values *before* the structural→micro conversion, so the solver, exporters, and validation keep their single canonical representation.

- **Two band kinds, never conflated:**
  - A **confidence band** answers *"where is the typical curve, given how well the model is estimated?"* — outer loop draws θ\* from the estimate covariance; the *typical* trajectory is simulated for each θ\* (no η, no ε). Quantiles over θ\* draws.
  - A **prediction band** (v0.2) answers *"where will an individual / a measured sample land?"* — η draws (and optionally ε).
  - A **combined band** nests them: outer θ\* (estimation) × inner η (BSV) × optional ε (residual). This is the honest full predictive interval; its three layers are recorded so the §7 decomposition can read them apart.

- **Determinism is mandatory**, inherited verbatim from v0.2 §6: every band call takes an explicit integer `seed`; identical `(seed, dataset_version, request)` → byte-identical quantiles. The nested loop seeds a child generator per outer draw (`rng.spawn`-style) so the inner BSV stream is reproducible *and* independent across outer draws. No band is ever drawn from an unseeded generator. Memory stays O(grid × quantiles), not O(grid × outer × inner), by streaming quantile accumulation.

---

## 7. The reducible/irreducible decomposition — the headline feature

This is where v0.3 earns its keep. v0.2 §7.2 decomposed the total predictive variance into three time-resolved shares: **structural** (between-model), **BSV** (within-model η), **residual** (Σ). v0.3 splits the picture along the axis that matters for *what to do about it*:

```
Var_total = Var_structural + Var_estimation + Var_BSV + Var_residual
            \_____________ reducible ______________/   \___ irreducible ___/
              (more models)    (more data)              (the population is the limit)
```

- **Var_estimation** — the new component: average within-model spread of the *typical* curve across θ\* draws. It is **reducible**: a larger or better-designed study shrinks it.
- **Var_structural** and **Var_estimation** are both reducible in principle (curate/fit more), but by *different* actions — more *models* vs. more *data per model*.
- **Var_BSV** and **Var_residual** are **irreducible**: no amount of data makes people identical or assays noiseless.

`compare()`'s `divergence["cp"]`/`["ce"]` (built in `_band_divergence`, simulate.py §527, whose `variance_share` key v0.2 already populates) gains a fourth share and a reducibility rollup:

```jsonc
"variance_share": { "structural": 0.61, "estimation": 0.18, "bsv": 0.17, "residual": 0.04 },
"reducibility":   { "reducible": 0.79, "irreducible": 0.21 },   // (structural+estimation) vs (bsv+residual)
"estimation_band_tier": "C"
```

The punchline a single curve can never give: **for this patient at this instant, 79% of the predictive uncertainty could be bought down with more research (more models, more data) and 21% is the patient and the assay — irreducible.** That tells a researcher *whether the answer is "fund a bigger study" or "no model will ever pin this person down."* It is a genuinely novel, publishable decomposition that the curated dataset — now carrying all three uncertainty layers with consistent tiering — makes possible for free.

> **Honesty note (inherited from v0.2 §7).** The decomposition is only as complete as the curated layers. A model with `uncertainty_status: none` contributes no estimation component; the readout is computed over the eligible subset and **names the excluded models**, exactly as v0.1 names envelope-greyed models and v0.2 names no-BSV models, rather than silently dropping them.

---

## 8. Export — projecting estimation uncertainty to each target's native convention

The estimation layer is, more than any other, where the export targets *diverge* on what they can represent. The PharmML **Standard Output (SO)** is the one format designed precisely for this; the rest carry it as annotation or companion.

| Format | v0.3 behavior | Status |
| --- | --- | --- |
| **PharmML + SO** | **First-class, and the durable anchor.** The SO's `Estimation/PopulationEstimates` carries `MLE` + `StandardError` / `RelativeStandardError`, and `CovarianceMatrix` / `CorrelationMatrix` natively. Hypnos emits the curated SEs, RSEs, CIs, method, and any `estimate_covariance` into exactly these elements. The whole estimation object round-trips. | planned |
| **NONMEM** | The estimate covariance is `$COV` *output* (a `.cov`/`.cor`/`.ext` file), not a control-stream input, so it cannot live in the `.mod`. Hypnos emits a companion `<id>.cov`-style block plus a control-stream comment, and the same SEs/correlations as `hypnos:` RDF in the export sidecar. `covariance_step_succeeded: false` is stated explicitly. | planned |
| **nlmixr2 / rxode2** | The structural `_pop` companion (v0.2) gains a curated `theta_cov` matrix usable to seed SIR/parametric resampling of θ in a simulation script; emitted as a `lotri`-style block beside the v0.2 `omega` matrix, clearly labeled *estimation* covariance, not Ω. | planned |
| **Pumas** | A companion estimate-covariance object suitable for `simobs`-style parametric uncertainty propagation; documented as estimation uncertainty, distinct from the `@random` η block. | planned |
| **SBML (L3v2)** | SBML core is deterministic; the SEs/RSEs/estimate correlations attach as MIRIAM/`hypnos:` RDF (`hypnos:parameterStandardError` / `hypnos:estimateCorrelation`), with an explicit comment that a deterministic consumer sees the point estimate only. (`distrib`-package `UncertML` objects are a possible future refinement; the RDF carrier is lossless and tool-portable today.) Same stance as v0.1/v0.2 SBML notes. | planned |
| **TCI-sim JSON** | The `estimation_uncertainty` / `estimate_covariance` blocks pass through verbatim (it is JSON; lossless). | planned |

Every estimation-carrying export inherits the universal `clinicalUse = "PROHIBITED"` annotation and, additionally, a `hypnos:estimationBandTier` and `hypnos:uncertaintyStatus`, so a consumer can pin reproducibility *and* know how much of the estimation object is curated.

---

## 9. Validation & verification

**Round-trip (automated, CI).** Exported SEs/RSEs/correlations must equal the curated values (algebraic, exact). Where a toolchain permits a parametric-uncertainty simulation (SIR / `$SIMULATION` with parameter resampling), the moments of its θ\*-driven typical-curve spread must match the reference confidence band within a distributional tolerance — the estimation-layer analog of v0.2 §9's prediction-band round-trip.

**Internal consistency (`hypnos validate`).** Extends v0.2's `_check_variability` (validate.py §122) with: (1) `rse_percent` recomputes from `se` and `value.central` on the declared `scale`, within tolerance (Trap 2); (2) `ci95` is consistent with `se` when `method = asymptotic_covariance` (Trap 3); (3) every `estimation_uncertainty` / `estimate_covariance` `primary_citation` resolves to a real citation record; (4) `uncertainty_status` matches the contents (no `correlated` without an `estimate_covariance`).

**Human verification — the estimation checklist.** `review_status` moves `unverified → verified` only when a human confirms, against the source table, as explicit yes/no items extending the v0.2 §9 list (the `_checklist_for` builder, verification.py §54, gains an `estimation` group beside `structural`/`covariate`/`population`/`envelope`/`citation`):

1. Is this number an **estimation RSE/SE** or a **between-subject CV**? (Trap 1 — the cardinal disambiguation; read the column header.)
2. On what **scale** is the SE (natural / log)? (Trap 2)
3. Is it an **SE**, a **CI**, or an **RSE**, and by what **method** (asymptotic / bootstrap / SIR)? (Trap 3)
4. Did the **covariance step succeed**? Are absent SEs *missing* or *zero*? (Trap 4)
5. For any `omega2_se`/`sigma_se`: RSE on the **variance** or the **SD**? (Trap 5)

As in v0.1/v0.2, **LLMs assist but never promote**: resolving an RSE-vs-CV ambiguity is precisely the human-judgment line item, not an LLM call.

---

## 10. Safety & scope guardrails (non-negotiable, inherited and extended)

The estimation layer is forward-only and changes none of v0.1 §10 / v0.2 §10. A *confidence* band, like a prediction band, is exactly the kind of artifact that invites a dosing question, so the v0.2 guardrails extend verbatim plus one clarification:

- **No inverse control through a confidence band.** Hypnos will not compute "the dose for which the typical-curve confidence interval stays below *X*," any more than it will quantile-target a prediction band (v0.2 §10). Bands describe a *given* forward dose history; they never invert one.
- **A confidence band is a statement about estimation precision, not a per-patient guarantee, and not the same as a prediction band.** The UI and every export label the two band kinds distinctly: a confidence ribbon says *"this is how well the typical curve is pinned down by the data the model was fit on"*; a prediction ribbon says *"this is how individuals scatter."* Reading one as the other is a category error the labeling exists to prevent.
- The `clinicalUse = "PROHIBITED"` annotation remains universal. As before, the more the output *looks* like a clinical tool, the tighter the guardrails — by design.

---

## 11. Worked example (illustrative — numbers pending verification)

The intended end-to-end behavior on v0.1/v0.2's elderly-patient case (72 y, 60 kg, F):

```text
$ hypnos compare --drug propofol --age 72 --weight 60 --height 162 --sex F \
                 --bands prediction,confidence --percentile 5,95 --samples 2000 --seed 7
included (3):
  - propofol.eleveld_2018   tier A  pred-band-tier B  est-band-tier C
        Ce peak 3.57   prediction [0.87, 6.90]   confidence [3.31, 3.86]
  - propofol.marsh_1991     tier B  (no published BSV; no published estimation uncertainty; line only)
  - propofol.schnider_1998  tier B  (BSV unresolved RSE/CV; no confidence band drawn — Trap 1)

effect-site divergence:
  variance share @ t*=3.33 min:
        structural 0.58 | estimation 0.05 | BSV 0.33 | residual 0.04
  reducibility:  reducible 0.63  |  irreducible 0.37
  reading: most of the spread here is WHICH-MODEL (structural, 0.58) — curating/validating
           models helps; Eleveld's own typical curve is tightly estimated (estimation 0.05),
           so 'more data for Eleveld' would barely move it.
  note: marsh_1991 & schnider_1998 excluded from confidence-band math (uncertainty_status: none / unresolved)
```

The reading: Eleveld's *confidence* band `[3.31, 3.86]` is tight (it is a large pooled fit — its typical curve is well pinned), while its *prediction* band `[0.87, 6.90]` is wide (people genuinely scatter). The decomposition makes the difference between "the model is uncertain" and "the population is wide" a number, not a vibe. Schnider is honestly excluded from the confidence math because its tabulated variability is the very RSE-vs-CV ambiguity (Trap 1) that v0.2 refused to guess at — shown as a reason, never as a fabricated ribbon.

(Every bracketed number above is illustrative. Per the project ethos, the actual SE/covariance values are curated `unverified` and await human PDF confirmation before any of these bands is presented as authoritative.)

---

## 12. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **E0 — Schema + Eleveld** | 🟢 machinery shipped | schema extension (§3) ✅; `validate` consistency checks (Traps 2–4) ✅; `models.py` parsing + `uncertainty_status`/`estimation_band_tier` ✅; `verify` estimation group ✅; the §7 reducible/irreducible split in `compare` ✅. **Pending:** curate per-θ SEs/RSEs + method for **Eleveld propofol** — the human-PDF transcription (the RSE-vs-CV call is a human line item; the Eleveld full text is paywalled to automated fetch). | Schema enforces the new traps and the RSE-vs-CV separation by construction (done); one model carries a verified-pathway estimation block (awaits PDF curation). |
| **E1 — Confidence bands** | 🟢 shipped | `sample_parameter_vector` kernel; seeded confidence bands in `simulate`; the never-synthesize rule; the two-band-kinds distinction; estimation-band-tier propagation. | `simulate(..., bands="confidence", seed=...)` returns reproducible quantiles distinct from the prediction band; no-SE models correctly draw none. |
| **E2 — Reducible/irreducible decomposition** | 🟢 shipped | the four-way `variance_share` + `reducibility` rollup in `compare()`; nested combined band; the headline figure (`docs/images/reducibility.png`); the dashboard readout naming reducible vs irreducible shares and the excluded models. | The divergence view answers "how much of this is buy-down-able with more data vs. the patient?" |
| **E3 — Exports + correlation breadth** | 🟢 NONMEM estimation block shipped; PharmML-SO + cross-model SE backfill remain | PharmML/SO first-class estimation object; NONMEM `.cov` companion + RDF; nlmixr2/Pumas SIR-seed covariance; **backfill** SEs + (where published) estimate correlation matrices for Schnider, the remifentanil trio, dexmedetomidine — *resolving the RSE-vs-CV ambiguity v0.2 left open per drug.* | The estimation object round-trips through the formats; the Schnider RSE/CV question is closed by curation, not guesswork. |

E0 alone is a useful, self-contained increment: it makes Hypnos the first open resource to curate estimation uncertainty *separately and explicitly* from between-subject variability — closing, by construction, the conflation that has dogged the field's parameter tables.

---

## 13. Cheat sheet (target API)

```python
import numpy as np, hypnos
ds = hypnos.load()

m = ds["hypnotics_iv.propofol.eleveld_2018"]
m.uncertainty_status                       # "marginal" | "correlated" | "none"
m.param("Cl").estimation_uncertainty.rse_percent   # estimation RSE (NOT the BSV cv_percent)
m.param("Cl").variability.cv_percent               # between-subject CV — a different number
m.estimation_band_tier                     # worst of structural/estimation tiers

# Forward simulation with a CONFIDENCE band (typical curve, estimation uncertainty only)
patient = dict(age=72, weight=60, height=162, sex="F")
schedule = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
traj = hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018",
                       patient=patient, schedule=schedule, t=np.linspace(0, 60, 600),
                       bands="confidence", percentile=(5, 95), samples=2000, seed=7)
traj.ce_confidence      # {5: array, 50: array, 95: array}; None if uncertainty_status == "none"
traj.ce_quantiles       # the v0.2 PREDICTION band — distinct object, distinct meaning

# The reducible/irreducible decomposition — the v0.3 headline
cmp = hypnos.compare(ds, drug="propofol", purpose="effect_site",
                     patient=patient, schedule=schedule, bands="combined", seed=7)
cmp.divergence["ce"]["variance_share"]    # {structural, estimation, bsv, residual}
cmp.divergence["ce"]["reducibility"]      # {reducible, irreducible}
```

```bash
hypnos compare  --drug propofol --age 72 --weight 60 --height 162 --sex F \
                --bands prediction,confidence --percentile 5,95 --samples 2000 --seed 7
hypnos export   --format pharmml --output exports/pharmml/   # SO now carries SEs/RSEs/estimate covariance
hypnos validate                                              # adds RSE<->SE<->CI checks + RSE-vs-CV separation
```

---

## 14. Open questions & explicit deferrals

- **Full Bayesian posteriors vs. frequentist SEs.** Asymptotic SEs are (loosely) a Laplace approximation to a posterior under a flat prior; some modern fits publish full MCMC posteriors. v0.3 curates the SE/covariance summary; ingesting a published posterior sample as the estimation object is a v0.4-or-later refinement, kept conceptually distinct so a posterior is never silently treated as a frequentist SE.
- **Reconstructing covariance from marginals.** When only marginal RSEs + a correlation matrix are published (common), the full covariance is recoverable; when *only* marginals are published, the off-diagonals are unknown (`complete: false`) — we never assume independence as truth. SIR from such partial information is a documented-caveat approximation, not a curated covariance.
- **Estimation uncertainty on covariate coefficients.** This spec curates SEs on the disposition θ; SEs on covariate-equation coefficients (the age/weight/LBM slopes) are in-scope-in-principle and curated where the source gives them, but their propagation interacts with v0.2 §14's deferred *covariate uncertainty* — the two are tracked together when that layer lands.
- **Relationship to external validation.** Estimation uncertainty (this spec) is *internal* — how well the data the model was fit on pins it down. *External* predictive performance against independent data (MDPE/MDAPE) is a different, complementary check; that empirical layer is the subject of [v0.4](../v0.4/external_validation.md). A model can be tightly estimated internally yet predict poorly externally; Hypnos keeps both, never conflated.
