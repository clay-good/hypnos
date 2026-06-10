# Hypnos — design spec (v0.7): the covariate-model-uncertainty layer

**Curate the *covariate sublayer* of each model as a first-class, named, cited object — the auxiliary body-size/composition equations a model is *derived with* (lean body mass, fat-free mass, body-surface area, ideal/normal-fat-mass, allometric size) — because the choice of, e.g., the James-1976 vs. Janmahasatian-2005 lean-mass equation silently reparameterizes a model and is the single most common source of "same model, different pump, different concentration" divergence in the field; propagate the *uncertainty in the input covariates themselves* (an estimated weight, a rounded height) through forward simulation as a seeded *covariate band*; and add a fourth, separately-tiered component — `covariate` — to the v0.2/v0.3 variance decomposition, splitting it honestly into its reducible (measure-better / pick-the-right-equation) and structural (which-equation) parts.**

> **Status: C0+C1+C2 shipped (dataset/schema `0.7.0`); C3 proposed.** The covariate-equation library, the `covariate_model` bindings, `hypnos.covariates.evaluate`, the validator + verification group, the `covariate_equations.png` figure (C0); the equation-divergence view — `covariate_divergence` + `hypnos covariate-divergence` + the `covariate_divergence.png` figure (C1); and the covariate-value band + the fifth variance component — `sample_covariate_vector`, `simulate(..., bands=["covariate"])`, the `covariate`/`covariate_split` decomposition in `compare`, the CLI/dashboard surfaces, and the `covariate_band.png` figure (C2) — are implemented and CI-green (see the §12 roadmap). C3 (covariate-aware exports) remains design-only. It is written to be additive over [v0.1](../v0.1/spec.md)–[v0.6](../v0.6/local_anesthetics.md) exactly as every prior layer was additive: no existing record changes meaning, the new blocks are optional, and "not curated" renders as an explicit gap, never as a fabricated equation or band. Dataset/schema target `0.7.0`. This spec closes the *covariate-uncertainty* deferral the project has named twice and never resolved: v0.2 §14 ("**Covariate uncertainty.** Bands here propagate Ω/Σ at a *fixed* covariate vector; uncertainty in the covariates themselves … is deferred.") and v0.3 §14 ("their propagation interacts with v0.2 §14's deferred *covariate uncertainty* — the two are tracked together when that layer lands"). This document *is* that layer. It also promotes to first-class status a field v0.1 already sketched but never formalized: the model record's `covariates.lbm_equation: "james_1976"` (v0.1 §4) — v0.7 turns that lone string into a validated library of named covariate kernels with their own validity envelopes, tiers, and a divergence view.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input.* v0.1 made **model-selection** uncertainty first-class; v0.2 **between-subject**; v0.3 **estimation**; v0.4 measured error **against the world**. v0.7 makes first-class the uncertainty that sits *underneath* all of them — the **covariate transform** every covariate-scaled model silently depends on. A population model is not `θ`; it is `θ` *plus the equations that turn a patient into the numbers θ scales*. Curating θ while leaving those equations implicit is the same category of omission v0.2 named for random effects: a load-bearing part of the model, dropped.

---

## 1. The problem this layer solves — the fourth hidden uncertainty

A v0.1–v0.3 forward simulation fixes a covariate vector — age, weight, height, sex — and propagates structural, between-subject, and estimation uncertainty around it. That picture still hides a **fourth** uncertainty, conceptually distinct from the other three, and it hides in two places at once:

1. **Structural / model-selection** — *which model?* (v0.1)
2. **Parametric / between-subject (BSV)** — *which patient, within a model?* (v0.2)
3. **Estimation** — *how well is the model pinned down?* (v0.3)
4. **Covariate-model** — *which equation turns this patient into the model's inputs, and how well do we even know the patient's covariates?* (v0.7)

Most covariate-scaled anesthesia models do not scale on raw weight. They scale on a **derived** body-size or body-composition quantity computed by a *named auxiliary equation*: Schnider's clearances and volumes scale on **lean body mass** computed by the **James (1976)** formula; modern general-purpose models (Eleveld) scale on **fat-free mass**, frequently the **Janmahasatian (2005)** formula or the pediatric-aware **Al-Sallami (2015)** one; some volatile and dosing conventions use **body-surface area** (Du Bois, Mosteller) or **ideal body weight** (Devine). These equations are *part of the model*. Substitute a different one and you have silently fitted a different model — the typical-value curve moves with no change to a single published θ.

This is not hypothetical, and it is not small. It is the best-documented "same model, different answer" divergence in the field:

- **The James-LBM inversion.** James's lean-body-mass equation contains a `−k·(weight/height)²` term. As weight rises at fixed height, computed LBM eventually *peaks and then decreases* — non-physically. Schnider-on-James therefore misbehaves in obesity, which Hypnos already curates as a `known_failure_mode` and auto-tiers to D (v0.1 §4 worked example; the README divergence panel greys Schnider for the obese). v0.7's insight: **that failure mode is a property of the covariate equation, not of the disposition parameters** — and it is invisible unless the equation is a first-class, swappable, validity-bounded object.
- **Pump-to-pump disagreement.** Commercial TCI implementations of nominally "the Schnider model" have differed in their LBM handling (some cap it, some substitute Janmahasatian, some clamp the input weight), producing materially different predicted concentrations *for the identical published parameter set*. The field has written about this directly; Hypnos has the machinery to make it a measured divergence rather than folklore.
- **Covariate values are themselves uncertain.** Weight is often estimated, not measured; height is rounded or self-reported; age may be the only hard number. v0.1–v0.3 treat the covariate vector as exact. A 5–10 kg weight uncertainty propagates straight through a weight-scaled clearance into the prediction, and no current layer accounts for it.

> The two halves are genuinely distinct and v0.7 keeps them separate, the way v0.3 kept estimation distinct from BSV. **Covariate-equation uncertainty** is *which transform* (James vs. Janmahasatian) — a structural choice, reducible by curating/validating the right equation. **Covariate-value uncertainty** is *how well we know the input* (weight ± SD) — reducible by measuring better. Both are reducible (in the v0.3 sense), by different actions, and neither is the patient-scatter (BSV) or the assay noise (Σ) that no amount of care removes.

The headline deliverable (§7) is the machine-readable answer to: **for this patient and model, how much of the prediction is riding on the choice of body-composition equation — and is that choice moving the curve more than the model's own stated variability?** No open resource answers this today; most do not even record which equation a model was derived with.

---

## 2. What we curate — the covariate sublayer

A covariate-scaled model has three curatable covariate-layer facts, each tiered and cited like everything else.

### 2.1 The covariate equations the model was *derived with* (the binding fact)

The non-negotiable curation: each model names, per derived covariate it consumes, the **exact published equation** its authors used, by a stable identifier into a shared equation library (§6). Schnider → `lbm: james_1976`; a model that fits on FFM → `ffm: janmahasatian_2005` or `ffm: al_sallami_2015`; a BSA-scaled convention → `bsa: du_bois_1916`. This is the single most important field this spec adds, because it is the one whose silent substitution is the documented failure mode. It is a *binding* fact: the validator (§9) checks that the model is simulated through the named equation and no other.

### 2.2 Each equation's own validity envelope (equations have failure modes too)

A covariate equation is itself a model with a derivation population and a range outside which it misbehaves. James-LBM was derived in non-obese adults and inverts above BMI ≈ 37–42; Al-Sallami FFM extends validity into childhood where Janmahasatian was not derived; Devine IBW is undefined below a height threshold. We curate each equation's `validity_envelope` and its own `known_failure_modes` (the James inversion is the canonical entry), reusing the existing `Envelope` / `FailureMode` `$defs` (models.py §264, §233). An equation evaluated outside its envelope tiers the covariate layer — and therefore the composed simulation — down to D, exactly as a disposition envelope violation does.

### 2.3 Covariate-value uncertainty (optional, caller-supplied, never invented)

Where a caller supplies a distribution on an input covariate (e.g. `weight = {mean: 70, sd: 6, dist: "normal"}`), Hypnos propagates it. This is a *caller input*, not a dataset field — Hypnos curates no patient's weight. The dataset's only role is to declare, per model, **which covariates each is sensitive to** (its `covariates.required`, v0.1 §4, already records this), so the propagation knows what matters. No covariate-value distribution is ever assumed; absent a caller distribution, the covariate is treated as exact and the value band is simply not drawn (the never-invent rule, §5).

---

## 3. Schema extension

Additive and backward-compatible, same discipline as v0.2 §3 / v0.3 §3: every addition is an *optional* property; `additionalProperties: false` is preserved by enumerating new keys; no v0.1–v0.6 record changes meaning. The JSON-Schema additions live in `dataset/schema/model.schema.json` under new `$defs` (`covariate_equation_ref`, `covariate_model`) joining the existing library, plus a new top-level dataset artifact `dataset/covariate_equations/*.json` (the equation library, §6).

### 3.1 Model-level: `covariate_model` (which equations this model binds to)

Promotes and structures v0.1's lone `covariates.lbm_equation` string:

```jsonc
"covariate_model": {
  "derived_inputs": [
    {
      "quantity": "lbm",                       // lbm | ffm | bsa | ibw | nfm | allometric_size
      "equation": "james_1976",                // resolves to dataset/covariate_equations/james_1976.json
      "used_for": ["V1", "Cl", "Q2"],          // which θ scale on it (drives §7 sensitivity)
      "verbatim": true,                        // true => the model's authors used exactly this; false => a documented substitution
      "tier": "B",
      "extraction": {
        "review_status": "unverified",
        "tier_rationale": "Schnider 1998 Methods: LBM by James; confirm sex-specific coefficients and the (W/H)^2 term.",
        "source_locator": "Schnider 1998, p. 1503 (Methods, LBM)"
      },
      "primary_citation": "schnider-1998-propofol-pk"
    }
  ],
  "covariate_sensitivity_status": "declared"   // "declared" | "computed" | "none" — drives §7 & the UI badge
}
```

### 3.2 Equation-library record: `covariate_equation` (a new dataset artifact)

One JSON per equation in `dataset/covariate_equations/`, each a small structured object — the equations are *shared* across models (many models use James), so they are curated once, exactly as `dataset/citations/` are shared:

```jsonc
{
  "id": "james_1976",
  "quantity": "lbm",
  "form": "1.1*weight - 128*(weight/height_cm)^2",   // verbatim, machine-parseable, sex-branched below
  "sex_specific": {
    "male":   "1.1*weight - 128*(weight/height_cm)^2",
    "female": "1.07*weight - 148*(weight/height_cm)^2"
  },
  "units": { "weight": "kg", "height_cm": "cm", "result": "kg" },
  "validity_envelope": { "bmi_kg_m2": { "min": 18, "max": 37 } },
  "known_failure_modes": [
    {
      "condition": "bmi_kg_m2 > ~37 (sex-dependent)",
      "behavior": "the -(W/H)^2 term dominates; computed LBM peaks then DECREASES with weight (non-physical inversion)",
      "action": "tier_down_to_D and emit warning on simulate/export",
      "citation": "absalom-2009-schnider-obese"
    }
  ],
  "tier": "B",
  "primary_citation": "james-1976-lbm"
}
```

### 3.3 Schema-diff summary

| Location | New optional key | Holds |
| --- | --- | --- |
| model root | `covariate_model` | the `derived_inputs[]` binding each consumed derived covariate to a named equation (+ `used_for`, `verbatim`, tier, extraction) and a `covariate_sensitivity_status` rollup |
| **new** `dataset/covariate_equations/*.json` | (record) | a shared, citation-backed equation library: `form`, sex branches, units, `validity_envelope`, `known_failure_modes`, tier |
| `applicability_envelope.` | (reuse) | no new keys; equation `validity_envelope`s compose with the disposition envelope via worst-tier (§5) |

`required` lists are unchanged. A model with no `covariate_model` block is treated as *covariate-equation-unspecified* — an explicit gap surfaced in the UI, not an assumption that raw weight was used.

---

## 4. The canonical representation and the transcription traps

The covariate layer is *uniquely* error-prone because it is the part most often re-implemented from memory rather than copied from a table — the equations feel like "standard formulas," so they get silently swapped. The representation is pinned to the named-equation-library form and the traps enumerated for the verifier, extending v0.1 §9's "the covariate equations are the part most worth double-checking."

**Trap 1 — the silent equation substitution (the cardinal sin).** Implementing Schnider with Janmahasatian FFM instead of James LBM (or vice-versa) changes the prediction with no change to any published θ, and leaves no trace. `verbatim: true` asserts the curated equation *is* the authors' choice; a `false` value flags a deliberate, cited substitution. This is the covariate-layer analog of v0.3 Trap 1's RSE-vs-CV: a number/equation that looks interchangeable but is not, and the verifier must read the **Methods**, not guess.

**Trap 2 — units in the equation.** Height in cm vs. m flips the `(W/H)²` term by 10⁴; BSA equations mix cm and m by convention. Every equation record carries explicit `units` per variable, and the validator checks dimensional consistency against the model's covariate units.

**Trap 3 — sex coding.** James, Janmahasatian, and Al-Sallami all have sex-specific coefficients; a model that applied the male branch to all subjects (or used a different sex encoding) is silently wrong. The `sex_specific` branches are curated explicitly and the sex-coding convention recorded.

**Trap 4 — evaluating an equation outside its own validity envelope.** James above BMI ≈ 37 is not "a bit off," it is *qualitatively wrong* (inversion). The equation's `validity_envelope` is enforced independently of the disposition envelope; a patient inside the disposition envelope but outside the *equation's* envelope still tiers the covariate layer to D (§5). This makes v0.1's Schnider-obesity failure mode a property of the right object.

**Trap 5 — fixed vs. fitted exponents (allometric overlap with v0.8).** An allometric size term may use the theory-based ¾-power exponent (fixed) or an empirically fitted exponent; treating a fitted exponent as the fixed theoretical one (or vice-versa) misstates the scaling. Recorded explicitly; the developmental-pharmacology spec ([v0.8](../v0.8/developmental.md)) leans on the same field.

These five become explicit yes/no items in the verification checklist (§9), each confirmed against the source Methods section.

---

## 5. Confidence tiers for the covariate layer — and the never-invent rule

The covariate layer gets its **own** tier, independent of the structural/variability/estimation tiers, because a model with a Tier-A disposition fit can rest on a covariate equation (James) that is Tier-B at best and Tier-D outside its envelope. The propagation rule extends v0.1's "worst input wins" and v0.3 §5's band-tier chain:

> The tier of a composed simulation is `worst_tier([structural, variability, estimation, covariate_equation, envelope_floor])`. The covariate-equation tier enters the chain like any other component, and an equation evaluated outside its own `validity_envelope` forces the whole result to D — exactly as a disposition envelope violation does.

**Band eligibility, driven by `covariate_sensitivity_status`:**

| `covariate_sensitivity_status` | Meaning | Behavior |
| --- | --- | --- |
| `none` | the model scales on raw covariates only (no derived equation) | No covariate-equation divergence axis; a caller-supplied covariate-value band is still drawn if requested. |
| `declared` | the model names its derived-input equations (§3.1) | The equation-divergence view (§7.1) is available; swapping among *admissible* equations is shown, each greyed where outside its own envelope. |
| `computed` | declared **and** a covariate-value distribution was supplied | The covariate band (§7.2) is drawn and enters the variance decomposition (§7.3). |

**The never-invent rule (non-negotiable), restated for this layer.** Hypnos does **not** auto-select a covariate equation to flatter a model, impute a covariate-value distribution, or substitute a "better" equation silently. If a model's derivation equation is uncurated, the field is absent (an honest "covariate equation unspecified"); if no covariate-value distribution is supplied, no value band is drawn. A substituted equation is only ever shown as an explicitly-labeled, cited `verbatim: false` alternative in the *divergence* view — never as the model's own prediction. This is the v0.2 never-synthesize rule applied to the covariate transform.

---

## 6. Reference kernels — a shared, validity-bounded equation library

The v0.1 kernels already contain the seed of this: `reference.lbm_james` (reference.py §28) and `reference.bmi` (§36) exist today, used inside the Schnider covariate model. v0.7 generalizes that one function into a small, pure, individually-validated library and threads covariate-value sampling through the existing band machinery. Nothing here is inverse control (§10).

- **`covariates.evaluate(equation_id, patient) -> float`** — dispatch to the named, pure equation (James, Janmahasatian, Al-Sallami, Du Bois, Mosteller, Devine, allometric size), each a verbatim transcription of the `form` field, each checking its own `validity_envelope` and attaching the equation's `known_failure_mode` warning when violated. `lbm_james` (§28) becomes the first registered member, unchanged in behavior.
- **`covariates.sample_covariate_vector(patient, covariate_uncertainty, rng) -> patient*`** — draw a perturbed covariate vector from a caller-supplied distribution (respecting declared correlations, e.g. weight–height; independent-with-caveat otherwise), seeded. Composes with the existing pipeline by perturbing the covariate vector *before* `_covariate_env` (simulate.py §88) and the structural→micro conversion (`MicroParams.from_volumes_clearances`, reference.py §62), so the solver, exporters, and validation keep their single canonical representation — the same composition point `sample_individual` (reference.py §221) and v0.3's `sample_parameter_vector` use.
- **Three band kinds now compose**, never conflated: the v0.2 *prediction* band (η/ε), the v0.3 *confidence* band (θ\*), and the v0.7 *covariate* band (covariate-value draws). A fully-nested band is θ\* × covariate × η × ε; each layer is recorded so the §7.3 decomposition reads them apart.
- **Determinism is mandatory**, inherited verbatim from v0.2 §6 / v0.3 §6: every band call takes an explicit integer `seed`; identical `(seed, dataset_version, request)` → byte-identical quantiles; nested generators are spawned per outer draw so streams are reproducible and independent.

---

## 7. The headline features — equation divergence, covariate bands, the fifth variance component

### 7.1 The covariate-equation divergence view — divergence *within* a model

v0.1's divergence view overlays *different models* for one patient. v0.7 adds the orthogonal view: for **one** model and one patient, overlay the predicted curve under each **admissible covariate equation**, greying any equation evaluated outside its own validity envelope (James above BMI ≈ 37), and report the spread. The same `compare`/`_divergence` machinery (simulate.py §548/§655) drives it, parameterized over equations instead of models:

```text
$ hypnos covariate-divergence --model propofol.schnider_1998 \
                              --age 50 --weight 130 --height 170 --sex M
propofol.schnider_1998 — covariate-equation divergence (BMI 45.0)
  equation            LBM (kg)   Ce peak (ug/mL)   status
  james_1976            42.1*        (greyed)       OUTSIDE equation envelope (BMI>37) — INVERSION
  janmahasatian_2005    66.8         6.4            in-envelope (substitution; verbatim=false)
  note: the model's DERIVED equation is james_1976; at BMI 45 it has INVERTED.
        this divergence IS the documented Schnider-in-obesity failure mode (v0.1), shown at its source.
```

The teaching point — and the safety point — is that the model-selection risk the divergence view exists to surface sometimes lives *inside* a single model, in a covariate equation no one was looking at. Making it visible at its source is the v0.7 contribution.

### 7.2 The covariate band — value uncertainty propagated

Given a caller-supplied covariate-value distribution, draw the seeded covariate band beside the v0.2/v0.3 bands:

```text
traj.ce_covariate_band   # {5: array, 50: array, 95: array} from weight ~ N(70, 6^2); None if no distribution supplied
```

framed always as "*how much does not-knowing-the-weight-exactly move the prediction?*" — a research question — never "what weight should I enter?"

### 7.3 The fifth variance component — extending the reducible/irreducible split

v0.3 §7 decomposed the total predictive variance into structural / estimation / BSV / residual. v0.7 adds **covariate**, and refines the reducibility rollup:

```
Var_total = Var_structural + Var_estimation + Var_covariate + Var_BSV + Var_residual
            \________________ reducible ________________/      \___ irreducible ___/
              (more models)  (more data)  (better covariate / right equation)
```

`compare()`'s `divergence` dicts gain the fifth share and a finer reducibility breakdown:

```jsonc
"variance_share": { "structural": 0.48, "estimation": 0.12, "covariate": 0.21, "bsv": 0.15, "residual": 0.04 },
"reducibility":   { "reducible": 0.81, "irreducible": 0.19 },
"covariate_split": { "equation_choice": 0.17, "value_uncertainty": 0.04 },   // the two halves of §1, kept distinct
"covariate_band_tier": "B"
```

The punchline a single curve can never give: **for this obese patient, 21% of the predictive uncertainty is the body-composition equation — and almost all of that is *which equation*, not how well we know the weight.** That tells a researcher the lever is "agree on the covariate model," not "weigh the patient more carefully" — a genuinely novel, publishable decomposition the curated equation library makes possible for free.

> **Honesty note (inherited from v0.2 §7 / v0.3 §7).** The decomposition is only as complete as the curated layers. A model with `covariate_sensitivity_status: none` contributes no covariate component; the readout is computed over the eligible subset and **names** what was excluded, never silently dropping it.

---

## 8. Export — projecting the covariate layer to each target's native convention

The covariate equations are *part of the model definition*, so unlike the post-hoc validation of v0.4 they belong inside every export, not as an annotation afterthought.

| Format | v0.7 behavior |
| --- | --- |
| **NONMEM** | The derived-covariate equations become explicit `$PK` (or `$ABBR`/`$MIX`) computations with a comment naming the source equation id — so the control stream *computes LBM the way the model was derived*, not however the reader's template happens to. The most consequential export fix this layer makes. |
| **PharmML + SO** | First-class: derived covariates map to PharmML `Covariate` / `IndividualParameter` transformation blocks with the equation made explicit and the source cited. The durable anchor carries the transform, not just θ. |
| **nlmixr2 / rxode2 · Pumas** | The `_pop` companion's covariate block emits the named equation verbatim with a comment; the covariate-value distribution (if supplied) seeds a parametric covariate resampling, clearly labeled covariate uncertainty, distinct from Ω and from estimate covariance. |
| **SBML (L3v2)** | The derived covariate is an `AssignmentRule` computing LBM/FFM from the input covariates, with the equation id and validity envelope as `hypnos:` RDF (`hypnos:covariateEquation`, `hypnos:covariateValidityEnvelope`). A deterministic consumer sees the right transform; the envelope/failure-mode warning rides as annotation. |
| **TCI-sim JSON** | The `covariate_model` block passes through verbatim (lossless), so a downstream simulator knows *which* LBM equation to use — directly addressing the documented pump-to-pump substitution problem (§1). |

Every covariate-carrying export inherits the universal `clinicalUse = "PROHIBITED"` annotation and, additionally, a `hypnos:covariateBandTier` and a `hypnos:covariateSensitivityStatus`, plus — where an equation is evaluated outside its envelope — the `hypnos:` failure-mode warning, so a consumer cannot reproduce the James inversion unknowingly.

---

## 9. Validation & verification

**Round-trip (automated, CI).** Each exported derived-covariate computation is re-evaluated against the pure-Python equation library and checked to algebraic tolerance (1e-6, v0.1 §6): the exported NONMEM/SBML LBM computation must equal `covariates.evaluate("james_1976", patient)` exactly. An export that quietly computes LBM differently is a CI failure — the strongest possible guard against Trap 1 leaking into artifacts.

**Equation-library self-tests.** Each equation is checked against published worked examples where the source gives them (James's own tabulated LBM values; Janmahasatian's), and its `known_failure_mode` is checked to actually trigger (James LBM must be verified to peak-then-decline above its envelope) — the failure mode is a *tested* property, not an assertion, in the spirit of v0.4 §7.3.

**Internal consistency (`hypnos validate`).** Extends the existing checks (validate.py §30) with: (1) every `covariate_model.derived_inputs[].equation` resolves to a library record; (2) every `used_for` symbol resolves to a real parameter; (3) units are dimensionally consistent between model and equation (Trap 2); (4) `covariate_sensitivity_status` matches contents; (5) every covariate-layer `primary_citation` resolves (the guardrail v0.1 added for `predictive_performance` and v0.2/v0.3 for their blocks).

**Human verification — the covariate checklist.** `review_status` moves `unverified → verified` only when a human confirms, against the source **Methods**, the five traps of §4 as explicit yes/no items (the `_checklist_for` builder, verification.py §54, gains a `covariate_equation` group beside `structural`/`covariate`/`population`/`envelope`/`citation`):

1. Which **named equation** did the authors use for each derived covariate? Is `verbatim: true` truthful? (Trap 1 — the cardinal disambiguation; read the Methods.)
2. Are the equation's **units** (esp. height cm vs. m) consistent with the model? (Trap 2)
3. Is the **sex coding** correct, and the right sex branch applied? (Trap 3)
4. Does the patient fall outside the **equation's own** validity envelope, independently of the disposition envelope? (Trap 4)
5. Are allometric exponents **fixed (theory)** or **fitted**, and labeled as such? (Trap 5)

As in every prior layer, **LLMs assist but never promote**: confirming that Schnider used James and not Janmahasatian is precisely a human-reads-the-Methods line item.

---

## 10. Safety & scope guardrails (non-negotiable, inherited and extended)

The covariate layer is forward-only and changes none of v0.1 §10 / v0.2 §10 / v0.3 §10. It adds guardrails specific to covariate bands and equation choice:

- **No covariate inverse control.** Hypnos will not compute "the weight to enter so the predicted Cp hits *X*," or search over covariates to achieve a target. That is inverse control with a covariate disguise — the regulated-device function the project forbids everywhere. Covariate bands describe a *given* covariate distribution; they never invert one.
- **A covariate band is a statement about input uncertainty, not a per-patient guarantee, and distinct from a prediction band.** The UI and exports label all band kinds distinctly: a covariate ribbon says *"this is how much not knowing the weight exactly moves the curve,"* not *"your patient will be in this range."*
- **The model's derived equation is never silently substituted.** A `verbatim: false` alternative appears only as an explicitly-labeled, cited entry in the *divergence* view — it is never presented as the model's own output. Auto-"fixing" James to Janmahasatian to suppress the inversion would *erase* a documented failure mode the project exists to surface; the inversion is shown, with its citation, as the honest result.
- The `clinicalUse = "PROHIBITED"` annotation remains universal. Adding a covariate band makes the output look more configurable, not more clinical; the guardrails tighten in proportion, by design.

---

## 11. Worked example (illustrative — numbers pending verification)

The intended end-to-end behavior on an obese patient, where the covariate layer dominates:

```text
$ hypnos compare --drug propofol --age 50 --weight 130 --height 170 --sex M \
                 --bands prediction,confidence,covariate --percentile 5,95 --seed 7
included:
  - propofol.eleveld_2018   tier A  covariate: ffm=al_sallami_2015  (in equation envelope)
  - propofol.schnider_1998  tier B  covariate: lbm=james_1976       (OUTSIDE equation envelope — INVERTED)
  - propofol.marsh_1991     tier B  covariate: none (weight-scaled)

covariate-equation note:
  schnider_1998: james_1976 LBM has inverted at BMI 45.0 -> result auto-tiered to D,
                 failure mode CONFIRMED at source (v0.1 known_failure_mode is a covariate-layer property)

effect-site divergence:
  variance share @ t*:  structural 0.41 | estimation 0.10 | covariate 0.34 | bsv 0.12 | residual 0.03
  covariate split:      equation_choice 0.31 | value_uncertainty 0.03
  reducibility:         reducible 0.85 | irreducible 0.15
  reading: in this obese patient ONE THIRD of the spread is the body-composition equation,
           almost all of it WHICH-equation (0.31) not weight-uncertainty (0.03) -> the lever is
           "agree on / fix the covariate model," and Schnider's number here rests on an inverted LBM.
```

The reading: the largest single reducible component is not "which disposition model" but "which body-composition equation" — and Schnider's prediction in this patient is resting on a covariate equation that has gone non-physical. v0.1 already greyed Schnider here; v0.7 explains *why*, at the level of the actual offending object, and quantifies how much of the disagreement it accounts for.

(Every number above is illustrative. Per the project ethos, the equation-library records and `covariate_model` bindings are curated `unverified` and await human Methods-section confirmation before any of this is presented as authoritative.)

---

## 12. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **C0 — Equation library + bindings** | **shipped (0.7.0)** | the `dataset/covariate_equations/` artifact + `covariate_model` schema (§3); generalized `lbm_james` into `covariates.evaluate` with per-equation validity envelopes; curated James/Janmahasatian/Al-Sallami (the equations the curated models are actually derived with — Du Bois/Mosteller/Devine deferred until a model binds one, per the never-invent rule); bound Schnider/Minto→James, Kim→Janmahasatian, Eleveld→Al-Sallami; `validate` resolution + input-availability + sensitivity-status checks; the covariate verification group; the `hypnos covariates` CLI + `covariate_equations.png` figure. | ✅ Every covariate-scaled model names its derived equation; the validator enforces it; James's inversion is a tested property of the library record. |
| **C1 — Equation-divergence view** | **shipped (0.7.0)** | `covariate_divergence` + `hypnos covariate-divergence` overlay one model under each admissible body-size equation (via a backward-compatible kernel override whose verbatim curve reproduces a plain `simulate()` exactly), greying equations outside their own envelope; reuses the v0.1 `_divergence` machinery over equations; the Schnider-obesity failure mode rendered at its covariate source; the `covariate_divergence.png` figure (obese-patient curves + a BMI sweep showing the divergence open at the BMI≤37 boundary). | ✅ The "divergence within a model" view exists and shows the James inversion as the documented failure mode it is. |
| **C2 — Covariate-value bands + the fifth variance component** | **shipped (0.7.0)** | `sample_covariate_vector` + `point_patient`; the seeded covariate band in `simulate(..., bands=["covariate"])` (a distribution-valued covariate collapses to its mean so scalar calls stay byte-identical); the never-invent rule for value distributions; the `covariate` share + `covariate_split` (equation-choice vs value-uncertainty) + refined `reducibility` in `compare`, **gated** so the legacy `bands=True` path keeps its 3-way decomposition; the CLI (`--covariate-band --weight-sd`), the dashboard weight-uncertainty slider + split readout, and the `covariate_band.png` figure. | ✅ `simulate(..., bands=["covariate"], seed=...)` returns reproducible covariate quantiles; the decomposition separates equation-choice from value-uncertainty. |
| **C3 — Exports** | proposed | the derived-covariate computation made explicit in NONMEM `$PK`, PharmML `Covariate` transforms, SBML `AssignmentRule`, TCI-JSON passthrough; round-trip equality of the exported computation vs. the library; `hypnos:` covariate RDF + envelope warnings. | The covariate transform round-trips through the pharmacometric formats; no export can silently recompute LBM. |

C0 alone is a useful, self-contained increment — even before any band is drawn, recording *which* body-composition equation each model was derived with, and validating that exports compute it that way, fixes the single most consequential silent-substitution bug in the field's tooling.

---

## 13. Cheat sheet (target API)

```python
import numpy as np, hypnos
from hypnos.covariates import evaluate as cov_eval
ds = hypnos.load()

m = ds["hypnotics_iv.propofol.schnider_1998"]
m.covariate_model.derived_inputs[0].equation       # "james_1976"
m.covariate_sensitivity_status                      # "declared"
cov_eval("james_1976",  dict(weight=130, height=170, sex="M"))   # LBM via James (will warn: inverted, BMI 45)
cov_eval("janmahasatian_2005", dict(weight=130, height=170, sex="M"))  # FFM, well-behaved

# Forward simulation WITH a covariate band (seeded; still no inverse control)
patient = dict(age=50, weight={"mean": 130, "sd": 8}, height=170, sex="M")
schedule = [("bolus", 0.0, "2 mg/kg"), ("infusion", 0.0, "6 mg/kg/h")]
traj = hypnos.simulate(ds, "hypnotics_iv.propofol.schnider_1998",
                       patient=patient, schedule=schedule, t=np.linspace(0, 60, 600),
                       bands=["prediction", "confidence", "covariate"], percentile=(5, 95), seed=7)
traj.ce_covariate_band     # covariate-value quantiles; None if no covariate distribution supplied
traj.covariate_warnings    # e.g. "james_1976 LBM inverted at BMI 45.0 -> tier D"

# Covariate-equation divergence — divergence WITHIN one model
cd = hypnos.covariate_divergence(ds, "hypnotics_iv.propofol.schnider_1998", patient=patient)
cd.by_equation             # per-equation derived value + curve + in/out-of-equation-envelope status
cd.derived_equation        # "james_1976" — the model's own, flagged inverted here

# The fifth variance component
cmp = hypnos.compare(ds, drug="propofol", purpose="effect_site",
                     patient=patient, schedule=schedule, bands="combined", seed=7)
cmp.divergence["ce"]["variance_share"]    # {structural, estimation, covariate, bsv, residual}
cmp.divergence["ce"]["covariate_split"]   # {equation_choice, value_uncertainty}
```

```bash
hypnos covariate-divergence --model propofol.schnider_1998 --age 50 --weight 130 --height 170 --sex M
hypnos compare  --drug propofol --age 50 --weight 130 --height 170 --sex M --bands covariate --seed 7
hypnos export   --format nonmem --output exports/nonmem/   # $PK now computes LBM via the named equation
hypnos validate                                            # adds equation-resolves + units + verbatim checks
```

---

## 14. Open questions & explicit deferrals

- **Correlated covariates.** Weight and height (hence BMI) are correlated; a covariate-value band that perturbs them independently overstates the spread. Where a caller supplies a covariance we honor it; absent one, we perturb the supplied marginals independently *with a recorded caveat* — the same honest-default stance as v0.2's `omega_block.complete = false`. A curated population covariate-covariance is a possible future refinement, not assumed.
- **Covariate-coefficient estimation uncertainty.** v0.3 §14 deferred SEs on the *covariate-equation coefficients* (the age/weight/LBM slopes) and noted the interaction with this layer. With v0.7's library in place, those coefficient SEs now have a home (the equation record can carry them) and the two propagate together; curating them is a C-phase backfill, gated on source tables exactly as the v0.3 backfill is.
- **Allometric size as the bridge to developmental pharmacology.** The `allometric_size` quantity (weight^¾ for clearance) is curated here as one more named covariate equation, but its *maturation* counterpart — the post-menstrual-age clearance-maturation sigmoid — is the heart of [v0.8](../v0.8/developmental.md); the two specs share the equation library and Trap 5 (fixed vs. fitted exponents) deliberately.
- **Time-varying covariates.** Body weight and organ function can change over a long case (fluid shifts, blood loss); v0.7 propagates covariates fixed at *t = 0*. Time-varying covariate trajectories are a deeper extension deferred until sources and the use-case justify the added machinery.
- **Is equation-choice uncertainty reducible or structural?** It is recorded as *reducible* (a field could converge on one equation), but in practice the "right" body-composition descriptor for PK is itself contested — so the honest, possibly-permanent finding for some patients is "the covariate model is as much a modeling choice as the disposition model," and the §7.1 view is designed to show exactly that, not to resolve it.
