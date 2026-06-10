# Hypnos — design spec (v0.5): breadth — reversal kinetics & special-population envelopes

**Build the two coverage subsystems v0.1 declared in scope but deferred: (A) `reversal` — the kinetics of antagonists and encapsulators (sugammadex, neostigmine, naloxone, flumazenil) modeled as *mechanistic interactions* with the agonist they reverse, not as standalone PK; and (B) `special_populations` — disease-state and special-population applicability (hepatic / renal / cardiac impairment, pregnancy, critical illness) encoded as *covariate-envelope extensions and cited adjustment annotations*, Tier-D by construction unless a model was actually fitted in that group — so the dataset can say, honestly, where a model has any standing in an organ-failure patient and where it is pure extrapolation.**

> **Status:** **S0 implemented (the organ-function envelope); reversal (B0–B1) and disease modifiers (S1) proposed.** The physiological envelope now speaks: `applicability_envelope` carries optional organ-function ranges + a cited `organ_tolerance[]` block, `Envelope.organ_check` / `evaluate_safety` grey every model with no standing in a declared hepatic/renal/cardiac/albumin impairment (Tier-D + named extrapolation), the remifentanil esterase exception is encoded with citations (Dershwitz 1996 / Hoke 1997), and the overlay is live in `compare`, the CLI, and the dashboard. The `reversal` subsystem (Part A) and the `disease_state_modifier` records (S1) remain proposed. Additive over [v0.1](../v0.1/spec.md)–[v0.4](../v0.4/external_validation.md): no existing record changes meaning; both subsystems are new `subsystem`/`purpose` values plus optional blocks. Dataset/schema target `0.5.0`. This spec builds the two items v0.1 named and held back: the `reversal` subsystem (v0.1 §3 lists it; only `rocuronium_tof_pd` exists today, no reversal agent is curated) and disease-state PK, which v0.1 §2 explicitly placed out of v0.x scope as *"a future tier-D research envelope."* v0.5 designs that envelope so the extrapolation is *labeled and bounded*, never silent.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input.* Breadth is the most extrapolation-prone kind of growth, so v0.5 leans hardest on the project's two oldest guardrails — the **applicability envelope** (v0.1) and the **never-invent rule** (v0.2): a reversal agent simulates only the trajectory of a *given* dose, never the reversal dose; a disease-state adjustment is applied only with a citation, a loud Tier-D label, and a warning, or not at all.

---

## Part A — the `reversal` subsystem

### A1. The problem — antagonism is an interaction, not a solo drug

A reversal agent has no clinically interesting trajectory *alone*; its entire purpose is to act on **another drug already in the patient**. Modeling it as a standalone PK record (the way Hypnos models propofol) misses the physics. Reversal is mechanistically one of three *interactions*, each distinct from the propofol–remifentanil *synergy* surface Hypnos already curates (`propofol_remifentanil_surface`, the Greco kernel, reference.py §494):

1. **Encapsulation / binding** — **sugammadex** is a γ-cyclodextrin that binds *free* rocuronium (and vecuronium) in plasma at ~1:1 stoichiometry, forming an inert complex. It does not act at a receptor; it *removes the drug from the active pool*. The free (effect-driving) NMB concentration drops, the drug redistributes out of the effect compartment, and train-of-four recovers. This is a **PK-level** mechanistic coupling — a binding reaction layered onto the existing NMB PK.
2. **Competitive receptor antagonism** — **naloxone** (opioid) and **flumazenil** (benzodiazepine) compete with the agonist at the *receptor*. Neither changes the agonist's PK; both shift the agonist's apparent potency. This is a **PD-level** coupling — a concentration-dependent rightward shift of the agonist's Hill curve (an increase in apparent Ce50).
3. **Enzyme inhibition (indirect)** — **neostigmine** inhibits acetylcholinesterase, raising synaptic ACh to *out-compete* the NMB. It is indirect, has a pronounced **ceiling** (useless at deep block), and carries muscarinic effects requiring co-administration of an antimuscarinic — a documented context, not a free parameter.

These three mechanisms need three small kernels, not one. None is inverse control (§A6).

### A2. What we curate

For each reversal agent: its own disposition PK (so the antagonist's own concentration is known), **plus** a `reversal_mechanism` block linking it to the agonist model(s) it acts on and the parameters of the coupling:

- **encapsulation**: the target NMB model id(s), binding stoichiometry, and association/dissociation (`k_on`/`k_off` or an assumed-irreversible `K_d`) where published; otherwise the commonly-modeled instantaneous-1:1 simplification, *labeled as a simplification*.
- **competitive_receptor**: the target agonist PD model id(s) and the antagonist potency (the `K_i` / apparent-Ce50-shift relation).
- **enzyme_inhibition**: the indirect-effect parameters and the explicit **ceiling** and **depth-of-block dependence** as a first-class `known_failure_mode`.

### A3. Schema extension

New `subsystem: "reversal"` and `purpose: "reversal"`. New `$def` `reversal_mechanism`:

```jsonc
"reversal_mechanism": {
  "type": "encapsulation",                 // encapsulation | competitive_receptor | enzyme_inhibition
  "targets": ["nmb_agents.rocuronium.wierda_1991"],   // the agonist model(s) reversed — must resolve
  "stoichiometry": { "ratio": 1.0, "basis": "molar" },
  "binding": { "k_on": null, "k_off": null, "kd": null, "approximation": "instantaneous_1to1" },
  "depth_dependence": {                     // sugammadex dose depends on block depth; encode the regimes as ENVELOPE
    "regimes": [
      { "block_depth": "moderate (TOF reappearance)", "context": "..." },
      { "block_depth": "deep (PTC 1-2)",              "context": "..." }
    ]
  },
  "tier": "C",
  "extraction": { "review_status": "unverified", "source_locator": "..." },
  "primary_citation": "..."
}
```

For `competitive_receptor`, `binding` is replaced by `antagonism: { ki, ce50_shift_model }`; for `enzyme_inhibition`, by `indirect: { ceiling, depth_floor }`. `targets` is validated to resolve to existing model ids (the same resolve-or-error guardrail v0.1 added for citations). The off-target / failure contexts (renarcotization, ceiling) are curated as standard `known_failure_modes`, reusing the existing `$def`.

### A4. Reference kernels

A new `reversal.py` beside the existing kernels, each primitive forward-only and round-trip-validated like every other (v0.1 §6):

- **`reversal.encapsulation_ode(nmb_micro, sugammadex_conc, stoich, binding, t)`** — augments the NMB's linear PK with a binding reaction that depletes the *free* central-compartment NMB; the existing matrix-exponential solver (reference.py §147) handles the linear part, the binding term integrates numerically (reference.py §303 `simulate_numeric` already provides the ODE path). Output: free-NMB and complex trajectories → the existing TOF PD model (`rocuronium_tof_pd`) consumes the free concentration unchanged.
- **`reversal.competitive_shift(agonist_ce, antagonist_ce, ki)`** — returns the agonist's apparent Ce50 multiplier `(1 + C_antag/Ki)`, fed into `sigmoid_emax` (reference.py §357). The agonist PK is untouched; only its PD potency shifts.
- **`reversal.indirect_inhibition(...)`** — the neostigmine ceiling/indirect-effect form, with the depth floor enforced as a hard envelope.

### A5. Tiers & the renarcotization failure mode

Reversal records tier on the usual A–D spine, with one mechanism-specific honesty requirement that is itself a safety feature:

- **Renarcotization / recurarization is a mandatory `known_failure_mode`** where applicable. **Naloxone's** half-life is shorter than that of most opioids it reverses, so the antagonist can wash out *before* the agonist — re-sedation/re-narcotization. The kernel makes this *visible*: simulate naloxone reversing a long-acting opioid forward and the predicted effect transiently recovers, then **relapses** as naloxone clears. Encoding and rendering this relapse is the whole point — it is the kind of true, cited, safety-relevant behavior the project exists to surface. Flumazenil carries the analogous re-sedation caveat; neostigmine carries the deep-block ceiling.
- **Stoichiometric simplifications are labeled.** An instantaneous-1:1 sugammadex approximation is Tier-C-at-best until binding kinetics are curated and validated; the `approximation` field names it so no one mistakes the simplification for the measured mechanism.

### A6. Safety (Part A)

Reversal is where the dosing temptation is most acute — *"how much sugammadex / naloxone do I give?"* The line is identical to v0.1 §10 and absolute:

- **Hypnos simulates the trajectory of a *given* reversal dose; it never computes the reversal dose.** No "dose to achieve TOF ≥ 0.9," no "naloxone dose to reverse this fentanyl." That is inverse control wearing an antidote's coat.
- The **renarcotization** behavior is rendered as honest *forward* output (a relapse curve), never as a re-dosing schedule.
- `clinicalUse = "PROHIBITED"` is universal here as everywhere, and reinforced because antidotes are emergency drugs.

---

## Part B — the `special_populations` subsystem

### B1. The problem — where every model is, honestly, extrapolating

v0.1's envelope machinery already greys out a model when a patient's *demographic* covariates (age, weight, BMI) fall outside the derivation range. But the highest-stakes extrapolations are *physiological*: hepatic failure, renal failure, low cardiac output, pregnancy, critical illness. Almost no anesthetic PK model was fitted in these groups, yet they are exactly the patients in whom getting the dose wrong is most dangerous. Today Hypnos is silent on them — and silence reads as "fine." v0.5 makes the silence *speak*: it adds the organ-function dimensions to the envelope so an out-of-physiological-range simulation is auto-tiered to D and warned, exactly as a too-high-BMI simulation already is.

### B2. What we curate — two record kinds, never blurred

1. **Validated special-population models.** Where a population model was *actually fitted* in a disease group (e.g. a propofol PK study in cirrhosis), it is a normal Hypnos model record, tiered on its own merits, with the disease group named in its envelope/populations. These earn their tier.
2. **Disease-state adjustment annotations.** Far more common: the literature reports a *documented adjustment* — "in Child-Pugh C, clearance is reduced by ~X%" — without a full refitted model. These are curated as `disease_state_modifier` annotations attached to an existing model: cited, tiered, and applied **only** opt-in, **always** Tier-D, **always** warned. A modifier is an advisory, cited correction — never a silent reparameterization, and never present without a citation.

> The single most important distinction in this subsystem: **a fitted-in-disease model is evidence; an adjustment factor is an annotation.** Conflating them is the failure mode this design exists to prevent.

### B3. The physiology the envelope must carry

- **Hepatic impairment (Child-Pugh A/B/C):** reduced clearance for hepatically-cleared drugs (most IV hypnotics/opioids); remifentanil is the notable exception (esterase metabolism, organ-independent) — encoded so the envelope does *not* wrongly penalize it.
- **Renal impairment (CrCl / eGFR):** matters where active metabolites accumulate (e.g. morphine-6-glucuronide); encoded per-drug as a metabolite-accumulation caveat, not a blanket clearance scaling.
- **Cardiac (low cardiac output / low ejection fraction):** slowed front-end kinetics — reduced apparent central volume and inter-compartmental clearance, higher and faster peak after a bolus. A documented, cited failure mode for bolus dosing.
- **Pregnancy:** increased volume of distribution and cardiac output, altered protein binding — a distinct physiological envelope.
- **Critical illness / ICU:** capillary leak (expanded volume), **hypoalbuminemia** → raised *free fraction* of highly-bound drugs (propofol ~98% bound, fentanyl) → effect from a given *total* concentration is under-predicted. This is the protein-binding story, and it is a genuine, citable failure mode for the binding-sensitive drugs.

### B4. Schema extension

New `subsystem: "special_populations"`. The `applicability_envelope` `$def` gains optional organ-function dimensions (reusing the `range` `$def`): `crcl_ml_min`, `child_pugh` (categorical), `albumin_g_dl`, `ejection_fraction_pct`. New `$def` `disease_state_modifier`:

```jsonc
"disease_state_modifier": {
  "condition": { "dimension": "child_pugh", "value": "C" },
  "affected_parameter": "Cl",
  "adjustment": { "type": "scale_factor", "value": 0.6, "equation": null },  // or a covariate-style equation
  "binding_sensitive": false,            // true => the modifier acts via free-fraction change (§B3)
  "evidence_tier": "D",                  // adjustments default to D; a fitted model is a different record
  "applied_by_default": false,           // OPT-IN ONLY — never silently reparameterizes a model
  "caveat": "Single small cohort; directionally supported, magnitude uncertain.",
  "primary_citation": "..."
}
```

`evaluate_safety` (simulate.py §211) and `Envelope.check` (models.py §161) extend to the new dimensions: a patient with a curated `child_pugh` / `crcl_ml_min` / `albumin_g_dl` outside a model's range is greyed and auto-tiered to D, identical to today's BMI handling. A `disease_state_modifier` is applied only when the caller opts in, and applying one **forces** the result to Tier-D with a warning naming the modifier and its citation.

### B5. Tiers & the never-invent rule (Part B)

- **Disease modifiers are Tier-D by construction.** They describe a documented *direction and rough magnitude*, not a validated model. The median line, once a modifier is applied, is labeled Tier-D — you cannot get an A-looking number out of a disease adjustment.
- **No modifier without a citation** (the never-invent rule). Hypnos will not "scale clearance by 30% for cirrhosis" on general physiological intuition; either a source reports it or the field is absent (an honest gap, rendered as "no curated disease adjustment").
- **No borrowing across drugs.** A hepatic adjustment curated for one drug is never auto-applied to another, even within a class — the same refusal-to-synthesize v0.2 applied to BSV.

### B6. Headline feature — the envelope-shrinkage / extrapolation overlay

The divergence view (`compare`, simulate.py §589) gains a special-population overlay that makes the extrapolation *visible*: enter a cirrhotic or low-EF virtual patient and watch the eligible-model set **shrink** as physiological envelopes are violated, with every greyed model named and reasoned (exactly as the obese-patient panel already greys Schnider). Where disease modifiers exist, they fan the prediction out as an explicitly Tier-D *extrapolation band*, labeled as such. The honest message the view delivers: *"in this organ-failure patient, model X has some standing and everything else is extrapolation — here is how far the extrapolations spread, and not one of them is validated here."* That is a research instrument for exactly the patients the literature is quietest about.

### B7. Safety (Part B)

This is the most extrapolation-heavy subsystem in Hypnos, so it carries the strongest guardrails:

- **Modifiers off by default, opt-in, always Tier-D, always warned, always cited** (§B5).
- **The envelope speaks in organ-failure space** — silence is replaced by an explicit greyed-out + Tier-D + reason, so a disease-state simulation can never *look* validated.
- **Protein-binding caveats are surfaced, not silently modeled** — for binding-sensitive drugs the free-fraction shift is flagged as a failure mode with its citation, not folded invisibly into a clearance number.
- `clinicalUse = "PROHIBITED"` is universal and, for this subsystem, doubly load-bearing: these are the patients where a misread extrapolation would do the most harm.

---

## C. Validation & verification (both subsystems)

- **Round-trip (CI), unchanged discipline.** Every reversal kernel and every disease-modifier application is re-simulated through the exported artifact and checked against the pure-Python reference within tolerance (algebraic 1e-6 / ODE 1e-4, v0.1 §6).
- **Literature validation where supported.** Sugammadex recovery-time profiles and naloxone renarcotization case data, where a source publishes them, are compared to the kernel output and the agreement recorded — the breadth analog of v0.1's check-against-published-sims.
- **Human verification — extended checklist.** The `_checklist_for` builder (verification.py §54) gains a `reversal` group (target-model resolves; stoichiometry/binding scale; renarcotization/ceiling failure mode present) and a `special_population` group (modifier has a citation; `evidence_tier` is D unless a fitted model; `applied_by_default` is false; binding-sensitivity flagged correctly). **LLMs assist but never promote** (v0.1 §9), and for disease modifiers the human specifically confirms the *direction and magnitude* against the cited source.
- **`hypnos validate` additions:** `reversal_mechanism.targets` and `disease_state_modifier.primary_citation` resolve; every disease modifier with `evidence_tier != "D"` is justified by a linked fitted-in-disease model (else flagged); `applied_by_default` is never `true` without explicit reviewer sign-off recorded in `extraction`.

## D. Export

Both subsystems flow through the existing exporters with mechanism-appropriate handling: reversal records export their own PK as normal models, plus the `reversal_mechanism` as a `hypnos:` RDF block linking to the target model (SBML/NONMEM core cannot express a cross-model binding reaction natively; the RDF carrier is lossless and the TCI-JSON passthrough keeps it). Disease modifiers export as `hypnos:diseaseStateModifier` annotations carrying their citation, tier, and the opt-in/Tier-D flags — so a downstream consumer cannot apply one without seeing that it is an unvalidated, opt-in extrapolation. Every artifact keeps the universal `clinicalUse = "PROHIBITED"` annotation.

## E. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **B0 — Sugammadex encapsulation** | proposed | `reversal` subsystem + schema; sugammadex record + `encapsulation_ode` kernel reversing the existing rocuronium PK→TOF stack; renarcotization/recovery as forward output. | The canonical reversal pairing (rocuronium + sugammadex) simulates end-to-end, recovery rendered, no dose computed. |
| **B1 — Competitive & indirect antagonists** | proposed | naloxone (opioid) + flumazenil (benzodiazepine) competitive-shift kernel; neostigmine indirect/ceiling kernel; renarcotization & deep-block-ceiling failure modes. | The three reversal mechanisms are each modeled and each carries its safety-relevant failure mode. |
| **S0 — Organ-function envelope** | ✅ shipped | `applicability_envelope` gains organ-function ranges + cited `organ_tolerance[]`; `Envelope.organ_check`/`evaluate_safety` grey + Tier-D every model with no standing in a declared hepatic/renal/cardiac/albumin impairment (named extrapolation), remifentanil's esterase exception encoded (Dershwitz/Hoke); overlay live in `compare`, CLI (`--child-pugh/--crcl/--albumin/--ejection-fraction`), and dashboard. Staging cut-points are definitional (KDIGO/Child-Pugh), named in code. | An organ-failure virtual patient correctly greys out and Tier-D's the models with no standing; remifentanil keeps standing with a citation + metabolite caveat. |
| **S1 — Disease modifiers + binding** | proposed | `disease_state_modifier` records (opt-in, Tier-D, cited); the free-fraction / protein-binding failure mode for binding-sensitive drugs; the explicit Tier-D extrapolation band. | Documented disease adjustments are curated, applied only with a citation + warning, and the binding story is surfaced. |

B0 alone is a useful, self-contained increment: it gives Hypnos its first reversal pairing — the most-used one in modern practice — with the renarcotization honesty built in. S0 alone is independently valuable: even before any modifier is curated, making the *physiological* envelope speak turns Hypnos's silence on organ-failure patients into an explicit, honest "this is extrapolation."

## F. Open questions & explicit deferrals

- **Full binding kinetics vs. instantaneous 1:1.** Curating sugammadex's true `k_on`/`k_off` (vs. the common instantaneous-1:1 simplification) awaits source-table confirmation, exactly as the v0.2 BSV backfill does — the kernel supports both; the simplification is labeled until the kinetics are verified.
- **Active metabolites as their own compartments.** Renal accumulation of active metabolites (morphine-6-glucuronide, norketamine) is currently a *caveat*; modeling them as explicit metabolite compartments is a deeper extension deferred to a later breadth pass.
- **Disease severity as a continuum.** Child-Pugh is categorical here; continuous severity scaling is deferred until sources support it, to avoid implying precision the evidence lacks.
- **Interaction stacking.** Reversal-on-top-of-synergy (e.g. sugammadex in a patient also on a propofol–remifentanil surface) composes in principle; the worst-tier-wins rule already covers the labeling, but the validated stacking of three mechanisms at once is deferred.
- **Lipid-rescue and other therapeutics are out of scope** — Hypnos models the antagonist's reversal kinetics, not the management of toxicity; that boundary is reinforced in [v0.6](../v0.6/local_anesthetics.md).
