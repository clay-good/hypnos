# Hypnos — design spec (v0.8): the developmental-pharmacology subsystem

**Make the pediatric and neonatal envelope *speak* the way v0.5 made the organ-failure envelope speak — by curating the two mechanistic components developmental pharmacology actually turns on: (1) theory-based **allometric size** scaling (clearance ∝ weight^¾, volume ∝ weight¹) and (2) the **post-menstrual-age clearance-maturation** sigmoid (the Anderson–Holford function) — and by drawing the same hard line v0.5 drew between *evidence* and *annotation*: a model **fitted in children** (Kataria, Paedfusor) earns its tier, while an adult model carried into a neonate by allometry-plus-maturation is an explicit, bounded, **Tier-D extrapolation** — never a silent linear-per-kilogram rescaling, which is exactly the mental model that overdoses neonates and the one this subsystem exists to make visibly wrong.**

> **Status: proposed — design spec only; nothing in this document is implemented yet.** Additive over [v0.1](../v0.1/spec.md)–[v0.7](../v0.7/covariate_uncertainty.md): no existing record changes meaning; the developmental block is a new optional structure plus new envelope dimensions, and "not curated" renders as an explicit gap, never as a fabricated maturation curve. Dataset/schema target `0.8.0`. This subsystem builds out the item v0.1 named in its phased roadmap and only stubbed: v0.1 Phase C promised "pediatric models (Kataria/Paedfusor) with **explicit Tier-D extrapolation labeling**," and the codebase already carries the seed — `_classify_age_extrapolation` (simulate.py §174) flags out-of-age-range use today, and `kataria-1994-propofol-pediatric` / `absalom-2005-paedfusor` exist as citation records. v0.8 turns that flag into a *mechanistic* developmental layer: not merely "this patient is younger than the derivation range" but "here is the allometric + maturation extrapolation, here is how far it spreads, and here is why none of it is validated in this neonate." It is the developmental sibling of v0.5's `special_populations` — same evidence-vs-annotation discipline, same envelope-speaks design, a different physiological axis.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input.* Pediatric dosing is the single most extrapolation-prone and highest-consequence corner of anesthetic pharmacology, so v0.8 leans hardest of any breadth spec on the project's two oldest guardrails — the **applicability envelope** (v0.1) and the **never-invent rule** (v0.2): a maturation curve is applied only with a citation, a loud Tier-D label, and a warning, or not at all; and the dangerous-but-intuitive "just scale the adult dose per kilogram" extrapolation is rendered as honest, *visibly wrong* forward output rather than offered as a number.

---

## 1. The problem this subsystem solves — children are not small adults, and silence reads as "fine"

v0.1's envelope greys a model when a patient's demographic covariates fall outside the derivation range, and `_classify_age_extrapolation` (simulate.py §174) already names an out-of-age-range simulation. But a bare "out of range" flag is the *beginning* of the developmental story, not the end. Two mechanistic facts make pediatric PK qualitatively different from adult PK, and neither is captured by greying-out alone:

1. **Size is non-linear in weight.** Physiological rates do not scale linearly with body weight. The theory-based allometric relations — **clearance ∝ weight^0.75**, **volume ∝ weight^1.0** — mean that a 3 kg neonate does *not* clear a drug at 3/70 of an adult's rate; per-kilogram, small bodies clear *faster*. A naive linear-per-kg adult rescaling therefore **under**-doses on volume and mis-states clearance — and the most dangerous version, "scale the adult mg/kg down linearly," is still a real-world mental model.
2. **Clearance organs are immature at birth, and mature over months.** Even after allometric size is accounted for, hepatic and renal elimination pathways are *functionally immature* in the neonate and mature along a well-described sigmoidal trajectory in **post-menstrual age** (PMA = gestational age + postnatal age). The **Anderson–Holford maturation function** captures this:

   ```
   MF(PMA) = PMA^Hill / (TM50^Hill + PMA^Hill)
   ```

   A term neonate (PMA ≈ 40 weeks) may have only a fraction of the size-adjusted adult clearance; ignoring maturation **over**-doses the neonate exactly where the margin for error is smallest. Maturation acts in the *opposite* direction to the allometric per-kg speed-up, and the two must be composed, not confused.

Today Hypnos is mechanistically silent on both. It can grey an adult model in a child, but it cannot say *"here is what allometry-plus-maturation would extrapolate, and here is the spread of that extrapolation"* — and silence, as v0.5 §B1 observed, reads as "fine." v0.8 makes the silence speak in developmental space, with the same move v0.5 used for organ failure: an explicit, mechanistic, **labeled Tier-D** extrapolation instead of a blank.

> The clinical stakes are the sharpest in the dataset. The literature is quietest exactly where the patients are smallest, the doses least forgiving, and the temptation to extrapolate strongest. A research instrument that shows *how far the extrapolations spread and that not one is validated here* is most valuable precisely in this population.

The headline deliverable (§6) is the machine-readable, mechanistically-honest answer to: **in this neonate/infant, which models have actual pediatric standing, how does the allometric + maturation extrapolation of the adult models spread, and how visibly wrong is the linear-per-kg shortcut everyone is tempted by?**

---

## 2. What we curate — two record kinds, never blurred (the v0.5 discipline, developmental axis)

The cardinal distinction is v0.5 §B2's, transposed:

1. **Fitted-in-children models — evidence.** Where a population model was *actually fitted* in a pediatric cohort (propofol: **Kataria 1994**, children 3–11 y; **Paedfusor/Absalom 2005**, pediatric TCI; the Marsh pediatric variant), it is a normal Hypnos model record, tiered on its own merits, with the pediatric age band and its size/maturation structure named in its envelope. These earn their tier. A Kataria simulation in a 6-year-old is in-envelope and tiered normally; the developmental machinery only *labels* what is already a fitted model.
2. **Developmental extrapolation annotations — annotation.** Far more common: an adult model carried into a child by an explicit allometric **size model** and, where the source supports it, a **maturation model**. These are curated as a `developmental_model` block attached to an existing adult model: cited, **Tier-D by construction**, applied **only** opt-in, **always** warned. An extrapolation is a documented, cited *mechanism for extending reach*, never a silent reparameterization and never present without a citation for its components.

> The single most important rule in this subsystem: **a fitted-in-children model is evidence; an allometric-plus-maturation extrapolation is an annotation.** Conflating them — letting an adult model's neonatal extrapolation wear a fitted model's tier — is the failure mode this design exists to prevent, and it is the exact analog of v0.5's "a fitted-in-disease model is evidence; an adjustment factor is an annotation."

---

## 3. The developmental physiology the envelope and the blocks must carry

- **Allometric size (universal, theory-based).** Clearance and inter-compartmental clearance scale with weight^¾; volumes with weight¹. The exponents are *theory-based and fixed* by default (West/Brown allometry as applied to PK by Anderson & Holford); a model that *fitted* its exponent records that instead (the v0.7 Trap 5 — fixed vs. fitted — applies here verbatim, and the two specs share the field).
- **Maturation (the neonatal/infant load-bearing fact).** Size-adjusted clearance is further multiplied by the PMA-driven maturation factor `MF(PMA)` (§1). `TM50` (the PMA at half-maximal maturation, in weeks) and the `Hill` coefficient are the curated parameters; both are drug- and pathway-specific and curated only where a source publishes them.
- **Body composition shifts with age.** Neonates have higher total body water and lower fat and muscle mass; the fat-free-mass equation appropriate across the pediatric range is **Al-Sallami (2015)**, *not* the adult Janmahasatian — which is exactly why v0.7's covariate-equation library is a prerequisite: the right size descriptor for a child is itself an equation-choice question.
- **Protein binding differs in the neonate** (lower albumin and α1-acid glycoprotein, higher free fraction of highly-bound drugs) — the same free-fraction failure mode v0.5 §B3/§S1 and v0.6 LA3 curate, here flagged as a developmental caveat rather than re-modeled.
- **Post-menstrual vs. postnatal vs. chronological age.** Maturation is driven by **PMA**, not chronological age; a 2-month-old born at 28 weeks (PMA ≈ 36 weeks) is developmentally *behind* a 2-month-old born at term (PMA ≈ 48 weeks). Collapsing these is the cardinal pediatric trap (§4).

---

## 4. Schema extension & the developmental transcription traps

New `subsystem` value `developmental` is *not* introduced — developmental records are existing models plus blocks, exactly as v0.5 disease modifiers attach to existing models. The `applicability_envelope` `$def` gains optional developmental dimensions (reusing the `range` `$def`, models.py §40): `pma_weeks`, `postnatal_age_days`, `gestational_age_weeks`, alongside the existing `age_years`. New `$defs`: `size_model`, `maturation_model`, and the attachment block `developmental_model`.

```jsonc
"developmental_model": {
  "size": {                                  // allometric size — the universal part
    "cl_exponent": 0.75, "v_exponent": 1.0,
    "reference_weight_kg": 70,
    "exponent_basis": "theory_fixed",        // theory_fixed | fitted  (v0.7 Trap 5; "fitted" needs an SE home)
    "size_descriptor": "ffm:al_sallami_2015",// resolves to the v0.7 covariate-equation library
    "tier": "C", "primary_citation": "anderson-holford-2008"
  },
  "maturation": {                            // PMA-driven clearance maturation — neonatal/infant load-bearing
    "function": "sigmoidal_pma",
    "tm50_weeks": 47.7, "hill": 3.4,
    "driver": "pma_weeks",                   // MUST be PMA, not chronological age (Trap 1)
    "affected_parameter": "Cl",
    "tier": "D", "primary_citation": "..."   // maturation params are rarely Tier > C for anesthetic drugs
  },
  "extrapolation_basis": "allometry_plus_maturation",  // allometry_only | allometry_plus_maturation | fitted_pediatric
  "evidence_tier": "D",                      // extrapolation defaults to D; a fitted-in-children model is a DIFFERENT record
  "applied_by_default": false,               // OPT-IN ONLY — never silently reparameterizes an adult model
  "caveat": "Adult model extrapolated by allometry + maturation; not validated below the derivation age.",
  "primary_citation": "..."
}
```

`evaluate_safety` (simulate.py §211) and `Envelope.check` (models.py §278) extend to the new dimensions: a patient with a curated `pma_weeks` / chronological age below a model's range is greyed and auto-tiered to D, identical to today's BMI and v0.5's organ-function handling. A `developmental_model` extrapolation is applied only when the caller opts in, and applying one **forces** Tier-D with a warning naming the size/maturation components and their citations.

**Developmental transcription traps (additions to the v0.2/v0.3/v0.7 trap lists):**

1. **PMA vs. postnatal vs. chronological age (the cardinal pediatric sin).** A maturation function driven by the wrong age clock is silently, dangerously wrong — most so in the premature neonate. The `driver` field is mandatory and validated to be a PMA-class clock whenever a `maturation` block is present; the verifier confirms which clock the source used.
2. **Allometric exponent: fixed vs. fitted, and on the right parameter.** A ¾ exponent belongs on clearances, a 1.0 on volumes; swapping them, or treating a fitted exponent as the theoretical one, misstates scaling. Shared with v0.7 Trap 5; `exponent_basis` records it.
3. **Reference weight / reference individual.** Allometry is *relative to a reference* (commonly 70 kg adult); an extrapolation that drops or mis-states the reference weight rescales everything. Curated explicitly.
4. **Maturation applied to the wrong pathway.** Renal and hepatic pathways mature on different trajectories; a single `MF` applied to a drug whose clearance is multi-pathway is a simplification — labeled as such, never presented as exact.
5. **Body-composition equation mismatch.** Using the adult Janmahasatian FFM in a child instead of Al-Sallami is a v0.7 Trap-1 substitution with developmental teeth; the `size_descriptor` must resolve to a pediatric-valid equation across the target age range.

These become explicit yes/no items in the verification checklist (§9).

---

## 5. Reference kernels & tiers

A new `developmental.py` beside the existing kernels, each primitive forward-only and round-trip-validated (v0.1 §6). None is inverse control (§10).

- **`developmental.allometric_scale(theta_ref, weight, exponents, reference_weight) -> theta`** — scale the reference (adult) disposition parameters to the patient's size by the curated exponents, applied *before* the structural→micro conversion (`MicroParams.from_volumes_clearances`, reference.py §62), so the matrix-exponential solver (reference.py §147) and exporters keep their single canonical representation — the same composition point v0.7's covariate sampling and v0.2's `sample_individual` use.
- **`developmental.maturation_factor(pma_weeks, tm50, hill) -> float`** — the Anderson–Holford sigmoid, a verbatim transcription of the `function` field, multiplying the size-scaled clearance. Returns 1.0 (no effect) only when explicitly mature; never silently defaults to 1.0 when maturation is *uncurated* (that is a gap, not a mature patient).
- **Tiers on the usual A–D spine, with developmental floors:**
  - **Allometric size extrapolation is Tier-D by construction** when carrying an adult model below its derivation age. Theory-based allometry is a strong, cited *direction*, but an adult model's *parameters* extrapolated to a neonate are not validated there — and the project does not let a strong mechanism launder an unvalidated extrapolation into a higher tier.
  - **Maturation parameters rarely exceed Tier-C** for anesthetic drugs (sparsely published, drug-specific); the design *embraces* that, exactly as v0.6 embraced Tier-C toxicity thresholds — a Tier-D extrapolation band shown honestly is the point.
  - **A fitted-in-children model keeps its own tier** (§2); the developmental machinery only annotates it.
- **The never-invent rule (non-negotiable), restated.** Hypnos does not invent a `TM50`, borrow a sibling drug's maturation curve, or assume allometric scaling stands in for maturation. A missing maturation block is a true statement ("no maturation model is curated for this drug"); a borrowed one is a lie with a sigmoid. Where only size is curated and maturation is not, the extrapolation is labeled `allometry_only` with the explicit caveat that **un-modeled maturation means the neonatal clearance is over-stated** — the failure mode named, not hidden, exactly as v0.6 LA names the binding-saturation gap.

---

## 6. Headline feature — the developmental extrapolation overlay (and the visibly-wrong linear shortcut)

The divergence view (`compare`, simulate.py §655) gains a developmental overlay that makes the extrapolation, and its danger, *visible* — the v0.5 §B6 envelope-shrinkage overlay, transposed to age:

- Enter a neonate/infant virtual patient (PMA + weight). The eligible-model set **shrinks** to those with actual pediatric standing (Kataria/Paedfusor in their bands), with every greyed adult model named and reasoned, exactly as the obese panel greys Schnider.
- For the greyed adult models, where a `developmental_model` exists, fan them out as an explicitly **Tier-D extrapolation band**, labeled as such — `allometry_plus_maturation` where maturation is curated, `allometry_only` (with the over-dose caveat) where it is not.
- Overlay, as a deliberately-labeled reference line, the **naive linear-per-kg adult rescaling** — and show it diverging from both the fitted pediatric models and the allometry+maturation extrapolation. The teaching point is the safety point: *the linear shortcut is visibly, mechanistically wrong*, over-dosing the neonate by ignoring maturation while mis-stating volume by ignoring allometry. Rendering that divergence is the whole purpose — a true, cited, safety-relevant behavior surfaced as honest forward output, never as a corrected dose.

The honest message the view delivers: *"in this neonate, model X (fitted) has standing; everything else is extrapolation; here is how far allometry-plus-maturation spreads; and here is how badly the per-kilogram shortcut misses — not one of these is validated in this patient."* That is a research instrument for exactly the population the literature is quietest about.

---

## 7. Export

Developmental records flow through the existing exporters with mechanism-appropriate handling, reusing v0.7's covariate-transform export path (the size descriptor *is* a covariate equation):

| Surface | v0.8 behavior |
| --- | --- |
| **NONMEM** | the allometric size scaling and the maturation `MF(PMA)` become explicit `$PK` computations with comments naming the source equations, so the control stream extrapolates the way the curated mechanism says — not however a template guesses. The maturation driver is emitted as PMA, with a comment forbidding its substitution by chronological age. |
| **PharmML + SO** | size and maturation map to PharmML covariate-transformation blocks (allometry) and a parameter-level maturation function, both cited; the durable anchor carries the developmental transform. |
| **nlmixr2 / rxode2 · Pumas** | the `_pop` companion's `@pre`/covariate block emits allometry + maturation verbatim with comments; clearly labeled extrapolation, distinct from the fitted-pediatric records. |
| **SBML (L3v2)** | size scaling and `MF(PMA)` as `AssignmentRule`s, with the developmental envelope, the Tier-D flag, and the `allometry_only` over-dose caveat as `hypnos:` RDF (`hypnos:developmentalModel`, `hypnos:maturationFunction`). |
| **TCI-sim JSON** | the `developmental_model` block passes through verbatim (lossless), carrying the opt-in/Tier-D flags so a downstream simulator cannot apply a neonatal extrapolation without seeing it is unvalidated. |

Every developmental artifact inherits the universal `clinicalUse = "PROHIBITED"` annotation and, additionally, a `hypnos:developmentalExtrapolation` flag and `hypnos:evidenceTier`, doubly load-bearing here because these are the patients in whom a misread extrapolation does the most harm.

---

## 8. Validation & verification

- **Round-trip (CI), unchanged discipline.** Every allometric/maturation kernel and its exported computation are re-simulated and checked against the pure-Python reference within tolerance (algebraic 1e-6 / ODE 1e-4, v0.1 §6). The exported maturation factor must equal `developmental.maturation_factor(...)` exactly — an export that quietly drives maturation off chronological age is a CI failure (Trap 1).
- **Literature validation where supported.** Where a source publishes pediatric concentration-time or maturation-trajectory data (Kataria's own profiles; published neonatal clearance-vs-PMA curves), the kernel output is compared and the agreement recorded — the developmental analog of v0.1's check-against-published-sims, and what lets the maturation-parameter tier be partly numeric.
- **Human verification — the developmental checklist.** The `_checklist_for` builder (verification.py §54) gains a `developmental` group beside the existing ones, confirming the five traps of §4 as explicit yes/no items: (1) maturation driven by **PMA**, not chronological age; (2) allometric exponents fixed-vs-fitted and on the right parameters; (3) reference weight stated; (4) maturation pathway appropriate / labeled if simplified; (5) the size descriptor is a pediatric-valid covariate equation (v0.7). For a *fitted-in-children* record, the human additionally confirms it is curated as evidence (its own tier), not as an annotation. **LLMs assist but never promote** (v0.1 §9), and for maturation parameters the human specifically confirms `TM50`/`Hill` and the age clock against the cited source.
- **`hypnos validate` additions:** the developmental envelope dimensions resolve; every `developmental_model` component carries a citation; `evidence_tier` is D unless backed by a linked fitted-in-children model; `applied_by_default` is never `true` without reviewer sign-off recorded in `extraction`; a `maturation` block's `driver` is a PMA-class clock; the `size_descriptor` resolves in the v0.7 equation library and is valid across the target age range.

---

## 9. Safety & scope guardrails (non-negotiable, inherited and extended)

This is, with [v0.9](../v0.9/pharmacogenomics.md), among the most extrapolation-heavy subsystems in Hypnos, so it carries the strongest guardrails — all of v0.1 §10 / v0.5 §B7 apply, plus pediatric-specific hard lines:

- **No pediatric dose, ever.** The temptation here is the most acute in the project — *"what dose for this neonate?"* Hypnos simulates the forward trajectory of a *given* dose and shows the extrapolation spread; it computes no pediatric dose, no per-kg scaling to a target, no maturation-corrected dose. That is inverse control in its most dangerous form, and it is absolute.
- **Extrapolations off by default, opt-in, always Tier-D, always warned, always cited** (§5). The envelope speaks: a below-derivation-age simulation is greyed + Tier-D + reasoned, so a pediatric extrapolation can never *look* validated.
- **The linear-per-kg shortcut is rendered as visibly wrong forward output, never offered as a number.** Showing the naive rescaling diverging from the mechanistic extrapolation (§6) is the safety message; the view will not present the shortcut, or any extrapolation, as a recommended dose.
- **The `allometry_only` over-dose caveat is always surfaced** when maturation is uncurated, because un-modeled neonatal immaturity over-states clearance exactly where the margin is smallest — a documented failure mode a naive size-only extrapolation would hide.
- `clinicalUse = "PROHIBITED"` is universal and, for this subsystem, doubly load-bearing.

> If a future contributor feels the pull to add "enter the child's weight and age and get the dose," that pull is the exact signal v0.1 §10 named — the subsystem is designed, deliberately, so that feature cannot be added without removing the guardrails that define it.

---

## 10. Worked example (illustrative — numbers pending verification)

```text
$ hypnos compare --drug propofol --pma-weeks 40 --postnatal-days 3 --weight 3.4 --sex M \
                 --developmental --bands --seed 7
term neonate (PMA 40 wk, 3.4 kg)
included (pediatric standing):
  - propofol.kataria_1994     tier C  GREYED: derivation age 3-11 y; neonate is BELOW range -> Tier D
  - propofol.paedfusor        tier C  GREYED: derivation age >1 y;   neonate is BELOW range -> Tier D
  (no model has fitted standing at PMA 40 wk -> everything here is extrapolation)

extrapolation overlay (Tier-D, opt-in):
  - eleveld_2018  + allometry_plus_maturation   Ce peak band [ .. , .. ]  (TM50 cited; maturation curated)
  - schnider_1998 + allometry_only              Ce peak band [ .. , .. ]  CAVEAT: maturation un-modeled
                                                -> neonatal clearance OVER-stated, exposure UNDER-predicted
  - reference line: naive linear-per-kg adult   -> diverges HIGH from both; visibly wrong (ignores maturation)

reading: not one model is validated in this neonate. The mechanistic extrapolations (allometry+maturation)
         spread widely and are all Tier-D; the per-kg shortcut misses badly. The honest output is the SPREAD
         and the labels, never a dose.
```

The reading: v0.1 could only grey the adult and pediatric models here. v0.8 adds the *mechanism* — it shows the allometry+maturation extrapolation, names where maturation is missing (and which direction that errs), and renders the per-kg shortcut as the visibly-wrong line it is. Every curve is Tier-D and labeled; the instrument's value is making the extrapolation honest, not resolving it.

(Every number above is illustrative. Per the project ethos, the size/maturation components are curated `unverified` and await human source confirmation — the PMA-clock and `TM50`/`Hill` line items specifically — before any of this is presented as authoritative.)

---

## 11. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **D0 — Developmental envelope + fitted-pediatric labeling** | proposed | `pma_weeks`/`postnatal_age_days`/`gestational_age_weeks` envelope dimensions; `Envelope.check`/`evaluate_safety` grey + Tier-D below-derivation-age use (promoting `_classify_age_extrapolation`, simulate.py §174, from a flag to a mechanism hook); curate Kataria/Paedfusor as fitted-in-children **evidence** records with their pediatric bands; the developmental verification group. | A neonate correctly greys + Tier-D's the out-of-age models; fitted-pediatric models are labeled as evidence, distinct from extrapolation. |
| **D1 — Allometric size** | proposed | `size_model` schema + `developmental.allometric_scale`; the universal ¾/1.0 scaling as a cited, opt-in, Tier-D extrapolation of adult models; reuse of the v0.7 size-descriptor equation (Al-Sallami across the pediatric range); the visibly-wrong linear-per-kg reference line. | Adult models extrapolate to children by allometry with the per-kg shortcut shown diverging; all Tier-D, all labeled. |
| **D2 — Maturation** | proposed | `maturation_model` schema + `developmental.maturation_factor` (Anderson–Holford PMA sigmoid); the PMA-clock trap enforced; the `allometry_only` over-dose caveat surfaced where maturation is uncurated; curate `TM50`/`Hill` where sources support it. | The neonatal/infant extrapolation composes allometry + maturation; un-modeled maturation is named as an over-dose failure mode. |
| **D3 — Overlay + exports** | proposed | the developmental extrapolation overlay in `compare`/CLI/dashboard; the headline figure (`docs/images/developmental.png`); size + maturation made explicit in NONMEM `$PK`, PharmML, SBML, TCI-JSON, with the PMA-driver comment and round-trip equality. | The overlay tells "fitted vs. extrapolated vs. visibly-wrong-shortcut" for a neonate; the developmental transform round-trips through the formats. |

D0 alone is a useful, self-contained increment — even before any allometric curve is drawn, making the *developmental* envelope speak (and separating fitted-pediatric evidence from extrapolation) turns Hypnos's bare age-range flag into an honest "this is extrapolation, and here is the standing each model actually has."

---

## 12. Cheat sheet (target API)

```python
import numpy as np, hypnos
from hypnos.developmental import allometric_scale, maturation_factor
ds = hypnos.load()

m = ds["hypnotics_iv.propofol.kataria_1994"]
m.applicability_envelope.age_years          # {min: 3, max: 11} — fitted-pediatric EVIDENCE
maturation_factor(pma_weeks=40, tm50=47.7, hill=3.4)   # Anderson-Holford MF for a term neonate

# Neonatal simulation: an adult model extrapolated (opt-in, Tier-D), with maturation
neonate = dict(pma_weeks=40, postnatal_age_days=3, weight=3.4, sex="M")
schedule = [("bolus", 0.0, "2 mg/kg")]
traj = hypnos.simulate(ds, "hypnotics_iv.propofol.eleveld_2018",
                       patient=neonate, schedule=schedule, t=np.linspace(0, 30, 600),
                       developmental=True)        # opt-in; forces Tier-D + warnings
traj.tier                                          # "D"
traj.warnings                                      # ["developmental extrapolation: allometry+maturation; not validated < derivation age", ...]

# The developmental extrapolation overlay — the headline
cmp = hypnos.compare(ds, drug="propofol", purpose="effect_site",
                     patient=neonate, schedule=schedule, developmental=True, seed=7)
cmp.excluded                # adult/pediatric models greyed for age, each reasoned
cmp.extrapolations          # allometry_plus_maturation vs allometry_only bands, all Tier-D, labeled
cmp.linear_per_kg_reference # the deliberately-labeled, visibly-wrong shortcut line
```

```bash
hypnos compare --drug propofol --pma-weeks 40 --postnatal-days 3 --weight 3.4 --sex M --developmental --seed 7
hypnos export  --format nonmem --output exports/nonmem/   # $PK computes allometry + MF(PMA), PMA-driven
hypnos validate                                            # adds PMA-clock + exponent + reference-weight checks
```

---

## 13. Open questions & explicit deferrals

- **Pathway-specific maturation.** Hepatic (CYP, UGT) and renal pathways mature on different PMA trajectories; v0.8 curates a single per-parameter `MF` and *labels* multi-pathway clearances as a simplification. Modeling separate maturation per pathway is a deeper extension deferred until drug-specific sources support it — and it connects to the active-metabolite compartments deferred in v0.5 §F / v0.6 §11.
- **Adolescence, puberty, and the size–maturation handoff.** The transition from the maturation-dominated infant regime to the size-dominated older child is continuous and drug-specific; v0.8 curates the components and lets them compose, but does not assert a universal handoff age — deferred to source-driven curation.
- **Obesity in children.** Pediatric obesity stacks the v0.7 covariate-equation problem (which FFM equation) on the developmental one (which size/maturation model); the two specs share machinery so the stack is *labeled* by worst-tier, but a validated pediatric-obesity model is its own evidence record, not an extrapolation.
- **Maturation of the PD side and of protein binding.** v0.8 is written around PK clearance maturation; PD-parameter maturation (Ce50, γ in children) and the neonatal free-fraction shift are flagged as caveats (the latter shared with v0.5 §S1 / v0.6 LA3) and curated only where published — the never-invent rule again.
- **Pharmacogenetics in the neonate.** Genotype and maturation interact (an immature pathway whose mature capacity is also genetically reduced); v0.8 and [v0.9](../v0.9/pharmacogenomics.md) compose by worst-tier and never-invent, but the validated *interaction* of a maturation curve with a genotype modifier is deferred — each is labeled, their stack is bounded, but their joint quantitative model is not asserted.
