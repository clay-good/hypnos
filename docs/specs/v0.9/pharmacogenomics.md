# Hypnos — design spec (v0.9): the pharmacogenomic-envelope subsystem

**Add the *genetic* special-population axis the way v0.5 added the disease-state axis — as a cited, opt-in, Tier-D `pharmacogenomic_modifier` for the genotypes that actually change anesthetic PK (butyrylcholinesterase variants prolonging succinylcholine/mivacurium and ester-local-anesthetic metabolism; CYP2D6/CYP3A for the opioids and midazolam) — while drawing one new, sharp distinction this subsystem exists to enforce: a pharmacogenetic *susceptibility* is **not** a dose change. Malignant-hyperthermia susceptibility (RYR1/CACNA1S) and atypical-BCHE prolonged block are surfaced as `safetyCritical` *avoidance/awareness flags*, never modeled as a clearance scale-factor — because in anesthesia the honest pharmacogenetic message is overwhelmingly "avoid this trigger / expect a prolonged effect," not "titrate the dose," and saying so is the safety content.**

> **Status: proposed — design spec only; nothing in this document is implemented yet.** Additive over [v0.1](../v0.1/spec.md)–[v0.8](../v0.8/developmental.md): no existing record changes meaning; the pharmacogenomic blocks are new optional structures plus new (declared-only) envelope dimensions, and "not curated" renders as an explicit gap, never as an inferred genotype or a fabricated modifier. Dataset/schema target `0.9.0`. This subsystem is the direct genetic analog of v0.5's `special_populations`: it reuses that spec's evidence-vs-annotation discipline (§B2), its opt-in/Tier-D/always-cited modifier pattern (§B5), and its envelope-speaks design (§B6), and it ties the butyrylcholinesterase story to the ester local anesthetics curated in [v0.6](../v0.6/local_anesthetics.md). It also inherits v0.6's `safetyCritical` flag for the one pharmacogenetic fact that is purely a contraindication, not a kinetic change — malignant hyperthermia. No new mechanism is invented; the machinery is the v0.5 modifier and the v0.6 safety flag, pointed at a genetic covariate.

> This spec inherits the project stance verbatim: *infrastructure, not a simulator; honest about uncertainty by default; a simulation is only as trustworthy as its weakest, least-validated input; no inverse control, ever.* Pharmacogenomics is where the precision-dosing temptation is strongest — *"sequence the patient and compute the personalized dose"* — so v0.9 leans hardest on the never-invent rule (a genotype is never inferred, a modifier never applied without a declared genotype and a citation) and on the inverse-control boundary (a genotype shifts the *forward* prediction of a given dose, and surfaces a contraindication; it never yields a dose). The subsystem's central honesty move, like v0.6's, is to make the *true science* — that most actionable anesthetic pharmacogenetics is avoidance and prolonged-effect awareness, not titration — *be* the safety posture.

---

## 1. The problem — genotype changes the prediction, but mostly it changes the *decision to avoid*

Anesthetic pharmacology has a small set of genetic facts that genuinely matter, and they split cleanly into two kinds that the field — and any naive "pharmacogenomic dosing" tool — routinely blur:

1. **Genotype that changes PK (a kinetic modifier).** A handful of anesthetic and perioperative drugs are cleared by enzymes with clinically important genetic variation:
   - **Butyrylcholinesterase (BCHE, "plasma/pseudo-cholinesterase").** Hydrolyzes **succinylcholine** and **mivacurium**, and the **ester** local anesthetics (procaine, 2-chloroprocaine, tetracaine). Atypical (dibucaine-resistant) variants — heterozygous, and especially homozygous — reduce hydrolysis, so a *given* dose produces a markedly **prolonged** neuromuscular block (succinylcholine apnea lasting hours in the homozygous atypical patient). The dibucaine number is the classic phenotype.
   - **CYP2D6** activates the opioid prodrugs **codeine** and **tramadol** to their active metabolites: poor metabolizers get little analgesia, ultrarapid metabolizers risk toxicity from a standard dose.
   - **CYP3A4/5** metabolizes **fentanyl, alfentanil, sufentanil, midazolam** — genetic variation here is real but generally smaller-effect and less actionable.
2. **Genotype that changes *risk*, not kinetics (a susceptibility).** The single most important anesthetic pharmacogenetic fact is **malignant hyperthermia** susceptibility (variants in **RYR1**, less often **CACNA1S**): exposure to a **volatile agent** or **succinylcholine** can trigger a life-threatening hypermetabolic crisis. This is **not** a clearance change and **not** a dose to adjust — it is an absolute **trigger to avoid**. Treating it as a PK modifier would be a category error with lethal implications.

The crucial, load-bearing observation: **in anesthesia, the genetic facts that matter are mostly of the second kind, or behave like it.** Even the BCHE story is, clinically, less "compute a smaller succinylcholine dose" and more "expect and prepare for prolonged block / consider an alternative relaxant." A pharmacogenomic module that modeled everything as a clearance scale-factor and emitted "personalized doses" would manufacture exactly the false precision — and exactly the forbidden dosing question — that v0.6 was built to refuse for local anesthetics. So, as in v0.6, the deferral-and-careful-framing is the correct path: design the subsystem so the honest science (avoidance and prolonged-effect awareness dominate; kinetic modifiers are few, cited, and wide) *is* the safety message.

Today Hypnos is silent on genotype, and — as v0.5 §B1 said of organ failure — silence reads as "fine." A patient who is homozygous atypical for BCHE, or MH-susceptible, looks to Hypnos exactly like anyone else. v0.9 makes the silence speak in genetic space: it adds genotype as a declared envelope axis so that a relevant genotype either surfaces a cited, opt-in, Tier-D kinetic modifier *or* raises a hard `safetyCritical` avoidance flag — and it keeps those two responses rigorously distinct.

The headline deliverable (§6): the machine-readable, honestly-bounded answer to **"given this declared genotype, which drugs carry a cited kinetic modifier (Tier-D, opt-in, forward-only) and which carry a hard safety/avoidance flag — and which, the honest majority, carry neither because no actionable evidence exists?"**

---

## 2. What we curate — two record kinds, never blurred (the v0.5 discipline, genetic axis)

The cardinal distinction, transposed from v0.5 §B2 and sharpened by §1:

1. **Pharmacogenomic kinetic modifiers — annotations on a PK model.** Where a source documents that a genotype/phenotype changes a drug's disposition (atypical BCHE → reduced succinylcholine/mivacurium hydrolysis; CYP2D6 PM → reduced codeine activation), it is curated as a `pharmacogenomic_modifier` attached to the relevant model: cited, tiered, and applied **only** opt-in, **always** Tier-D unless a model was actually *fitted by genotype*, **always** warned. Like a disease modifier, it is an advisory, cited correction — never a silent reparameterization, never present without a citation.
2. **Pharmacogenomic safety flags — `safetyCritical` susceptibilities, not PK.** Where a genotype confers an adverse-event *susceptibility* that is not a kinetic change (MH from RYR1/CACNA1S; the avoidance dimension of atypical BCHE), it is curated as a `pharmacogenomic_safety_flag` carrying v0.6's `hypnos:safetyCritical` annotation and rendered as a contraindication-style warning. It is **never** expressed as a clearance scale-factor. This is the new discipline v0.9 adds, and the validator enforces the separation (§9).

> The single most important rule in this subsystem: **a kinetic modifier scales a prediction; a safety flag forbids a trigger.** They are different objects with different renderings, and the design refuses to let a susceptibility (MH) be smuggled in as a PK number — the exact failure mode that would turn an avoidance fact into a deceptive dose.

A *fitted-by-genotype* population model (rare; e.g. a PK study stratified by BCHE phenotype) is, as always, a normal Hypnos model record tiered on its own merits — evidence, not annotation.

---

## 3. The pharmacogenetics the envelope and blocks must carry

- **BCHE phenotype** (`normal | heterozygous_atypical | homozygous_atypical`, characterized by dibucaine number): prolonged hydrolysis-dependent effect for **succinylcholine, mivacurium, ester LAs**. Both a *kinetic modifier* (reduced hydrolysis clearance, where quantified) **and** a *safety/awareness flag* (prepare for prolonged block) — curated as both kinds, explicitly, because it is genuinely both.
- **CYP2D6 metabolizer status** (`PM | IM | EM | UM`): activation of **codeine/tramadol** prodrugs; a kinetic modifier on the active-metabolite formation (which connects to the active-metabolite compartments deferred in v0.5 §F — a CYP2D6 modifier is most honestly expressed on a metabolite model when one exists).
- **CYP3A4/5** for **fentanyl/alfentanil/sufentanil/midazolam**: real but generally smaller-effect; curated only where a source quantifies it, and tiered accordingly (often no actionable modifier — an honest gap).
- **RYR1 / CACNA1S → MH susceptibility**: a hard `safetyCritical` avoidance flag on the **volatile agents** and **succinylcholine** — *not* a PK modifier.
- **The substrate-specificity guardrails (the anti-trap facts):**
  - **Remifentanil is hydrolyzed by nonspecific blood/tissue esterases, *not* by BCHE** — so BCHE deficiency does **not** prolong remifentanil. Curated explicitly so the esterase story (already encoded for organ-independence in v0.5 §B3 / Dershwitz/Hoke) is not wrongly generalized to a BCHE modifier.
  - **Amide local anesthetics (lidocaine, bupivacaine, ropivacaine — the v0.6 set) are hepatically metabolized, *not* by BCHE** — only the *ester* LAs are BCHE substrates. The modifier attaches to esters only.
  - **Rocuronium/vecuronium are not cholinesterase substrates** — BCHE does not prolong them (and sugammadex, v0.5 Part A, is their reversal path).

---

## 4. Schema extension & the pharmacogenomic transcription traps

No new `subsystem` value — pharmacogenomic records attach to existing models, as v0.5 modifiers do. The `applicability_envelope` `$def` gains optional, **declared-only** genotype/phenotype dimensions (categorical, reusing the pattern of v0.5's `child_pugh`): `bche_phenotype`, `cyp2d6_phenotype`, `cyp3a_phenotype`, `mh_susceptibility`. New `$defs`: `pharmacogenomic_modifier` (kinetic) and `pharmacogenomic_safety_flag` (susceptibility).

```jsonc
"pharmacogenomic_modifier": {                // KINETIC — scales a prediction
  "gene": "BCHE",
  "phenotype": { "dimension": "bche_phenotype", "value": "homozygous_atypical" },
  "phenotype_basis": "dibucaine_number < 30 (dibucaine-resistant)",   // the genotype->phenotype mapping is lossy (Trap 1)
  "affected_parameter": "Cl",
  "substrate_scope": ["succinylcholine", "mivacurium", "ester_la"],   // NOT remifentanil/amide-LA/rocuronium (§3)
  "adjustment": { "type": "scale_factor", "value": 0.2, "equation": null },  // directionally strong; magnitude wide
  "evidence_tier": "D",                       // modifiers default to D; a fitted-by-genotype model is a DIFFERENT record
  "applied_by_default": false,                // OPT-IN ONLY
  "caveat": "Prolonged block well-documented; the SCALE is illustrative and individual — phenotype, not a precise CL.",
  "primary_citation": "..."
}

"pharmacogenomic_safety_flag": {             // SUSCEPTIBILITY — forbids a trigger; NEVER a PK number
  "gene": "RYR1",
  "phenotype": { "dimension": "mh_susceptibility", "value": true },
  "trigger_agents": ["volatiles", "succinylcholine"],
  "consequence": "malignant_hyperthermia (life-threatening hypermetabolic crisis)",
  "action": "AVOID trigger agents; this is NOT a dose adjustment",
  "safety_critical": true,                    // -> hypnos:safetyCritical (v0.6) in every export
  "evidence_tier": "C",                       // the gene-trigger link; not a kinetic claim
  "primary_citation": "..."
}
```

`evaluate_safety` (simulate.py §211) and `Envelope.check` (models.py §278) extend to the new declared dimensions: a patient who **declares** a relevant genotype/phenotype either (a) triggers a cited, opt-in, Tier-D kinetic modifier (applied only on opt-in, forcing Tier-D + a warning naming the modifier and citation), or (b) raises a `safetyCritical` avoidance flag rendered as a contraindication warning with **no** numeric dose effect. A patient who declares nothing sees no genotype effect (the never-infer rule, §5).

**Pharmacogenomic transcription traps (additions to the trap lists):**

1. **Genotype vs. phenotype vs. metabolizer status (the cardinal mapping sin).** Star-allele diplotypes map to a metabolizer phenotype through a lossy, periodically-revised translation; curating a raw genotype as if it were a phenotype (or vice-versa) is silently wrong. We curate the **phenotype** as the actionable unit and record the genotype/diclplotype basis and the translation source in `phenotype_basis`. The verifier confirms the mapping against the cited source.
2. **Susceptibility mis-modeled as a kinetic scale-factor (the category-error sin).** MH is not a clearance change; modeling it as one would be both wrong and dangerous. The schema *separates* the two `$defs`, and `hypnos validate` rejects a `safety_flag` carrying an `adjustment` and a `modifier` carrying a `trigger_agents` list (§9).
3. **Population allele frequency ≠ individual genotype.** A genotype is never inferred from ancestry or prevalence; modifiers apply only to a **declared** phenotype. (This is also an ethics line — §10.)
4. **Substrate specificity.** BCHE acts on succinylcholine/mivacurium/ester-LAs only — **not** remifentanil, amide-LAs, or rocuronium (§3). A modifier mis-scoped across a class is silently wrong; `substrate_scope` is mandatory and validated against the drug's metabolic route.
5. **Phenotype dose-direction.** A CYP2D6 *prodrug* (codeine) and a CYP2D6 *substrate-to-inactive* behave oppositely under the same phenotype (PM → less active drug for a prodrug, more parent for a directly-active drug); the direction is curated explicitly, never assumed from "PM = less drug."

These become explicit yes/no items in the verification checklist (§9).

---

## 5. Tiers & the never-invent / never-infer rules

- **Kinetic modifiers are Tier-D by construction** (§B5 transposed): they describe a documented *direction and rough magnitude* by phenotype, not a validated by-genotype model. The median line, once a modifier is applied, is labeled Tier-D — you cannot get an A-looking number out of a genotype adjustment. A *fitted-by-genotype* model is a different, separately-tiered record.
- **Safety flags carry their own tier on the gene–trigger/gene–phenotype evidence**, independent of any kinetic claim, and never imply a kinetic tier. An MH flag being Tier-A on the trigger link says nothing about a (non-existent) dose effect.
- **The never-invent rule (non-negotiable).** No modifier without a citation; no borrowing a modifier across drugs or across genes; no synthesized scale-factor on physiological intuition. A missing modifier is an honest gap ("no curated pharmacogenomic adjustment").
- **The never-infer rule (new, genetic-specific, and an ethics line).** Hypnos never infers a genotype — not from ancestry, not from population frequency, not from phenotype-looks-likely. A modifier or flag fires **only** on a genotype/phenotype the caller explicitly declares. Absent a declaration, there is no genetic effect, and that absence is shown as "genotype not declared," never as "wild-type assumed."

---

## 6. Headline feature — the pharmacogenomic safety overlay (avoidance first, kinetics second)

The divergence/safety view (`compare`/`evaluate_safety`, simulate.py §655/§211) gains a pharmacogenomic overlay that makes the genetic facts visible *in the right register* — and the ordering is the message:

- **Avoidance flags first, and loudest.** Declare MH susceptibility → the volatile agents and succinylcholine are rendered with a hard `safetyCritical` **AVOID** flag and the consequence named; no curve is dose-adjusted, because there is nothing to adjust — the honest output is a contraindication, not a number. Declare homozygous-atypical BCHE → succinylcholine/mivacurium/ester-LAs carry the **prolonged-block awareness** flag.
- **Kinetic modifiers second, as opt-in Tier-D extrapolation bands.** For the drugs that *do* carry a cited kinetic modifier, the genotype fans the forward prediction into an explicitly Tier-D band (the prolonged-effect trajectory for a *given* dose), labeled as a phenotype extrapolation, with the substrate scope and the wide-magnitude caveat shown. For atypical BCHE + succinylcholine the band's punchline is the recovery time stretching out — rendered as honest forward output, exactly as v0.5's naloxone renarcotization relapse is.
- **The honest majority: neither.** For most drugs and most genotypes, the overlay says *"no actionable pharmacogenomic evidence curated"* — and that, the view makes clear, is the truthful state of anesthetic pharmacogenetics for the bulk of the formulary, not a gap to apologize for.

The honest message the view delivers: *"the genetic facts that matter here are mostly about what to avoid and what prolonged effect to expect, not a dose to compute — and for this genotype, here are the few cited kinetic modifiers (all Tier-D, all opt-in) and the safety flags (rendered as avoidance, never as numbers)."*

---

## 7. Export

Both record kinds flow through the existing exporters with the v0.5/v0.6 handling:

| Surface | v0.9 behavior |
| --- | --- |
| **TCI-sim JSON / dataset** | `pharmacogenomic_modifier` / `pharmacogenomic_safety_flag` blocks pass through verbatim (lossless), carrying the opt-in/Tier-D flags and the `safety_critical` marker, so a downstream consumer cannot apply a modifier without seeing it is unvalidated, or read a safety flag as a kinetic effect. |
| **PharmML / SO · NONMEM · nlmixr2 / Pumas · SBML** | kinetic modifiers attach as `hypnos:pharmacogenomicModifier` MIRIAM/RDF on the model (genotype is a covariate the core formats describe only weakly); safety flags attach as `hypnos:safetyCritical` + `hypnos:pharmacogenomicSafetyFlag`, **explicitly not** as any parameter change — the RDF carrier keeps the avoidance fact lossless and tool-portable without ever encoding it as a number a deterministic consumer might apply. |

Every artifact inherits the universal `clinicalUse = "PROHIBITED"` annotation; pharmacogenomic records additionally carry `hypnos:safetyCritical` where a safety flag applies, doubly warning a downstream consumer — the same belt-and-braces v0.6 gave LA records.

---

## 8. Validation & verification

- **Round-trip (CI), unchanged discipline.** Every kinetic-modifier application is re-simulated through the export and checked against the reference within tolerance (v0.1 §6). A *safety flag* has nothing to round-trip kinetically — and the test asserts exactly that: a safety flag must produce **no** change in the simulated trajectory, only a warning (the machine-checkable form of Trap 2).
- **Literature validation where supported.** Where a source publishes recovery-time data by BCHE phenotype, the kinetic-modifier band is compared to it and the agreement recorded — the pharmacogenomic analog of v0.5's sugammadex-recovery check.
- **Human verification — the pharmacogenomic checklist.** The `_checklist_for` builder (verification.py §54) gains a `pharmacogenomics` group confirming the five traps of §4 as explicit yes/no items: (1) genotype→phenotype mapping correct and sourced; (2) susceptibility curated as a **flag**, not a modifier; (3) phenotype **declared**, never inferred; (4) substrate scope correct (BCHE excludes remifentanil/amide-LA/rocuronium); (5) phenotype dose-direction correct (prodrug vs. directly-active). **LLMs assist but never promote** (v0.1 §9); for a kinetic modifier the human confirms the direction and magnitude against the cited source, and for a safety flag confirms it carries no numeric effect.
- **`hypnos validate` additions:** every modifier/flag `primary_citation` resolves; a `pharmacogenomic_safety_flag` may **not** carry an `adjustment` and a `pharmacogenomic_modifier` may **not** carry `trigger_agents` (Trap 2, enforced by construction); `evidence_tier` is D for a kinetic modifier unless backed by a linked fitted-by-genotype model; `applied_by_default` is never `true`; `substrate_scope` is consistent with the drug's curated metabolic route (so a BCHE modifier cannot attach to remifentanil or an amide LA).

---

## 9. Safety & scope guardrails (non-negotiable, inherited and extended)

With [v0.8](../v0.8/developmental.md), this is among the most guardrail-heavy subsystems in Hypnos — all of v0.1 §10 / v0.5 §B7 / v0.6 §7 apply, plus genetic-specific hard lines:

- **No genotype-guided dose, ever.** "Sequence the patient, compute the personalized dose" is the precision-medicine framing of the exact inverse-control function the project forbids. A genotype shifts the *forward* prediction of a *given* dose and may raise an avoidance flag; it never yields a dose, a per-genotype titration, or a target. Absolute.
- **A susceptibility is an avoidance flag, never a number.** MH (and the avoidance dimension of atypical BCHE) is rendered as a contraindication-style warning with no kinetic effect (§2, Trap 2). Encoding a susceptibility as a clearance change would convert a "don't use this trigger" into a deceptive "use less" — the category error the schema separation exists to make impossible.
- **Never infer a genotype** (§5). Modifiers and flags fire only on a **declared** genotype/phenotype; nothing is assumed from ancestry, frequency, or likelihood. This is an ethics line as much as a scientific one — Hypnos does not profile.
- **Modifiers off by default, opt-in, always Tier-D, always cited, always warned**, and **substrate-scoped** so a genetic effect never leaks to a non-substrate drug (no BCHE effect on remifentanil/amide-LA/rocuronium).
- **The honest-majority message is preserved**: for most of the formulary, the truthful pharmacogenomic answer is "no actionable evidence," and the view says so rather than manufacturing a modifier to look thorough.
- `clinicalUse = "PROHIBITED"` is universal; `hypnos:safetyCritical` rides on every safety-flag record.

> If a future contributor feels the pull to add "enter the genotype and get the personalized dose," that pull is the exact signal v0.1 §10 named — the subsystem is designed, deliberately, so that feature cannot be added without removing the schema separation and the never-infer rule that define it.

---

## 10. Worked example (illustrative — numbers pending verification)

```text
$ hypnos pgx --bche homozygous_atypical --mh-susceptible true
declared genotype overlay (research/education; NOT a dosing tool):

  SAFETY / AVOIDANCE FLAGS (not dose adjustments):
    RYR1 (MH-susceptible)  -> AVOID volatiles + succinylcholine   [safetyCritical, Tier C link]
                              consequence: malignant hyperthermia. This is a contraindication, NOT a dose.
    BCHE homozygous atypical -> succinylcholine / mivacurium / ester-LAs: EXPECT PROLONGED BLOCK [awareness]

  KINETIC MODIFIERS (opt-in, Tier D, forward-only):
    succinylcholine + BCHE(hom. atypical)  -> recovery time band stretches markedly (illustrative scale; wide)
    mivacurium      + BCHE(hom. atypical)  -> prolonged hydrolysis-dependent recovery (Tier D, cited)
    NOT applied: remifentanil (nonspecific esterase, not BCHE) · lidocaine/bupivacaine (amide, hepatic)
                 · rocuronium (not a cholinesterase substrate)  <- substrate-scope guardrail

  NO ACTIONABLE PGx EVIDENCE CURATED: propofol, dexmedetomidine, sevoflurane-PK, ...
    (the honest majority — most anesthetic PGx is avoidance/awareness, not titration)
```

The reading: the overlay leads with what to *avoid* (MH, the highest-stakes fact, rendered as a flag with no number), then the prolonged-effect *awareness* for BCHE substrates, then the few opt-in Tier-D kinetic bands for a *given* dose — and it is explicit that remifentanil, the amide LAs, and rocuronium are *not* affected (the substrate guardrail), and that most of the formulary has no actionable evidence. Nothing is a dose; the susceptibility is never a number.

(Every scale above is illustrative. Per the project ethos, the modifiers and flags are curated `unverified` and await human source confirmation — the genotype→phenotype mapping and the substrate scope specifically — before any of this is presented as authoritative.)

---

## 11. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **G0 — Safety flags first (MH + BCHE awareness)** | proposed | the declared genotype envelope dimensions; the `pharmacogenomic_safety_flag` `$def` + `hypnos:safetyCritical` rendering; RYR1/CACNA1S → volatile/succinylcholine **AVOID** flag; atypical-BCHE prolonged-block **awareness** flag; the `validate` rule that a flag carries no kinetic effect; the pharmacogenomics verification group. | The highest-stakes, non-kinetic genetic facts are surfaced as avoidance/awareness flags — and the schema makes it impossible to encode them as numbers. Safest possible entry point. |
| **G1 — BCHE kinetic modifier** | proposed | the `pharmacogenomic_modifier` `$def`; the atypical-BCHE reduced-hydrolysis modifier on succinylcholine/mivacurium/ester-LAs as opt-in, Tier-D, substrate-scoped; the prolonged-recovery forward band; the substrate-scope guardrail (excludes remifentanil/amide-LA/rocuronium). | The canonical anesthetic pharmacogenetic *kinetic* fact simulates as a forward, Tier-D, opt-in prolonged-effect band — no dose computed, scope enforced. |
| **G2 — Opioid CYP modifiers** | proposed | CYP2D6 (codeine/tramadol activation; the prodrug dose-direction trap) and CYP3A (fentanyl family/midazolam, where quantified) modifiers; the link to active-metabolite compartments (deferred v0.5 §F) where a metabolite model exists; honest "no actionable modifier" where evidence is thin. | The opioid/benzodiazepine CYP story is curated where actionable and honestly absent where not. |
| **G3 — Overlay + exports** | proposed | the pharmacogenomic safety overlay in `compare`/CLI (`hypnos pgx`)/dashboard, avoidance-first; the headline figure; modifiers/flags exported as `hypnos:` RDF (never as parameter changes for flags); round-trip null-effect test for flags. | The overlay tells "avoid / expect-prolonged / opt-in-Tier-D-band / no-evidence" for a declared genotype; the records round-trip with the flag/modifier distinction intact. |

G0 alone is a useful, self-contained increment **and** the safest entry point: it surfaces the avoidance facts (MH especially) that matter most, *before* any kinetic modifier draws a single band — so the subsystem's first release teaches the most important pharmacogenetic idea (this is about avoidance, not titration) without yet computing anything that could be misread as a personalized dose.

---

## 12. Cheat sheet (target API)

```python
import hypnos
ds = hypnos.load()

# Declared genotype overlay — avoidance flags + opt-in Tier-D kinetic modifiers
genotype = dict(bche_phenotype="homozygous_atypical", mh_susceptibility=True)
pgx = hypnos.pgx_overlay(ds, genotype=genotype)
pgx.safety_flags        # [{gene: "RYR1", action: "AVOID volatiles + succinylcholine", safety_critical: True}, ...]
pgx.modifiers           # opt-in Tier-D kinetic modifiers, substrate-scoped (succinylcholine/mivacurium/ester-LA)
pgx.no_evidence         # the honest majority of the formulary

# Forward simulation WITH a declared genotype (opt-in modifier; forward-only; no dose computed)
patient = dict(age=40, weight=70, height=175, sex="M", bche_phenotype="homozygous_atypical")
traj = hypnos.simulate(ds, "nmb_agents.succinylcholine.<model>",
                       patient=patient, schedule=[("bolus", 0.0, "1 mg/kg")],
                       pharmacogenomics=True)     # opt-in; forces Tier-D + prolonged-block warning
traj.tier               # "D"
traj.warnings           # ["BCHE homozygous atypical: prolonged block; scale illustrative, individual", ...]

# A safety flag changes NO number — it forbids a trigger
mh = hypnos.pgx_overlay(ds, genotype=dict(mh_susceptibility=True))
mh.safety_flags[0].affects_kinetics   # False  (by construction — Trap 2)
```

```bash
hypnos pgx --bche homozygous_atypical --mh-susceptible true        # avoidance-first overlay
hypnos simulate --model succinylcholine.<model> --bche homozygous_atypical --pharmacogenomics
hypnos validate                                                     # flag-has-no-adjustment + substrate-scope checks
```

---

## 13. Open questions & explicit deferrals

- **Genotype as a continuum / activity scores.** Metabolizer phenotype is a coarse bin over a continuous activity-score landscape (CYP2D6 especially); v0.9 curates the phenotype bin and defers continuous activity-score modeling until sources support it without implying false precision — the same stance v0.5 took on Child-Pugh-as-categorical.
- **Active-metabolite compartments.** CYP2D6 modifiers act most honestly on *metabolite formation* (codeine→morphine), which needs the explicit active-metabolite compartments deferred in v0.5 §F and v0.6 §11; until those land, a CYP2D6 modifier is a caveat on the parent model, labeled as a simplification.
- **Gene–gene and gene–environment interactions.** A patient may carry multiple relevant variants, or a genotype plus an enzyme-inhibiting co-medication; the worst-tier/never-invent rules cover labeling, but validated *interaction* of modifiers (or with the v0.8 maturation curve, §13 there) is deferred — each is labeled, the stack is bounded by worst-tier, the joint quantitative model is not asserted.
- **Pharmacodynamic pharmacogenetics.** v0.9 is written around PK (metabolism) genotypes; PD-side genetic variation (opioid-receptor variants and analgesic response) is real but far less settled, and is deferred to source-driven curation rather than asserted.
- **Ethics and data handling.** Genotype is sensitive personal data. v0.9's never-infer rule (§5) is partly an ethics guardrail: Hypnos curates *what a declared genotype implies for a model*, never a person's genotype, never an inference from ancestry — consistent with v0.4 §3's "derived facts only, never raw subject data" discipline, here applied to the most sensitive covariate of all.
