# Hypnos — design spec (v0.6): the local-anesthetic systemic-toxicity (LAST) subsystem

**Bring local anesthetics (lidocaine, bupivacaine, levobupivacaine, ropivacaine, mepivacaine) into Hypnos with the safety framing v0.1 explicitly deferred — by curating (1) *site-specific systemic-absorption* PK (the rate of rise of plasma concentration depends far more on *where* the block is placed than on the milligrams), (2) documented systemic-toxicity *threshold ranges* (CNS then cardiovascular) as exactly that, *ranges* with large, honest uncertainty, never a single "toxic at X," (3) the protein-binding / free-fraction story including its saturation nonlinearity, and (4) the stereochemistry that makes bupivacaine more cardiotoxic than ropivacaine/levobupivacaine — and by reframing the obvious-but-forbidden "margin to toxicity" view as a *double-uncertainty* honesty instrument that discourages over-trust rather than a dosing calculator that invites it.**

> **Status:** **LA0 + LA1 + LA2 shipped (disposition + site absorption + binding + the toxicity-threshold *ranges* and the double-uncertainty view + the stereochemistry/cardiotoxicity differentiation); LA3 proposed.** The `local_anesthetics` subsystem is live: an additive `absorption` schema block (site rates with a robust `rank` + Tier-C `ka`), three curated LA models (lidocaine/bupivacaine/ropivacaine, Tucker & Mather 1979, `unverified`), the forward-only Bateman kernel + `site_comparison` in `hypnos.la`, the `hypnos la` CLI surfacing the **site-of-injection dominance**, and drug-level binding curated and correctly kept off the v0.5 albumin axis (LA binding is α1-AG-driven). **LA1** adds the `toxicity_threshold` block — always a **range**, never a line (the schema forbids a single-value threshold by construction; `hypnos validate` enforces `low < high`, a mandatory `basis`, and a free-fraction `saturation_caveat` whenever a saturable drug carries a total-plasma threshold), the `la.free_concentration` transform (linear + the saturation under-prediction caveat surfaced, never hidden), the headline **double-uncertainty view** (`la.double_uncertainty` / `hypnos la --site …`) whose honest punchline is that the threshold uncertainty *dwarfs* the PK spread so no single safe-concentration line is defensible, the LA verification checklist (§8), and the `hypnos:safetyCritical` flag + threshold-range RDF carried through every export (thresholds export **as ranges**; TCI-JSON verbatim). Curated CNS/cardiovascular ranges for lidocaine (Tucker & Mather 1979), bupivacaine and ropivacaine (Knudsen 1997, Scott 1989), all `unverified` Tier-C. **LA2** adds the stereochemistry/cardiotoxicity differentiation: a drug-level `cardiotoxicity_class` (rank + stereochemistry + the qualitative CNS-to-CVS margin, curated like `protein_binding` — see the design note below), a new **levobupivacaine** drug + model (Bardsley 1998) anchoring the intermediate point, the `la.cardiotoxicity_comparison` agent-choice axis (`hypnos la --agents`) whose numeric CNS-to-CVS fold-margin is monotone with the qualitative ranking (racemic bupivacaine 1.4× < levobupivacaine 1.7× < ropivacaine 2.0× < lidocaine 5.0×), the class surfaced in the double-uncertainty view + a narrow-margin cardiotoxicity warning, and `hypnos:cardiotoxicityClass` carried through every export. (Design note: `cardiotoxicity_class` is curated at the **drug** level beside `protein_binding` — both are intrinsic chemistry, shared across a drug's models — rather than as the model-schema `$def` §4 sketched; this matches the established `protein_binding` precedent and is validated for citation/vocabulary in `hypnos validate`.) Still proposed: the nonlinear binding-saturation free-fraction model (LA3). Additive over [v0.1](../v0.1/spec.md)–[v0.5](../v0.5/breadth.md). Dataset/schema target `0.6.0`. This subsystem is the one v0.1 twice held back on purpose: v0.1 §2 lists `local_anesthetics` ⚠️ as **deferred** ("the safety framing needs its own pass before shipping"), and v0.1 §10 repeats *"local-anesthetic systemic-toxicity thresholds are deferred until the safety framing for that subsystem is designed separately."* This document **is** that separate safety-framing pass. It is written safety-first: the guardrails (§7) are not an appendix but the reason the subsystem exists in this shape.

> This spec inherits the project stance verbatim and leans on it harder than any prior spec: *infrastructure, not a simulator; honest about uncertainty by default; no inverse control, ever.* LAST is a leading cause of regional-anesthesia mortality, and a plasma-concentration-vs-threshold view is *precisely* the artifact that tempts a reader to ask "so what's the max safe dose?" — the one question Hypnos must never answer. The design's central move is to make the honest scientific content (absorption is site-driven; thresholds are wide and individual; total concentration hides free-drug risk) *be* the safety message.

---

## 1. The problem — and why it needed its own spec

Every other Hypnos subsystem answers "what concentration/effect results from this dose?" For local anesthetics that forward question has a dangerous shadow: **"what is the maximum safe dose?"** The classic mg/kg "maximum recommended doses" are themselves scientifically shaky — they are not derived from a clean concentration-toxicity model, they ignore that **site of injection dominates systemic absorption**, and they collapse a wide, individual, method-dependent toxicity distribution into a single deceptively precise number. A naive Hypnos LA module that plotted predicted plasma concentration against a single toxicity line would manufacture exactly that false precision and invite exactly that forbidden question.

So the deferral in v0.1 was correct, and the resolution is not "add LA like any other drug" but "design the framing such that the *honest science* and the *safety posture* are the same thing." Three facts about LA pharmacology make that possible, because each one, stated truthfully, *undercuts* the max-dose framing:

1. **Absorption is site-driven, not milligram-driven.** The same dose produces wildly different peak plasma concentrations depending on the vascularity of the injection site (the well-documented rank order: intercostal > caudal > epidural > brachial plexus > sciatic/femoral > subcutaneous tissue infiltration). A single mg/kg ceiling is therefore unreliable on its face — and saying so *is* the safety message.
2. **Toxicity thresholds are wide, individual, and method-dependent.** The plasma concentration at which CNS symptoms (perioral numbness → tinnitus → seizures) and then cardiovascular toxicity (conduction block → arrhythmia → arrest) appear is a *range*, varies between people, depends on speed of rise, acid-base status, and the assay/basis used to report it. Curating it as a range with explicit uncertainty *prevents* the single-line reading.
3. **Toxicity tracks *free* drug, and binding can saturate.** LAs are highly protein-bound (bupivacaine especially); only the *free* fraction is pharmacologically active and toxic. At high total concentrations binding **saturates**, so the free fraction rises *nonlinearly* — total concentration under-predicts risk exactly when risk is highest. This is a genuine, citable failure mode, and it is invisible to any "total concentration vs. line" view.

The deliverable (§6) is therefore not a margin calculator. It is a *double-uncertainty* view: the predicted plasma-concentration band (with all its site/PK uncertainty) shown against the *threshold band* (with all its individual/method uncertainty), whose honest punchline is usually *"these uncertainties are so wide that no single safe-dose number is defensible"* — the v0.2/v0.3 variance-decomposition ethos applied to a safety-critical question.

---

## 2. Scope — what this subsystem does and does not contain

**In scope:** systemic-absorption PK by injection site; disposition PK; protein binding and its saturation; documented systemic-toxicity threshold *ranges* (CNS and cardiovascular) with their uncertainty and basis; the stereochemistry that differentiates cardiotoxicity (racemic bupivacaine vs. levobupivacaine vs. ropivacaine).

**Explicitly out of scope (declared, like every Hypnos exclusion):**

- **Any maximum-safe-dose output, margin-to-toxicity number presented as a safety guarantee, or "is this dose safe?" answer.** Hard line (§7). This is the inverse-control boundary in its sharpest form.
- **Local (compartment) anesthetic *efficacy* / block onset-and-duration at the nerve.** This subsystem is about *systemic* concentration and *systemic* toxicity, not whether the block works. (Block efficacy is a different, deferred modeling problem.)
- **LAST treatment / resuscitation.** Lipid emulsion ("intralipid") rescue, ACLS modification, and airway management are clinical therapeutics, not modeled here. They are named in §7 only to mark the boundary.
- **Regional-technique selection or guidance** of any kind.

---

## 3. What we curate — the LA model record

A LA record extends the standard Hypnos model object with LA-specific blocks, each tiered and cited like everything else.

### 3.1 Site-specific systemic absorption

The key LA-specific structure: a first-order (or biexponential, where published) systemic absorption with a **site-dependent** rate, so the *rate of rise* of plasma concentration is a curated function of injection site, not an afterthought:

```jsonc
"absorption": {
  "model": "first_order",                  // first_order | biexponential
  "site_rates": [
    { "site": "intercostal",        "ka": 0.20, "bioavailability": 1.0, "rank": 1, "citation": "..." },
    { "site": "caudal_epidural",    "ka": 0.12, "bioavailability": 1.0, "rank": 2, "citation": "..." },
    { "site": "lumbar_epidural",    "ka": 0.08, "bioavailability": 1.0, "rank": 3, "citation": "..." },
    { "site": "brachial_plexus",    "ka": 0.05, "bioavailability": 1.0, "rank": 4, "citation": "..." },
    { "site": "subcutaneous",       "ka": 0.02, "bioavailability": 1.0, "rank": 5, "citation": "..." }
  ],
  "tier": "C",
  "extraction": { "review_status": "unverified", "source_locator": "..." }
}
```

The `rank` field encodes the *direction* (the absorption rank order is more robustly established than any single `ka`), so even where the absolute `ka` is Tier-C, the qualitative site-dependence is curated and surfaced.

### 3.2 Protein binding (and its saturation)

```jsonc
"protein_binding": {
  "fraction_bound": 0.95,                  // typical total binding
  "saturable": true,                       // bupivacaine: binding saturates -> free fraction rises nonlinearly
  "free_fraction_model": "...",            // where published; else flagged as a documented failure mode only
  "binding_protein": "AAG_and_albumin",
  "tier": "C", "extraction": { "review_status": "unverified" }, "primary_citation": "..."
}
```

### 3.3 Toxicity thresholds — *ranges*, never lines

The safety-critical block. A new `$def` `toxicity_threshold`, always a **range with uncertainty and an explicit basis**:

```jsonc
"toxicity_thresholds": [
  {
    "endpoint": "cns_first_symptoms",      // cns_first_symptoms | cns_seizure | cardiovascular
    "concentration_range": { "low": 2.0, "high": 4.0, "units": "ug/mL" },
    "basis": "total_plasma",               // total_plasma | free_plasma  (drives the §3.2 free-fraction caveat)
    "individual_variability": "high",      // qualitative; thresholds are person- and rate-of-rise-dependent
    "method_caveat": "Volunteer infusion studies; speed of rise and acid-base status shift thresholds.",
    "tier": "C", "primary_citation": "..."
  }
]
```

Thresholds are **always** ranges; the schema forbids a single-value threshold (§4). CNS-before-cardiovascular ordering (the clinical warning sequence) is curated as the endpoint progression, with the explicit caveat that for the most cardiotoxic agents the CNS-warning margin can be narrow or absent (bupivacaine).

### 3.4 Stereochemistry / cardiotoxicity differentiation

Racemic **bupivacaine** is markedly more cardiotoxic than its S-enantiomer **levobupivacaine** or **ropivacaine** ("fast-in, slow-out" avid cardiac sodium-channel binding). Curated as a per-drug `cardiotoxicity_class` with citation, so the divergence view can show *why* agent choice matters for the cardiovascular threshold even at similar CNS thresholds.

---

## 4. Schema extension & the canonical-representation traps

New `subsystem: "local_anesthetics"`. New `$defs`: `absorption` (with `site_rates`), `protein_binding`, `toxicity_threshold`, `cardiotoxicity_class`. The `applicability_envelope` reuse from v0.5 (organ-function dimensions) applies — hepatic impairment and acidosis shift LA disposition and binding. `required` lists for existing records unchanged.

**Hard schema constraints, because this is safety-critical:**

- A `toxicity_threshold` **must** be a `concentration_range` with `low < high`; a single-value threshold is a schema error. (No false-precision lines, by construction.)
- Every threshold **must** declare its `basis` (`total_plasma` vs `free_plasma`); a `total_plasma` basis on a `saturable: true` drug **must** carry the free-fraction failure-mode caveat. `hypnos validate` enforces both.
- `absorption.site_rates` must include the `rank` even where `ka` is null (direction curated even when magnitude is not).

**Transcription traps (the LA-specific additions to the v0.2/v0.3 trap lists):**

1. **Total vs free concentration.** A threshold reported on *free* drug compared against a *total* prediction (or vice-versa) is silently wrong — the binding fraction is the conversion. The `basis` field is mandatory precisely to prevent this.
2. **Salt vs base, and concentration units.** LA doses are given as the salt (mg of hydrochloride); plasma assays may report base. Units and salt form are curated explicitly.
3. **Speed-of-rise dependence.** A threshold from a slow volunteer infusion does not transfer to a rapid intravascular injection; the `method_caveat` records the derivation context so the number is never read as context-free.

---

## 5. Reference kernels & tiers

- **`la.absorption_pk(dose, site, disposition_micro, t)`** — site-selected first-order absorption feeding the standard linear disposition PK (the matrix-exponential solver, reference.py §147), forward-only. Output: total plasma concentration trajectory.
- **`la.free_concentration(c_total, protein_binding)`** — the bound→free transform, *with the saturation nonlinearity* where curated; where only the linear binding fraction is known, it returns the linear free concentration **and** attaches the saturation failure-mode caveat so the under-prediction-at-high-concentration risk is never hidden.
- Tiers on the usual A–D spine. Toxicity thresholds will rarely exceed Tier-C (they are wide, old, method-dependent) — and the design *embraces* that: a Tier-C threshold range shown honestly is the point. The never-invent rule (v0.2) is absolute here: **no synthesized threshold, no borrowed binding model, no imputed site rate.** A missing number is a stated gap.

---

## 6. The headline feature — the double-uncertainty view (honesty, not advice)

The obvious LA visualization is "plasma concentration vs. toxicity threshold." Done naively it is a dosing calculator. Done in Hypnos's idiom it is the opposite — a **double-uncertainty** instrument:

- Overlay the **predicted total-plasma-concentration band** (PK + site + v0.2/v0.3 uncertainty, where curated) against the **toxicity-threshold band** (the curated range, with its individual/method uncertainty), for a *given* forward dose and site.
- Add the **free-concentration trace** beside the total, so the binding-saturation gap is visible whenever it matters.
- Apply the v0.2/v0.3 variance-decomposition ethos: report *which uncertainty dominates*. For LA the honest finding is usually that the **threshold uncertainty dwarfs the PK uncertainty** — i.e., the answer to "how safe is this?" is *"the threshold itself is too uncertain to draw a line."* That conclusion is the safety message, generated by the science.
- The view is framed throughout as a **research/education** question — *"how do site, agent, and binding move the predicted concentration relative to the range of published toxicity thresholds, and how uncertain is all of it?"* — **never** *"is this dose safe?"* It will not report a margin number as a guarantee, will not rank doses, and will not compute a maximum.

The same `compare`/band machinery (simulate.py §527/§589) drives it, so it never drifts from the CLI — and the agent-choice axis lets it show, e.g., bupivacaine's narrower cardiovascular margin vs. ropivacaine for the same block, the genuinely useful educational comparison.

---

## 7. Safety & scope guardrails (the reason this subsystem has its own spec)

These are not an appendix; they are the design. All of v0.1 §10 / v0.2 §10 / v0.5 §A6,§B7 apply, plus LA-specific hard lines:

- **No maximum-safe-dose. No margin-as-guarantee. No "is this safe?" output. Ever.** Hypnos simulates the forward concentration trajectory of a *given* dose at a *given* site and shows it against a *range* of thresholds. It computes no dose, no ceiling, no probability-of-toxicity target. This is the inverse-control boundary (v0.1 §10) in its most safety-critical form, and it is absolute.
- **Thresholds are ranges, individual, and method-dependent — never a single line.** Enforced in the schema (§4). The view's design exists to *prevent* the single-line reading, not enable it.
- **The free-fraction saturation nonlinearity is always surfaced** for binding-sensitive agents, because total concentration under-predicts risk exactly when risk is highest — a documented failure mode that a naive view would hide.
- **Site-of-injection dominance is the headline, not a footnote** — because "the mg/kg max-dose tables are unreliable since absorption is site-driven" is itself the safety message this subsystem most needs to deliver.
- **LAST treatment (lipid rescue, ACLS) is out of scope** and named as such, so no one mistakes a concentration model for a management tool.
- `clinicalUse = "PROHIBITED"` is universal, and LA records additionally carry a `hypnos:safetyCritical = true` annotation so a downstream consumer is doubly warned. The README/CLI/dashboard boundary language is reinforced specifically for this subsystem.

> If a future contributor feels the pull to add "just enter the block and dose and tell me if it's safe," that pull is the exact signal v0.1 §10 named — the subsystem has been designed, deliberately, to make that feature impossible to add without removing the guardrails that define it.

---

## 8. Validation & verification

- **Round-trip (CI):** absorption + disposition + free-fraction kernels re-simulated through every export and checked against the reference within tolerance (v0.1 §6).
- **Literature validation where supported:** published plasma-concentration-time profiles after defined blocks (the site-absorption studies) are compared to the kernel and the agreement recorded; this is what lets the site-`ka` tier be partly numeric.
- **Human verification — the LA checklist** (a new group in `_checklist_for`, verification.py §54): (1) threshold basis total vs free, and conversion correct (Trap 1); (2) salt vs base and units (Trap 2); (3) speed-of-rise/method context recorded (Trap 3); (4) saturation caveat present for binding-sensitive agents; (5) threshold curated as a *range*, not a line. **LLMs assist but never promote** (v0.1 §9), and given the safety stakes, threshold ranges specifically require human source confirmation before `verified`.
- **`hypnos validate` additions:** thresholds are ranges with `low < high`; `basis` present; saturable+total → caveat present; site rates carry `rank`; all citations resolve.

## 9. Export

LA records flow through the existing exporters as standard PK models plus `hypnos:` RDF for the LA-specific blocks (`hypnos:siteAbsorption`, `hypnos:proteinBinding`, `hypnos:toxicityThresholdRange`, `hypnos:cardiotoxicityClass`), carrying the `hypnos:safetyCritical` flag and the universal `clinicalUse = "PROHIBITED"`. TCI-JSON passes the blocks verbatim (lossless). Threshold ranges export *as ranges* — there is no format projection that collapses them to a single value.

## 10. Phased roadmap

| Phase | Status | Content | Done = |
| --- | --- | --- | --- |
| **LA0 — Disposition + site absorption + binding** | 🟢 shipped | `local_anesthetics` subsystem + additive `absorption` schema; one-compartment disposition + site-specific absorption (rank robust, `ka` Tier-C) + protein binding for lidocaine, bupivacaine, ropivacaine (Tucker & Mather 1979, `unverified`); the `hypnos.la` Bateman kernel + `site_comparison`; the `hypnos la` CLI; **no thresholds**. | The systemic concentration trajectory by site is simulable (`hypnos la`); the site-rank dominance is curated and surfaced as the headline. |
| **LA1 — Toxicity thresholds + the double-uncertainty view** | 🟢 shipped | `toxicity_threshold` ranges (CNS first-symptoms / seizure / cardiovascular) with mandatory basis + uncertainty; the double-uncertainty view (§6, `la.double_uncertainty` + `hypnos la --site`); the free-concentration trace (`la.free_concentration`, linear + saturation caveat); schema range-not-line + saturable→caveat enforcement in `hypnos validate`; the LA verification checklist (§8); `hypnos:safetyCritical` + threshold-range RDF in every export (TCI-JSON verbatim). Lidocaine (Tucker & Mather 1979), bupivacaine + ropivacaine (Knudsen 1997 / Scott 1989), `unverified` Tier-C. | The honest concentration-vs-threshold-*band* view exists, framed as research/education, with the threshold uncertainty made explicit — and the dominant uncertainty named (it is the threshold's). |
| **LA2 — Stereochemistry / cardiotoxicity** | 🟢 shipped | Drug-level `cardiotoxicity_class` (rank + stereochemistry + CNS-to-CVS margin) for lidocaine/bupivacaine/ropivacaine + a new **levobupivacaine** drug & model (Bardsley 1998); the `la.cardiotoxicity_comparison` agent-choice axis (`hypnos la --agents`) with a numeric CNS-to-CVS fold-margin monotone with the ranking; the class surfaced in the double-uncertainty view + a narrow-margin warning; `hypnos:cardiotoxicityClass` in every export; a `la_cardiotoxicity.png` figure. | The view shows *why* agent choice changes the cardiovascular margin (bupivacaine 1.4× < levobupivacaine 1.7× < ropivacaine 2.0× < lidocaine 5.0×), the key educational comparison. |
| **LA3 — Binding saturation failure mode** | proposed | the nonlinear free-fraction model where published; the documented saturation failure mode surfaced wherever a total-basis threshold meets a saturable drug. | Total concentration's under-prediction of free-drug risk at high concentration is made visible, not hidden. |

LA0 alone is a useful, self-contained increment **and** the safest possible entry point: it adds the systemic-absorption science (the site-dominance message) *before* any toxicity threshold, so the subsystem's first release teaches the most important safety idea — that the milligram ceiling is the wrong mental model — without yet drawing a single threshold that could be misread.

## 11. Open questions & explicit deferrals

- **Continuous / catheter infusions and cumulative toxicity.** Repeated-dose and continuous-catheter LA (with metabolite accumulation, e.g. for lidocaine) is a richer absorption/disposition problem deferred to a later phase; LA0 covers single defined doses.
- **Acid-base and ion-trapping effects on free drug.** pH-dependent shifts in ionization and free fraction (clinically important in acidosis/arrest) are curated as caveats first; explicit modeling is deferred until sources support it without implying false precision.
- **Combined-agent / additive systemic load.** Mixtures and field blocks with multiple agents create an additive systemic toxicity question; the worst-tier/never-invent rules cover labeling, but validated additive-toxicity modeling is deferred.
- **Block efficacy at the nerve** remains explicitly out of scope (§2) — this subsystem is systemic concentration and systemic toxicity only.
- **The threshold-uncertainty problem may be irreducible.** Unlike PK uncertainty (which more data narrows), the individual variability of toxicity thresholds may be a genuine floor — in which case the honest, permanent conclusion of the §6 view is "no defensible single safe-dose number exists," and that conclusion is a feature, not a gap to be closed.
