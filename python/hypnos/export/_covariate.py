"""Project the v0.7 covariate-equation layer into each export's native convention (C3).

The derived-covariate equations (LBM/FFM) are *part of the model definition*, so they
belong inside every export — not as an afterthought. This module turns a curated
``covariate_equations`` record's verbatim ``sex_specific`` form into:

* explicit **NONMEM ``$PK``** computations (so the control stream computes LBM the way the
  model was derived, named by equation id — the most consequential export fix, v0.7 §8);
* **``hypnos:`` RDF** predicates (``hypnos:covariateEquation`` + validity envelope) for the
  SBML annotation, plus an SBML ``assignmentRule`` that *computes* the derived covariate,
  so a deterministic consumer sees the right transform and the envelope warning rides along;
* a first-class **PharmML** ``CovariateTransformation`` block, and a named **R / Julia
  function** for the nlmixr2/rxode2 and Pumas companions.

The expression translator is a faithful, tested token substitution of the same library
form ``hypnos.covariates.evaluate`` uses, so the exported computation round-trips against
the kernel (v0.7 §8: "round-trip equality of the exported computation vs. the library").
"""
from __future__ import annotations

import re
import ast
from typing import Any, Dict, Iterator, List, Optional, Tuple
from xml.sax.saxutils import escape, quoteattr

# Library covariate variable -> the NONMEM data column it maps to.
_VAR_TO_COL = [("height_cm", "HT"), ("weight", "WT"), ("age", "AGE"), ("bmi", "BMI"), ("sex", "SEX")]
# Which derived quantities become which NONMEM-friendly symbol.
_QUANTITY_SYM = {"lbm": "LBM", "ffm": "FFM", "ibw": "IBW", "nfm": "NFM", "bsa": "BSA"}


def to_nonmem_expr(expr: str) -> str:
    """Translate a library covariate expression to NONMEM/Fortran infix.

    ``weight``->``WT``, ``height_cm``->``HT``, ``age``->``AGE``, ``bmi``->``BMI``, ``^``->``**``.
    Word-boundary substitution (longest names first) so ``height_cm`` is replaced before any
    bare token could shadow it. The arithmetic is otherwise identical to the curated form."""
    s = expr
    for var, col in _VAR_TO_COL:
        s = re.sub(rf"\b{re.escape(var)}\b", col, s)
    return s.replace("^", "**")


def _sym_for(quantity: str) -> str:
    return _QUANTITY_SYM.get(quantity, quantity.upper())


def nonmem_covariate_pk(model, ds) -> Tuple[List[str], List[str]]:
    """Return ``(pk_lines, extra_input_cols)`` computing each derived covariate in ``$PK``.

    Empty when the model declares no ``covariate_model`` (an explicit gap, never an assumed
    equation). The lines are a named, sex-conditional computation of LBM/FFM/etc. from the
    input columns, with the equation id, validity envelope, and ``used_for`` parameters in
    comments — so the control stream cannot silently recompute the body-size descriptor a
    different way (v0.7 §1)."""
    cm = model.covariate_model
    if cm is None or ds is None:
        return [], []
    eq_lib = getattr(ds, "covariate_equations", {})
    lines: List[str] = []
    needed_cols: set = set()
    bmi_emitted = False
    for di in cm.derived_inputs:
        rec = eq_lib.get(di.equation)
        if not rec:
            continue
        sexspec = rec.get("sex_specific") or {}
        male = sexspec.get("male")
        female = sexspec.get("female")
        if not (male and female):
            continue
        sym = _sym_for(di.quantity)
        env = rec.get("validity_envelope") or {}
        env_str = "; ".join(f"{ax} [{r.get('min')}, {r.get('max')}]"
                            for ax, r in env.items() if isinstance(r, dict)) or "none"
        nm_male, nm_female = to_nonmem_expr(male), to_nonmem_expr(female)
        for col in ("WT", "HT", "AGE"):
            if col in nm_male or col in nm_female:
                needed_cols.add(col)
        needs_bmi = "BMI" in nm_male or "BMI" in nm_female
        lines.append(f"; --- derived covariate {sym} via {di.equation} "
                     f"(v0.7 covariate equation; tier {di.tier}, verbatim={str(di.verbatim).lower()}) ---")
        lines.append(f";     validity envelope: {env_str}; OUTSIDE -> documented failure mode "
                     "(e.g. James LBM inversion) -> Tier D")
        if needs_bmi and not bmi_emitted:
            lines.append("  BMI = WT/(HT/100)**2")
            needed_cols.update(("WT", "HT"))
            bmi_emitted = True
        lines.append("  IF (SEX.EQ.1) THEN")
        lines.append(f"    {sym} = {nm_male}   ; male ({di.equation})")
        lines.append("  ELSE")
        lines.append(f"    {sym} = {nm_female}   ; female ({di.equation})")
        lines.append("  ENDIF")
        lines.append(f";     used_for: {', '.join(di.used_for)}  (these parameters scale on {sym})")
    return lines, sorted(needed_cols)


def _derived(model, ds) -> Iterator[Tuple[Any, Dict[str, Any], str, str, str]]:
    """Yield ``(binding, record, male_expr, female_expr, envelope_text)`` for each derived
    covariate the model binds to a library equation with a sex-specific form."""
    cm = model.covariate_model
    if cm is None or ds is None:
        return
    eq_lib = getattr(ds, "covariate_equations", {})
    for di in cm.derived_inputs:
        rec = eq_lib.get(di.equation)
        sexspec = (rec or {}).get("sex_specific") or {}
        if not (sexspec.get("male") and sexspec.get("female")):
            continue
        env = rec.get("validity_envelope") or {}
        env_str = "; ".join(f"{ax} [{r.get('min')}, {r.get('max')}]"
                            for ax, r in env.items() if isinstance(r, dict)) or "none"
        yield di, rec, sexspec["male"], sexspec["female"], env_str


def covariate_status(model, ds, patient: Optional[Dict[str, Any]]) -> Tuple[Optional[str], List[str]]:
    """``(covariate band tier, warnings)`` for the model's bound equations AT ``patient``.

    The tier is the worst the bound equations contribute for this individual (an equation
    outside its own validity envelope forces D); the warnings are the equation's own
    failure-mode text, so an export instantiated for such a patient cannot reproduce the
    James inversion unknowingly (v0.7 §8). ``(None, [])`` without a covariate_model/patient."""
    from ..covariates import EQUATIONS, covariate_layer_tier, evaluate

    cm = model.covariate_model
    if cm is None or ds is None or not patient:
        return None, []
    try:
        tier = covariate_layer_tier(model, patient, ds)
        warns = [w for d in cm.derived_inputs if d.equation in EQUATIONS
                 for w in evaluate(d.equation, patient, ds=ds).warnings]
    except (KeyError, TypeError, ValueError):      # patient lacks a covariate the equation needs
        return None, []
    return tier, warns


def _function_expr(expr: str) -> str:
    """Library expression with the data-column names (``^`` kept — R and Julia both use it)."""
    s = expr
    for var, col in _VAR_TO_COL:
        s = re.sub(rf"\b{re.escape(var)}\b", col, s)
    return s


def function_covariate_block(model, ds, name: str, lang: str) -> List[str]:
    """A named R (``lang="r"``) or Julia (``lang="julia"``) function per derived covariate.

    The nlmixr2/rxode2 and Pumas companions are instantiated at one patient, so the
    transform would otherwise be invisible: this emits it as a callable — the equation id,
    tier, and validity envelope in comments — so a script that re-instantiates the model for
    another patient computes LBM/FFM the way the model was derived (v0.7 §8). Empty when the
    model declares no covariate_model."""
    lines: List[str] = []
    for di, _rec, male, female, env_str in _derived(model, ds):
        fn = f"{name}_{di.quantity}"
        m, f = _function_expr(male), _function_expr(female)
        if not lines:
            lines += ["", "# ── covariate equations (v0.7): the body-size transform this model was derived with ──",
                      "#    Covariate-model facts — distinct from Omega (between-subject) and from the",
                      "#    estimate covariance. SEX: 1 = male, otherwise female. WT kg, HT cm, AGE years."]
        lines.append(f"# {_sym_for(di.quantity)} via {di.equation} (tier {di.tier}, "
                     f"verbatim={str(di.verbatim).lower()}); used for: {', '.join(di.used_for)}")
        lines.append(f"#    validity envelope: {env_str}; OUTSIDE it -> documented failure mode -> Tier D")
        if lang == "r":
            lines += [f"{fn} <- function(WT, HT, SEX, AGE = NA) {{",
                      "  BMI <- WT / (HT / 100)^2",
                      f"  if (SEX == 1) {m} else {f}",
                      "}"]
        else:
            lines += [f"function {fn}(WT, HT, SEX; AGE = NaN)",
                      "    BMI = WT / (HT / 100)^2",
                      f"    return SEX == 1 ? ({m}) : ({f})",
                      "end"]
    return lines


def pharmml_covariate_xml(model, ds, patient: Optional[Dict[str, Any]] = None) -> str:
    """First-class PharmML ``CovariateTransformation`` blocks (v0.7 §8): each derived
    covariate with its named equation, both sex branches, its validity envelope, and its
    source — the durable anchor carries the transform, not just θ. ``""`` without a
    covariate_model."""
    parts: List[str] = []
    for di, rec, male, female, _ in _derived(model, ds):
        cit = ds.citation(rec.get("primary_citation")) or {}
        src = f' sourceDoi={quoteattr(cit["doi"])}' if cit.get("doi") else ""
        parts.append(
            f'      <CovariateTransformation quantity="{escape(di.quantity)}" '
            f'equation="{escape(di.equation)}" verbatim="{str(di.verbatim).lower()}" '
            f'usedFor={quoteattr(" ".join(di.used_for))} tier="{di.tier}"{src}>')
        parts.append(f'        <Branch sex="male">{escape(_function_expr(male))}</Branch>')
        parts.append(f'        <Branch sex="female">{escape(_function_expr(female))}</Branch>')
        for ax, r in (rec.get("validity_envelope") or {}).items():
            if isinstance(r, dict):
                parts.append(f'        <ValidityEnvelope axis="{escape(ax)}" '
                             f'min="{r.get("min")}" max="{r.get("max")}"/>')
        parts.append('      </CovariateTransformation>')
    if not parts:
        return ""
    tier, warns = covariate_status(model, ds, patient)
    head = (f'      <CovariateLayer covariateSensitivityStatus='
            f'{quoteattr(str(model.covariate_sensitivity_status))}'
            + (f' covariateBandTier="{tier}"' if tier else "") + '/>')
    tail = [f'      <CovariateWarning>{escape(w)}</CovariateWarning>' for w in warns]
    return "\n" + "\n".join([head] + parts + tail)


def _mathml(node: ast.AST, pad: str) -> List[str]:
    """Content-MathML for an arithmetic expression AST (numbers, names, + - * / **, unary -)."""
    if isinstance(node, ast.Constant):
        return [f"{pad}<cn>{node.value!r}</cn>"]
    if isinstance(node, ast.Name):
        return [f"{pad}<ci>{node.id}</ci>"]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return [f"{pad}<apply><minus/>"] + _mathml(node.operand, pad + "  ") + [f"{pad}</apply>"]
    op = {ast.Add: "plus", ast.Sub: "minus", ast.Mult: "times", ast.Div: "divide",
          ast.Pow: "power"}[type(node.op)]
    return ([f"{pad}<apply><{op}/>"] + _mathml(node.left, pad + "  ")
            + _mathml(node.right, pad + "  ") + [f"{pad}</apply>"])


def sbml_covariate_rules(model, ds, patient: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """``(parameter_lines, rule_lines)`` making each derived covariate an SBML
    ``assignmentRule`` over the input covariates (v0.7 §8), so a deterministic consumer sees
    the transform itself — a sex ``piecewise`` over the same library form the kernel uses.

    The inputs (WT/HT/AGE/SEX) are emitted as constant parameters at the instantiated patient;
    the structural rate constants stay instantiated numbers, so the dynamics are unchanged.
    Empty when the model declares no covariate_model or the patient lacks weight/height."""
    derived = list(_derived(model, ds))
    if not derived or patient.get("weight") is None or patient.get("height") is None:
        return [], []
    w, h = float(patient["weight"]), float(patient["height"])
    env = {"WT": w, "HT": h, "AGE": float(patient.get("age", 35.0)), "BMI": w / (h / 100.0) ** 2,
           "SEX": 1.0 if str(patient.get("sex", "M")).upper().startswith("M") else 0.0}
    labels = {"WT": "weight (kg)", "HT": "height (cm)", "AGE": "age (years)",
              "BMI": "body-mass index (kg/m2)", "SEX": "sex (1 = male, 0 = female)"}
    params = [f'      <parameter id="{k}" name="{labels[k]}" value="{v:.10g}" constant="true"/>'
              for k, v in env.items()]
    rules: List[str] = []
    for di, _rec, male, female, _ in derived:
        sym = _sym_for(di.quantity)
        branch = {s: ast.parse(to_nonmem_expr(e), mode="eval").body
                  for s, e in (("male", male), ("female", female))}
        value = evaluate_nonmem_expr(to_nonmem_expr(male if env["SEX"] == 1.0 else female), patient)
        params.append(f'      <parameter id="{sym}" name="{di.quantity} via {di.equation}" '
                      f'value="{value:.10g}" constant="false"/>')
        rules += [f'      <assignmentRule variable="{sym}">',
                  '        <math xmlns="http://www.w3.org/1998/Math/MathML">',
                  "          <piecewise>", "            <piece>"]
        rules += _mathml(branch["male"], "              ")
        rules += ["              <apply><eq/><ci>SEX</ci><cn>1</cn></apply>", "            </piece>",
                  "            <otherwise>"]
        rules += _mathml(branch["female"], "              ")
        rules += ["            </otherwise>", "          </piecewise>", "        </math>",
                  "      </assignmentRule>"]
    return params, rules


def covariate_rdf(model, ds, indent: str = "  ", patient: Optional[Dict[str, Any]] = None) -> str:
    """``hypnos:`` RDF predicates for the covariate-equation layer (v0.7 §8).

    Carries the named equation, its verbatim flag, its validity envelope, and the model's
    covariate-sensitivity rollup as metadata so they travel with an SBML model. With a
    ``patient`` it also states the ``covariateBandTier`` at that individual and any
    out-of-envelope failure-mode warning. Returns ``""`` when the model declares no
    covariate_model."""
    cm = model.covariate_model
    if cm is None or ds is None:
        return ""
    eq_lib = getattr(ds, "covariate_equations", {})
    pad = f"{indent}    "
    lines = [
        f"{pad}<hypnos:covariateSensitivityStatus>{model.covariate_sensitivity_status}"
        "</hypnos:covariateSensitivityStatus>",
    ]
    tier, warns = covariate_status(model, ds, patient)
    if tier:
        lines.append(f"{pad}<hypnos:covariateBandTier>{tier}</hypnos:covariateBandTier>")
    lines += [f"{pad}<hypnos:covariateWarning>{escape(w)}</hypnos:covariateWarning>" for w in warns]
    for di in cm.derived_inputs:
        rec = eq_lib.get(di.equation) or {}
        form = (rec.get("form") or "").replace("&", "and")
        lines.append(
            f'{pad}<hypnos:covariateEquation hypnos:quantity="{di.quantity}" '
            f'hypnos:equation="{di.equation}" hypnos:verbatim="{str(di.verbatim).lower()}" '
            f'hypnos:usedFor="{" ".join(di.used_for)}" hypnos:tier="{di.tier}">')
        if form:
            lines.append(f'{pad}  <hypnos:form>{form}</hypnos:form>')
        env = rec.get("validity_envelope") or {}
        for ax, r in env.items():
            if isinstance(r, dict):
                lines.append(
                    f'{pad}  <hypnos:covariateValidityEnvelope hypnos:axis="{ax}" '
                    f'hypnos:min="{r.get("min")}" hypnos:max="{r.get("max")}"/>')
        lines.append(f"{pad}</hypnos:covariateEquation>")
    return "\n".join(lines)


def evaluate_nonmem_expr(expr_nonmem: str, patient: Dict[str, Any]) -> Optional[float]:
    """Numerically evaluate a translated NONMEM covariate expression for a patient.

    Used by the round-trip test: the emitted ``$PK`` computation must equal the library's
    ``covariates.evaluate`` value. Maps the data columns back to patient covariates and
    evaluates the (arithmetic-only) infix with ``**`` for power."""
    w = float(patient["weight"])
    h = float(patient["height"])
    env = {"WT": w, "HT": h, "AGE": float(patient.get("age", 35.0)),
           "BMI": w / ((h / 100.0) ** 2)}
    py = expr_nonmem
    for col in ("BMI", "AGE", "HT", "WT"):
        py = re.sub(rf"\b{col}\b", repr(env[col]), py)
    return float(eval(py, {"__builtins__": {}}, {}))  # noqa: S307 - arithmetic-only, no names
