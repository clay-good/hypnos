"""Project the v0.7 covariate-equation layer into each export's native convention (C3).

The derived-covariate equations (LBM/FFM) are *part of the model definition*, so they
belong inside every export — not as an afterthought. This module turns a curated
``covariate_equations`` record's verbatim ``sex_specific`` form into:

* explicit **NONMEM ``$PK``** computations (so the control stream computes LBM the way the
  model was derived, named by equation id — the most consequential export fix, v0.7 §8);
* **``hypnos:`` RDF** predicates (``hypnos:covariateEquation`` + validity envelope) for the
  SBML/PharmML annotation, so a deterministic consumer sees the right transform and the
  envelope/failure-mode warning rides along.

The expression translator is a faithful, tested token substitution of the same library
form ``hypnos.covariates.evaluate`` uses, so the exported computation round-trips against
the kernel (v0.7 §8: "round-trip equality of the exported computation vs. the library").
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

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


def covariate_rdf(model, ds, indent: str = "  ") -> str:
    """``hypnos:`` RDF predicates for the covariate-equation layer (v0.7 §8).

    Carries the named equation, its verbatim flag, its validity envelope, and the model's
    covariate-sensitivity rollup as metadata so they travel with an SBML/PharmML model.
    Returns ``""`` when the model declares no covariate_model."""
    cm = model.covariate_model
    if cm is None or ds is None:
        return ""
    eq_lib = getattr(ds, "covariate_equations", {})
    pad = f"{indent}    "
    lines = [
        f"{pad}<hypnos:covariateSensitivityStatus>{model.covariate_sensitivity_status}"
        "</hypnos:covariateSensitivityStatus>",
    ]
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
