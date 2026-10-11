"""Project the v0.8 developmental extrapolation block into each export (v0.8 §7).

An adult model's ``developmental_model`` is two curated mechanisms — allometric size
scaling and (where curated) the Anderson–Holford PMA maturation of clearance. Exports
are instantiated at an adult reference patient and the block is **opt-in, Tier-D**
(``applied_by_default`` is always false), so every carrier here states the transform
without applying it:

* **NONMEM** — a commented ``$PK`` template (uncommenting it is the opt-in);
* **nlmixr2/rxode2 · Pumas** — named functions, labeled extrapolation;
* **PharmML** — a ``DevelopmentalModel`` block;
* **SBML** — ``assignmentRule``\\ s computing the scaling factors as outputs.

The expressions are the same closed forms :mod:`hypnos.developmental` evaluates
(``(WT/ref)^exp`` and ``PMA^Hill / (TM50^Hill + PMA^Hill)``), and the tests evaluate
each exported form against those kernels. The maturation clock is **PMA** everywhere,
with the substitution of chronological age forbidden in the artifact itself (Trap 1).
"""
from __future__ import annotations

import ast
from typing import Any, Dict, List, Optional, Tuple
from xml.sax.saxutils import escape, quoteattr

_PMA_RULE = ("PMA = post-menstrual age (weeks) — NEVER substitute chronological age "
             "(v0.8 Trap 1)")
_ALLOMETRY_ONLY = ("allometry_only — maturation UN-modeled: neonatal clearance OVER-stated, "
                   "exposure UNDER-predicted")


def developmental_terms(model) -> Optional[Dict[str, Any]]:
    """The curated size/maturation numbers as infix expressions, or None without a block.

    ``size_cl`` / ``size_v`` are in ``WT``; ``maturation`` (None for ``allometry_only``) is
    in ``PMA``. ``^`` is the power operator (R and Julia read it as written)."""
    dev = model.developmental_model
    if dev is None:
        return None
    terms: Dict[str, Any] = {"dev": dev, "size_cl": None, "size_v": None, "maturation": None}
    if dev.size is not None:
        ref = f"{dev.size.reference_weight_kg:g}"
        terms["size_cl"] = f"(WT/{ref})^{dev.size.cl_exponent:g}"
        terms["size_v"] = f"(WT/{ref})^{dev.size.v_exponent:g}"
    if dev.maturation is not None:
        tm50, hill = f"{dev.maturation.tm50_weeks:g}", f"{dev.maturation.hill:g}"
        terms["maturation"] = f"PMA^{hill}/({tm50}^{hill} + PMA^{hill})"
    return terms


def _header(dev) -> str:
    return (f"developmental extrapolation (v0.8): {dev.extrapolation_basis}; evidenceTier "
            f"{dev.evidence_tier}; appliedByDefault {str(dev.applied_by_default).lower()} — OPT-IN")


def nonmem_template(model) -> List[str]:
    """A commented ``$PK`` template for the developmental extrapolation (no leading ``;``
    on the returned lines — the caller comments them). Empty without a block."""
    t = developmental_terms(model)
    if t is None:
        return []
    dev = t["dev"]
    lines = [f"--- {_header(dev)} ---",
             "NOT applied: THETA below describe the adult model. To extrapolate to a child, add WT"
             + (" and PMA" if t["maturation"] else "") + " to $INPUT and uncomment."]
    cl_factor = []
    if t["size_cl"]:
        lines.append(f"  FSIZECL = {t['size_cl'].replace('^', '**')}   ; allometric CLEARANCE scaling")
        lines.append(f"  FSIZEV  = {t['size_v'].replace('^', '**')}   ; allometric VOLUME scaling")
        cl_factor.append("FSIZECL")
    if t["maturation"]:
        lines.append(f"  {_PMA_RULE}")
        lines.append(f"  FMAT    = {t['maturation'].replace('^', '**')}   ; maturation of "
                     f"{dev.maturation.affected_parameter} ({dev.maturation.function})")
        cl_factor.append("FMAT")
    else:
        lines.append(f"  CAVEAT: {_ALLOMETRY_ONLY}")
    if cl_factor:
        f = "*".join(cl_factor)
        lines.append(f"  CL = CL*{f} ; Q2 = Q2*{f} ; Q3 = Q3*{f}")
    if t["size_v"]:
        lines.append("  V1 = V1*FSIZEV ; V2 = V2*FSIZEV ; V3 = V3*FSIZEV")
    return lines


def function_block(model, name: str, lang: str) -> List[str]:
    """Named R (``lang="r"``) or Julia functions for the size + maturation factors."""
    t = developmental_terms(model)
    if t is None:
        return []
    dev = t["dev"]
    lines = ["", f"# ── {_header(dev)} ──",
             "#    A Tier-D EXTRAPOLATION of an adult model — not a fitted-pediatric model, and NOT",
             "#    applied to the parameters above. Multiply clearances by size_cl"
             + (" * maturation" if t["maturation"] else "") + ", volumes by size_v."]

    def define(suffix: str, arg: str, expr: str) -> str:
        return (f"{name}_{suffix} <- function({arg}) {expr}" if lang == "r"
                else f"{name}_{suffix}({arg}) = {expr}")

    if t["size_cl"]:
        lines.append(define("size_cl", "WT", t["size_cl"]))
        lines.append(define("size_v", "WT", t["size_v"]))
    if t["maturation"]:
        lines.append(f"#    {_PMA_RULE}")
        lines.append(define("maturation", "PMA", t["maturation"]))
    else:
        lines.append(f"#    CAVEAT: {_ALLOMETRY_ONLY}")
    return lines


def pharmml_xml(model) -> str:
    """A PharmML ``DevelopmentalModel`` block: the allometric transformation and the
    parameter-level maturation function, flagged Tier-D / opt-in (v0.8 §7)."""
    t = developmental_terms(model)
    if t is None:
        return ""
    dev = t["dev"]
    parts = [f'    <DevelopmentalModel extrapolationBasis="{escape(dev.extrapolation_basis)}" '
             f'evidenceTier="{dev.evidence_tier}" developmentalExtrapolation="true" '
             f'appliedByDefault="{str(dev.applied_by_default).lower()}">']
    if dev.size is not None:
        s = dev.size
        parts.append(
            f'      <AllometricSize clExponent="{s.cl_exponent:g}" vExponent="{s.v_exponent:g}" '
            f'referenceWeightKg="{s.reference_weight_kg:g}" exponentBasis="{escape(s.exponent_basis)}" '
            f'tier="{s.tier}" clearanceFactor={quoteattr(t["size_cl"])} '
            f'volumeFactor={quoteattr(t["size_v"])}/>')
    if dev.maturation is not None:
        m = dev.maturation
        parts.append(f'      <!-- {_PMA_RULE} -->')
        parts.append(
            f'      <MaturationFunction function="{escape(m.function)}" tm50Weeks="{m.tm50_weeks:g}" '
            f'hill="{m.hill:g}" driver="{escape(m.driver)}" '
            f'affectedParameter="{escape(m.affected_parameter)}" tier="{m.tier}" '
            f'factor={quoteattr(t["maturation"])}/>')
    else:
        parts.append(f'      <MaturationCaveat>{escape(_ALLOMETRY_ONLY)}</MaturationCaveat>')
    parts.append('    </DevelopmentalModel>')
    return "\n" + "\n".join(parts)


def _eval(expr: str, env: Dict[str, float]) -> float:
    return float(eval(expr.replace("^", "**"), {"__builtins__": {}}, dict(env)))  # noqa: S307


def sbml_rules(model, patient: Dict[str, Any], taken: Tuple[str, ...] = ()) -> Tuple[List[str], List[str]]:
    """``(parameter_lines, rule_lines)``: the scaling factors as SBML ``assignmentRule``
    outputs (v0.8 §7). They are computed, never wired into the rate constants — the dynamics
    stay the adult model's, which is what opt-in means. ``taken`` lists input ids another
    block already declared (``WT`` when the covariate rules ran)."""
    from ._covariate import _mathml

    t = developmental_terms(model)
    if t is None or patient.get("weight") is None:
        return [], []
    env = {"WT": float(patient["weight"])}
    labels = {"WT": "weight (kg)"}
    if t["maturation"]:
        pma = patient.get("pma_weeks")
        # an adult patient declares no PMA; 40 wk gestation + age is used only to give the
        # (fully mature, MF ~ 1) adult a number — a child must declare PMA explicitly.
        env["PMA"] = float(pma) if pma is not None else 40.0 + 52.0 * float(patient.get("age", 35.0))
        labels["PMA"] = "post-menstrual age (weeks) - NOT chronological age"
    params = [f'      <parameter id="{k}" name="{labels[k]}" value="{v:.10g}" constant="true"/>'
              for k, v in env.items() if k not in taken]
    rules: List[str] = []
    for pid, key, what in (("FSIZE_CL", "size_cl", "allometric clearance scaling"),
                           ("FSIZE_V", "size_v", "allometric volume scaling"),
                           ("FMAT", "maturation", "PMA clearance maturation")):
        expr = t[key]
        if expr is None:
            continue
        params.append(f'      <parameter id="{pid}" name="{what} (Tier-D extrapolation; opt-in, '
                      f'not applied)" value="{_eval(expr, env):.10g}" constant="false"/>')
        rules += [f'      <assignmentRule variable="{pid}">',
                  '        <math xmlns="http://www.w3.org/1998/Math/MathML">']
        rules += _mathml(ast.parse(expr.replace("^", "**"), mode="eval").body, "          ")
        rules += ["        </math>", "      </assignmentRule>"]
    return params, rules
