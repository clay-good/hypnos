"""v0.7 C3 — covariate-aware exports.

The derived-covariate equation is part of the model definition, so every export must
carry it: NONMEM computes LBM/FFM the model's way (named, enveloped), the emitted
computation round-trips against the library, SBML carries it as RDF, and TCI-JSON passes
the structured block through verbatim. A model with no covariate_model emits no block.
"""
import json
import xml.dom.minidom as minidom

import pytest

import hypnos
from hypnos.export import nonmem, sbml, tci_json
from hypnos.export._covariate import evaluate_nonmem_expr, to_nonmem_expr

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"   # LBM via james_1976 on Cl1
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"     # FFM via al_sallami_2015 on V3
MARSH = "hypnotics_iv.propofol.marsh_1991"         # no covariate_model


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_nonmem_computes_the_named_covariate(ds):
    ctl = nonmem.build(ds[SCHNIDER], ds)
    assert "LBM via james_1976" in ctl                  # named by equation id
    assert "IF (SEX.EQ.1) THEN" in ctl                  # sex-conditional, not a baked number
    assert "validity envelope: bmi_kg_m2" in ctl        # the envelope rides along
    assert "WT" in ctl.split("$INPUT")[1].split("\n")[0]  # input columns added


def test_no_covariate_model_emits_no_block(ds):
    ctl = nonmem.build(ds[MARSH], ds)
    assert "derived covariate" not in ctl               # an explicit gap, never a fabricated equation


@pytest.mark.parametrize("mid", [SCHNIDER, ELEVELD,
                                 "opioids.remifentanil.minto_1997",
                                 "opioids.remifentanil.kim_2017",
                                 "opioids.remifentanil.eleveld_2017"])
def test_emitted_expression_round_trips_against_library(ds, mid):
    # the exported NONMEM computation must equal hypnos.covariates.evaluate (v0.7 §8).
    m = ds[mid]
    patient = {"weight": 80.0, "height": 175.0, "age": 45.0, "sex": "M"}
    for di in m.covariate_model.derived_inputs:
        rec = ds.covariate_equations[di.equation]
        male = to_nonmem_expr(rec["sex_specific"]["male"])
        emitted = evaluate_nonmem_expr(male, patient)
        lib = hypnos.evaluate_covariate(di.equation, patient).value
        assert emitted == pytest.approx(lib, rel=1e-9), f"{mid}/{di.equation}"


def test_sbml_carries_covariate_rdf_and_stays_well_formed(ds):
    xml = sbml.build(ds[SCHNIDER], ds)
    assert 'hypnos:covariateEquation hypnos:quantity="lbm"' in xml
    assert "hypnos:covariateValidityEnvelope" in xml
    assert "hypnos:covariateSensitivityStatus" in xml
    minidom.parseString(xml)                            # must remain valid XML


def test_tci_json_passes_covariate_model_through(ds):
    doc = json.loads(tci_json.build(ds[SCHNIDER], ds)[1]) if isinstance(
        tci_json.build(ds[SCHNIDER], ds), tuple) else tci_json.build_dict(ds[SCHNIDER], ds)
    assert doc["covariate_model"] is not None           # lossless structured block
    assert doc["covariate_sensitivity_status"] == "declared"
    di = doc["covariate_model"]["derived_inputs"][0]
    assert di["equation"] == "james_1976" and di["quantity"] == "lbm"


# --------------------------------------------------------------------------- #
# C3 completion — PharmML, nlmixr2/rxode2, Pumas, and the SBML assignmentRule (v0.7 §8)
# --------------------------------------------------------------------------- #
import re  # noqa: E402

from hypnos.export import export_model  # noqa: E402
from hypnos.export._covariate import covariate_status  # noqa: E402

_COV_MODELS = [SCHNIDER, ELEVELD, "opioids.remifentanil.minto_1997",
               "opioids.remifentanil.kim_2017", "opioids.remifentanil.eleveld_2017"]
_PATIENTS = [{"weight": 80.0, "height": 175.0, "age": 45.0, "sex": "M"},
             {"weight": 62.0, "height": 160.0, "age": 71.0, "sex": "F"}]


def _eval_infix(expr, patient):
    """Evaluate an exported R/Julia/PharmML branch (``^`` power, data-column names)."""
    return evaluate_nonmem_expr(expr.replace("^", "**"), patient)


def _branches(text, lang):
    """quantity-function name -> (male_expr, female_expr) parsed back out of an R/Julia export."""
    if lang == "r":
        pat = r"^(\w+) <- function\(WT, HT, SEX, AGE = NA\) \{\n.*\n  if \(SEX == 1\) (.*) else (.*)\n\}"
    else:
        pat = r"^function (\w+)\(WT, HT, SEX; AGE = NaN\)\n.*\n    return SEX == 1 \? \((.*)\) : \((.*)\)\nend"
    return {m.group(1): (m.group(2), m.group(3)) for m in re.finditer(pat, text, re.M)}


@pytest.mark.parametrize("mid", _COV_MODELS)
@pytest.mark.parametrize("fmt,lang", [("rxode2", "r"), ("pumas", "julia")])
def test_r_and_julia_functions_round_trip_against_library(ds, fmt, lang, mid):
    m = ds[mid]
    got = _branches(export_model(fmt, m, ds)[1], lang)
    assert len(got) == len(m.covariate_model.derived_inputs)
    for di in m.covariate_model.derived_inputs:
        male, female = got[f"{mid.replace('.', '_')}_{di.quantity}"]
        for patient in _PATIENTS:
            emitted = _eval_infix(male if patient["sex"] == "M" else female, patient)
            lib = hypnos.evaluate_covariate(di.equation, patient).value
            assert emitted == pytest.approx(lib, rel=1e-9), f"{fmt}/{mid}/{di.equation}"


@pytest.mark.parametrize("mid", _COV_MODELS)
def test_pharmml_covariate_transformation_round_trips(ds, mid):
    m = ds[mid]
    dom = minidom.parseString(export_model("pharmml", m, ds)[1])
    blocks = {b.getAttribute("equation"): b for b in dom.getElementsByTagName("CovariateTransformation")}
    assert set(blocks) == {di.equation for di in m.covariate_model.derived_inputs}
    for di in m.covariate_model.derived_inputs:
        b = blocks[di.equation]
        assert b.getAttribute("usedFor").split() == list(di.used_for)
        expr = {br.getAttribute("sex"): br.firstChild.data for br in b.getElementsByTagName("Branch")}
        for patient in _PATIENTS:
            sex = "male" if patient["sex"] == "M" else "female"
            assert _eval_infix(expr[sex], patient) == pytest.approx(
                hypnos.evaluate_covariate(di.equation, patient).value, rel=1e-9)
    layer = dom.getElementsByTagName("CovariateLayer")[0]
    assert layer.getAttribute("covariateSensitivityStatus") == m.covariate_sensitivity_status


def _eval_mathml(node, env):
    """Evaluate the content-MathML subset the SBML exporter emits."""
    kids = [k for k in node.childNodes if k.nodeType == k.ELEMENT_NODE]
    tag = node.tagName
    if tag == "cn":
        return float(node.firstChild.data)
    if tag == "ci":
        return env[node.firstChild.data.strip()]
    if tag == "piecewise":
        piece = next(k for k in kids if k.tagName == "piece")
        value, cond = [k for k in piece.childNodes if k.nodeType == k.ELEMENT_NODE]
        if _eval_mathml(cond, env):
            return _eval_mathml(value, env)
        return _eval_mathml(next(k for k in kids if k.tagName == "otherwise"), env)
    if tag in ("otherwise", "math"):
        return _eval_mathml(kids[0], env)
    op, args = kids[0].tagName, [_eval_mathml(k, env) for k in kids[1:]]
    if op == "minus":
        return -args[0] if len(args) == 1 else args[0] - args[1]
    return {"plus": lambda a, b: a + b, "times": lambda a, b: a * b, "divide": lambda a, b: a / b,
            "power": lambda a, b: a ** b, "eq": lambda a, b: a == b}[op](*args)


@pytest.mark.parametrize("mid", _COV_MODELS)
@pytest.mark.parametrize("patient", _PATIENTS)
def test_sbml_assignment_rule_computes_the_library_value(ds, mid, patient):
    m = ds[mid]
    dom = minidom.parseString(export_model("sbml", m, ds, patient)[1])
    env = {p.getAttribute("id"): float(p.getAttribute("value"))
           for p in dom.getElementsByTagName("parameter")}
    rules = {r.getAttribute("variable"): r for r in dom.getElementsByTagName("assignmentRule")}
    assert env["SEX"] == (1.0 if patient["sex"] == "M" else 0.0)
    for di in m.covariate_model.derived_inputs:
        sym = di.quantity.upper()
        lib = hypnos.evaluate_covariate(di.equation, patient).value
        math = rules[sym].getElementsByTagName("math")[0]
        # the rule's MathML, evaluated over the exported input parameters, IS the library value
        assert _eval_mathml(math, env) == pytest.approx(lib, rel=1e-9)
        assert env[sym] == pytest.approx(lib, rel=1e-9)        # and so is the stated initial value


def test_sbml_assignment_rule_leaves_the_dynamics_alone(ds):
    # the structural parameters and rate rules are unchanged — the covariate rule is additive
    dom = minidom.parseString(export_model("sbml", ds[SCHNIDER], ds)[1])
    assert len(dom.getElementsByTagName("rateRule")) == 4
    derived = {"LBM", "FFM"}
    rules = [r.getAttribute("variable") for r in dom.getElementsByTagName("assignmentRule")]
    assert [v for v in rules if v in derived] == ["LBM"]
    marsh = minidom.parseString(export_model("sbml", ds[MARSH], ds)[1])   # no covariate_model
    assert not derived & {r.getAttribute("variable")
                          for r in marsh.getElementsByTagName("assignmentRule")}


def test_covariate_band_tier_and_warning_follow_the_instantiated_patient(ds):
    # v0.7 §8: an export instantiated outside the equation's envelope says so, in the artifact
    obese = {"age": 50, "weight": 140, "height": 165, "sex": "M"}       # BMI ~51 > James' 37
    lean = {"age": 50, "weight": 77, "height": 177, "sex": "M"}
    assert covariate_status(ds[SCHNIDER], ds, lean) == ("B", [])
    tier, warns = covariate_status(ds[SCHNIDER], ds, obese)
    assert tier == "D" and any("INVERTED" in w for w in warns)
    for fmt in ("sbml", "pharmml"):
        bad = export_model(fmt, ds[SCHNIDER], ds, obese)[1]
        minidom.parseString(bad)
        assert re.search(r'covariateBandTier\W+D\b', bad) and "INVERTED" in bad
        ok = export_model(fmt, ds[SCHNIDER], ds, lean)[1]
        assert re.search(r'covariateBandTier\W+B\b', ok) and "ovariateWarning" not in ok


@pytest.mark.parametrize("fmt", ["pharmml", "rxode2", "pumas", "sbml"])
def test_no_covariate_model_emits_no_transform_anywhere(ds, fmt):
    text = export_model(fmt, ds[MARSH], ds)[1]
    for marker in ("CovariateTransformation", "covariate equations (v0.7)", 'variable="LBM"',
                   'variable="FFM"', "covariateBandTier"):
        assert marker not in text


# --------------------------------------------------------------------------- #
# v0.6 §9 / v0.5 §B3 — the drug-level protein binding travels with the model
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fmt", ["sbml", "pharmml"])
def test_xml_exports_carry_protein_binding_flags(ds, fmt):
    bup = export_model(fmt, ds["local_anesthetics.bupivacaine.systemic"], ds)[1]
    minidom.parseString(bup)
    assert re.search(r'proteinBinding (hypnos:)?fractionBound="0\.95"', bup)
    assert re.search(r'saturable="true"', bup) and re.search(r'indingSensitive="false"', bup)
    prop = export_model(fmt, ds[SCHNIDER], ds)[1]          # the v0.5 hypoalbuminemia failure mode
    assert re.search(r'fractionBound="0\.98"', prop) and re.search(r'indingSensitive="true"', prop)


def test_tci_json_carries_protein_binding_for_any_drug_that_curates_it(ds):
    assert tci_json.build_dict(ds[SCHNIDER], ds)["protein_binding"]["binding_sensitive"] is True
    assert "protein_binding" not in tci_json.build_dict(ds["volatiles.sevoflurane.mac"], ds)


@pytest.mark.parametrize("mid", _COV_MODELS)
def test_nonmem_pk_reads_only_columns_it_declares(ds, mid):
    # the sex-conditional covariate block reads SEX; it used to be missing from $INPUT, so
    # the control stream referenced an undefined variable and could not run
    ctl = nonmem.build(ds[mid], ds)
    declared = set(ctl.split("$INPUT")[1].split("\n")[0].split())
    pk = ctl.split("$PK")[1].split("$ERROR")[0]
    code = [ln.split(";")[0] for ln in pk.splitlines() if ln.strip() and not ln.lstrip().startswith(";")]
    defined = {ln.split("=")[0].strip() for ln in code if "=" in ln and "IF" not in ln.split("=")[0]}
    used = set()
    for ln in code:
        rhs = ln.split("=", 1)[1] if ("=" in ln and "IF" not in ln.split("=")[0]) else ln
        used |= set(re.findall(r"\b[A-Z][A-Z0-9]*\b", rhs))
    used -= {"IF", "THEN", "ELSE", "ENDIF", "EQ", "EXP", "THETA", "ETA"}
    assert used <= declared | defined, f"undeclared in $PK: {sorted(used - declared - defined)}"
    assert "SEX" in declared
