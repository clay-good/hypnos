import pytest

import hypnos
from hypnos.cli import main
from hypnos.verification import (
    checklist_markdown,
    model_verification,
    next_to_verify,
    verification_summary,
)

SCHNIDER = "hypnotics_iv.propofol.schnider_1998"
ELEVELD = "hypnotics_iv.propofol.eleveld_2018"
FENTANYL = "opioids.fentanyl.shafer_1990"  # still kernel-pending


@pytest.fixture(scope="module")
def ds():
    return hypnos.load()


def test_summary_counts_all_models(ds):
    s = verification_summary(ds)
    assert s["n_models"] == len(ds)
    bs = s["by_review_status"]
    assert sum(bs.values()) == len(ds)
    # the new pending_human_review state is counted alongside the others
    assert "pending_human_review" in bs
    assert 0.0 <= s["verified_fraction"] <= 1.0
    # nothing is human-verified by an automated process (the governance invariant)
    assert bs.get("verified", 0) == 0


def test_checklist_covers_params_and_covariate_equations(ds):
    mv = model_verification(ds, SCHNIDER)
    groups = {it.group for it in mv.checklist}
    assert {"structural", "covariate", "envelope", "population", "citation"} <= groups
    # the James LBM equation (where transcription errors hide) is a checklist item
    assert any("LBM" in it.label or "lbm" in it.value.lower() for it in mv.checklist)
    # one covariate item per parameter that has an equation
    cov = [it for it in mv.checklist if it.group == "covariate"]
    assert any("0.0681" in it.value for it in cov)  # the Schnider Cl1 LBM coefficient


def test_checklist_surfaces_curated_source_locator(ds):
    # Schnider curates a per-parameter source locator; the checklist must surface it so
    # a verifier goes straight to the right table instead of hunting the whole PDF.
    mv = model_verification(ds, SCHNIDER)
    v1 = next(it for it in mv.checklist if it.label == "parameter V1")
    assert v1.locator == "Schnider 1998, Table 2"
    md = checklist_markdown(mv)
    assert "Schnider 1998, Table 2" in md


def test_checklist_flags_missing_source_locator(ds, capsys):
    # Where a structural parameter has no curated locator, both the API field and the
    # CLI render must flag it as a gap to fill (honest about the dataset's own provenance).
    mv = model_verification(ds, "hypnotics_iv.propofol.marsh_1991")
    structural = [it for it in mv.checklist if it.group == "structural"]
    assert any(it.locator is None for it in structural)
    assert main(["verify", "hypnotics_iv.propofol.marsh_1991"]) == 0
    out = capsys.readouterr().out
    assert "no source locator curated" in out
    assert "have no curated source locator" in out  # the summary nudge


def test_next_to_verify_prioritizes_pending_then_kernels(ds):
    ordered = next_to_verify(ds)
    # all returned models are not yet verified
    assert all(m.review_status != "verified" for m in ordered)
    # pending_human_review (the cheapest human win — confirm an automated check) leads
    statuses = [m.review_status for m in ordered]
    if "pending_human_review" in statuses and "unverified" in statuses:
        assert statuses.index("pending_human_review") < statuses.index("unverified")
    # within a status group, an implemented kernel ranks before a kernel-pending one
    pend = [m for m in ordered if m.review_status == "pending_human_review"]
    impl_idx = [i for i, m in enumerate(pend) if m.kernel_implemented]
    pending_idx = [i for i, m in enumerate(pend) if not m.kernel_implemented]
    if impl_idx and pending_idx:
        assert min(impl_idx) < min(pending_idx)


def test_kernel_pending_model_lists_blocking(ds):
    # rocuronium is still unverified AND kernel-pending: the kernel block must surface.
    mv = model_verification(ds, "nmb_agents.rocuronium.wierda_1991")
    assert mv.review_status == "unverified"
    assert any("kernel" in b.lower() for b in mv.blocking)


# --------------------------------------------------------------------------- #
# the pending_human_review mechanism (automated source cross-check, never promotion)
# --------------------------------------------------------------------------- #
from hypnos.models import Model  # noqa: E402
from hypnos.validate import _check_review_state  # noqa: E402


def test_pending_human_review_records_carry_auditable_provenance(ds):
    pend = [m for m in ds if m.review_status == "pending_human_review"]
    assert pend, "expected some models cross-checked against sources"
    for m in pend:
        sr = m.source_review
        assert sr is not None and sr["human_verified"] is False   # never human-verified
        assert sr.get("sources"), f"{m.id} pending review must cite the sources compared"
        assert sr.get("outcome") in ("match", "partial_match", "discrepancy")
        assert m.is_source_reviewed


def test_pending_review_blocking_names_the_automated_check(ds):
    pend = next(m for m in ds if m.review_status == "pending_human_review")
    mv = model_verification(ds, pend.id)
    assert any("AUTOMATED source cross-check" in b for b in mv.blocking)
    assert any("LLMs do not promote" in b or "promote to 'verified'" in b for b in mv.blocking)


def test_contested_record_documents_the_disagreement(ds):
    # the Schnider-1999 BIS attribution problem was surfaced, not silently fixed
    m = ds["pd_effect.propofol.bis_sigmoid"]
    assert m.review_status == "contested"
    assert m.source_review and m.source_review["human_verified"] is False


def _review_model(status, source_review=None):
    raw = {"id": "x.y.z", "drug": {"name": "propofol"}, "purpose": "pk", "tier": "C",
           "primary_citation": "marsh-1991-propofol-pk",
           "structure": {"compartments": 1, "parameterization": "volumes_clearances"},
           "parameters": [], "extraction": {"review_status": status}}
    if source_review is not None:
        raw["extraction"]["source_review"] = source_review
    return Model(raw)


def test_validate_rejects_pending_without_provenance():
    assert any("requires a source_review" in p
               for p in _check_review_state(_review_model("pending_human_review")))


def test_validate_rejects_machine_claiming_human_verification():
    sr = {"reviewed_by": "x", "method": "y", "outcome": "match", "sources": ["u"],
          "human_verified": True}
    assert any("human_verified=false" in p
               for p in _check_review_state(_review_model("pending_human_review", sr)))


def test_validate_rejects_pending_with_no_sources():
    sr = {"reviewed_by": "x", "method": "y", "outcome": "match", "sources": [],
          "human_verified": False}
    assert any("sources[] is empty" in p
               for p in _check_review_state(_review_model("pending_human_review", sr)))


def test_validate_accepts_clean_pending_record():
    sr = {"reviewed_by": "x", "method": "y", "outcome": "match", "sources": ["url"],
          "human_verified": False}
    assert _check_review_state(_review_model("pending_human_review", sr)) == []


def test_checklist_markdown_is_copy_pasteable(ds):
    md = checklist_markdown(model_verification(ds, SCHNIDER))
    assert md.startswith("# Verification checklist")
    assert "- [ ]" in md
    assert "review_status" in md  # tells the human what to edit


def test_cli_status(capsys):
    assert main(["status"]) == 0
    out = capsys.readouterr().out
    assert "verification coverage" in out and "start here" in out


def test_cli_verify_and_markdown(capsys):
    assert main(["verify", SCHNIDER]) == 0
    out = capsys.readouterr().out
    assert "checklist" in out and "do not promote" in out
    assert main(["verify", SCHNIDER, "--markdown"]) == 0
    md = capsys.readouterr().out
    assert md.startswith("# Verification checklist")


def test_cli_verify_unknown_model(capsys):
    assert main(["verify", "nope.nope.nope"]) == 2
