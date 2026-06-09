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
    assert sum(s["by_review_status"].values()) == len(ds)
    # every current model is unverified -> 0% coverage, and every model is "next"
    assert 0.0 <= s["verified_fraction"] <= 1.0


def test_checklist_covers_params_and_covariate_equations(ds):
    mv = model_verification(ds, SCHNIDER)
    groups = {it.group for it in mv.checklist}
    assert {"structural", "covariate", "envelope", "population", "citation"} <= groups
    # the James LBM equation (where transcription errors hide) is a checklist item
    assert any("LBM" in it.label or "lbm" in it.value.lower() for it in mv.checklist)
    # one covariate item per parameter that has an equation
    cov = [it for it in mv.checklist if it.group == "covariate"]
    assert any("0.0681" in it.value for it in cov)  # the Schnider Cl1 LBM coefficient


def test_next_to_verify_prioritizes_implemented_kernels(ds):
    ordered = next_to_verify(ds)
    # all returned models are not yet verified
    assert all(m.review_status != "verified" for m in ordered)
    # an implemented-kernel model ranks before a kernel-pending one of the same/better tier
    impl_idx = [i for i, m in enumerate(ordered) if m.kernel_implemented]
    pending_idx = [i for i, m in enumerate(ordered) if not m.kernel_implemented]
    if impl_idx and pending_idx:
        assert min(impl_idx) < min(pending_idx)


def test_kernel_pending_model_lists_blocking(ds):
    mv = model_verification(ds, FENTANYL)  # curated record, kernel still pending
    assert mv.review_status == "unverified"
    assert any("kernel" in b.lower() for b in mv.blocking)


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
