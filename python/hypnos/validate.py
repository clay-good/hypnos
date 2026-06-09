"""Validate the dataset against its JSON Schema and check cross-references.

Two layers:

* **Schema validation** — every model record conforms to ``model.schema.json``.
* **Integrity checks** — referential and semantic invariants the schema cannot
  express on its own (citations resolve, kernel bindings exist, record tier is
  the worst contributing parameter tier, envelopes are well-ordered).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

from .load import Dataset, find_dataset_dir, load
from .models import TIER_RANK, worst_tier


class ValidationError(Exception):
    pass


def _load_schema(root: Path) -> dict:
    schema_path = root / "schema" / "model.schema.json"
    with schema_path.open(encoding="utf-8") as fh:
        return json.load(fh)


def validate_dataset(ds: Optional[Dataset] = None) -> List[str]:
    """Return a list of problem strings. Empty list == dataset is valid."""
    if ds is None:
        ds = load()
    root = find_dataset_dir()
    problems: List[str] = []

    # --- schema layer -----------------------------------------------------
    try:
        import jsonschema

        schema = _load_schema(root)
        validator = jsonschema.Draft7Validator(schema)
        for m in ds:
            for err in sorted(validator.iter_errors(m.raw), key=lambda e: e.path):
                loc = "/".join(str(p) for p in err.path)
                problems.append(f"[schema] {m.id}: {loc}: {err.message}")
    except ImportError:  # pragma: no cover - jsonschema is a hard dependency
        problems.append("[schema] jsonschema not installed; skipped schema validation")

    # --- integrity layer --------------------------------------------------
    from .export.registry import (  # local import to avoid cycle
        INTERACTION_KERNELS,
        KERNELS,
        PD_KERNELS,
        VOLATILE_KERNELS,
    )

    known_kernels = set(KERNELS) | set(PD_KERNELS) | set(INTERACTION_KERNELS) | set(VOLATILE_KERNELS)

    known_citations = set(ds.citations.keys())
    for m in ds:
        # record-level citation resolves
        if m.primary_citation not in known_citations:
            problems.append(f"[cite] {m.id}: primary_citation '{m.primary_citation}' not in citations/")
        # per-parameter citations resolve when present
        for p in m.parameters:
            if p.primary_citation and p.primary_citation not in known_citations:
                problems.append(
                    f"[cite] {m.id}: parameter {p.symbol} cites unknown '{p.primary_citation}'"
                )
        # failure-mode citations resolve when present
        for fm in m.known_failure_modes:
            if fm.citation and fm.citation not in known_citations:
                problems.append(f"[cite] {m.id}: failure-mode cites unknown '{fm.citation}'")
        # predictive-performance citations resolve when present (a performance
        # number must be traceable to a source, never asserted bare — spec §5/§9)
        for pp in m.predictive_performance:
            cid = pp.get("citation")
            if cid and cid not in known_citations:
                problems.append(
                    f"[cite] {m.id}: predictive_performance ({pp.get('metric')}) cites unknown '{cid}'"
                )

        # record tier == worst contributing parameter tier (the "worst input wins" invariant)
        param_tiers = [p.tier for p in m.parameters]
        if param_tiers:
            expected = worst_tier(param_tiers)
            if TIER_RANK[m.tier] < TIER_RANK[expected]:
                problems.append(
                    f"[tier] {m.id}: record tier {m.tier} is better than worst "
                    f"parameter tier {expected} (worst input must win)"
                )

        # kernel binding resolves to a registered function
        if m.kernel_implemented:
            fn = m.kernel_function
            if fn not in known_kernels:
                problems.append(f"[kernel] {m.id}: kernel.function '{fn}' not registered")

        # envelope ranges well-ordered (min <= max)
        env = m.applicability_envelope
        for name in ("age_years", "weight_kg", "height_cm", "bmi_kg_m2"):
            rng = getattr(env, name)
            if rng.min is not None and rng.max is not None and rng.min > rng.max:
                problems.append(f"[envelope] {m.id}: {name} min {rng.min} > max {rng.max}")

        # id prefix matches declared subsystem
        if m.id.split(".")[0] != m.subsystem:
            problems.append(f"[id] {m.id}: id prefix does not match subsystem '{m.subsystem}'")

    return problems


def assert_valid(ds: Optional[Dataset] = None) -> None:
    problems = validate_dataset(ds)
    if problems:
        raise ValidationError(
            f"{len(problems)} dataset problem(s):\n  " + "\n  ".join(problems)
        )
