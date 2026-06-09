"""Export builders: deterministic projections of the dataset into the formats
the pharmacometrics and TCI communities run (NONMEM, PharmML, SBML, TCI-JSON).

Exports are generated artifacts, never hand-edited, and each carries the
propagated confidence tier plus the universal ``clinicalUse = PROHIBITED``
annotation (see :mod:`hypnos.export.annotate`).
"""
from __future__ import annotations

from . import annotate, bibtex, combine, nonmem, pharmml, pumas, rxode2, sbml, tci_json  # noqa: F401

# Per-model text exporters.
BUILDERS = {
    "nonmem": (nonmem.build, nonmem.filename),
    "pharmml": (pharmml.build, pharmml.filename),
    "sbml": (sbml.build, sbml.filename),
    "tci_json": (tci_json.build, tci_json.filename),
    "rxode2": (rxode2.build, rxode2.filename),
    "pumas": (pumas.build, pumas.filename),
    "bibtex": (bibtex.build_for_model, bibtex.filename),
}

# Binary / archive exporters handled separately (see hypnos.export.combine).
BINARY_FORMATS = {"omex"}

FORMATS = list(BUILDERS) + sorted(BINARY_FORMATS)


def export_model(fmt: str, model, ds=None, patient=None):
    """Return (filename, text) for one model in one text format."""
    build, fname = BUILDERS[fmt]
    return fname(model), build(model, ds, patient)


__all__ = ["BUILDERS", "BINARY_FORMATS", "FORMATS", "export_model", "annotate", "bibtex",
           "combine", "nonmem", "pharmml", "sbml", "tci_json", "rxode2", "pumas"]
