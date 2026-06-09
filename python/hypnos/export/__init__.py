"""Export builders: deterministic projections of the dataset into the formats
the pharmacometrics and TCI communities run (NONMEM, PharmML, SBML, TCI-JSON).

Exports are generated artifacts, never hand-edited, and each carries the
propagated confidence tier plus the universal ``clinicalUse = PROHIBITED``
annotation (see :mod:`hypnos.export.annotate`).
"""
from __future__ import annotations

from . import annotate, nonmem, pharmml, sbml, tci_json  # noqa: F401

BUILDERS = {
    "nonmem": (nonmem.build, nonmem.filename),
    "pharmml": (pharmml.build, pharmml.filename),
    "sbml": (sbml.build, sbml.filename),
    "tci_json": (tci_json.build, tci_json.filename),
}

FORMATS = list(BUILDERS)


def export_model(fmt: str, model, ds=None, patient=None):
    """Return (filename, text) for one model in one format."""
    build, fname = BUILDERS[fmt]
    return fname(model), build(model, ds, patient)


__all__ = ["BUILDERS", "FORMATS", "export_model", "annotate", "nonmem", "pharmml", "sbml", "tci_json"]
