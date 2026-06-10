"""Locate and load the Hypnos dataset.

The canonical dataset lives at the repository-root ``dataset/`` directory and is
the single source of truth. For installed wheels it may also be vendored into
``hypnos/dataset/`` by ``sync_dataset_into_package.py``. Resolution order:

1. ``$HYPNOS_DATASET`` environment variable (explicit override);
2. a packaged ``hypnos/dataset/`` copy, if present;
3. the repository-root ``dataset/`` found by walking up from this file.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

from .models import Model

_THIS = Path(__file__).resolve()


def find_dataset_dir(explicit: Optional[os.PathLike] = None) -> Path:
    """Return the directory that contains ``models/``, ``drugs/``, ``citations/``."""
    candidates: List[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit))
    env = os.environ.get("HYPNOS_DATASET")
    if env:
        candidates.append(Path(env))

    # The in-package copy (hypnos/dataset) is written by sync_dataset_into_package.py
    # for wheel builds. In a *source* checkout it can go stale, so it must NOT
    # shadow the repo-root dataset/ (the single source of truth). Prefer any
    # dataset/ found by walking up that is OUTSIDE the package; fall back to the
    # in-package copy only as a last resort (the normal case for installed wheels).
    packaged = _THIS.parent / "dataset"
    for parent in _THIS.parents:
        cand = parent / "dataset"
        if cand == packaged:
            continue
        candidates.append(cand)
    candidates.append(packaged)

    for c in candidates:
        if (c / "models").is_dir():
            return c
    raise FileNotFoundError(
        "Could not locate the Hypnos dataset directory (looked for a 'models/' "
        "subfolder). Set $HYPNOS_DATASET to point at it."
    )


def _load_json_dir(path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    if not path.is_dir():
        return out
    for f in sorted(path.glob("*.json")):
        with f.open(encoding="utf-8") as fh:
            data = json.load(fh)
        key = data.get("id", f.stem)
        out[key] = data
    return out


class Dataset:
    """In-memory view of the dataset: models, drugs, and citations."""

    def __init__(self, root: Path):
        self.root = root
        self._models_raw = _load_json_dir(root / "models")
        self.drugs = _load_json_dir(root / "drugs")
        self.citations = _load_json_dir(root / "citations")
        # Shared, citation-backed covariate-equation library (v0.7 §6); keyed by id.
        # Empty in a pre-v0.7 checkout — an absent library is not an error.
        self.covariate_equations = _load_json_dir(root / "covariate_equations")
        self._models = {mid: Model(raw) for mid, raw in self._models_raw.items()}
        self.version = self._read_version(root)

    @staticmethod
    def _read_version(root: Path) -> str:
        vfile = root / "VERSION"
        if vfile.is_file():
            return vfile.read_text(encoding="utf-8").strip()
        return "0.1.0"

    # --- dict-like access -------------------------------------------------
    def __getitem__(self, model_id: str) -> Model:
        return self._models[model_id]

    def __contains__(self, model_id: str) -> bool:
        return model_id in self._models

    def __iter__(self):
        return iter(self._models.values())

    def __len__(self) -> int:
        return len(self._models)

    @property
    def models(self) -> List[Model]:
        return list(self._models.values())

    @property
    def ids(self) -> List[str]:
        return list(self._models.keys())

    def citation(self, cid: str) -> Optional[Dict[str, Any]]:
        return self.citations.get(cid)

    def drug(self, name: str) -> Optional[Dict[str, Any]]:
        for d in self.drugs.values():
            if d.get("name") == name:
                return d
        return None

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<Dataset v{self.version} models={len(self._models)} drugs={len(self.drugs)}>"


@lru_cache(maxsize=8)
def _cached_load(root_str: str) -> Dataset:
    return Dataset(Path(root_str))


def load(path: Optional[os.PathLike] = None) -> Dataset:
    """Load the dataset. Results are cached per resolved root directory."""
    root = find_dataset_dir(path)
    return _cached_load(str(root))
