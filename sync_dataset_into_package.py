#!/usr/bin/env python3
"""Copy the root-level ``dataset/`` (single source of truth) into the package as
``python/hypnos/dataset/`` so wheels carry the data.

For editable installs and running from source this is unnecessary — the loader
walks up to find the repo-root ``dataset/`` automatically. Run this only when
building a distributable wheel.

    python sync_dataset_into_package.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "dataset"
DST = ROOT / "python" / "hypnos" / "dataset"


def main() -> None:
    if not (SRC / "models").is_dir():
        raise SystemExit(f"no dataset found at {SRC}")
    if DST.exists():
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    n = len(list((DST / "models").glob("*.json")))
    print(f"synced {SRC} -> {DST} ({n} model records)")


if __name__ == "__main__":
    main()
