"""Loaders for the SmartChemist SMARTS pattern library (``torbengutermuth/SmartChemist``)."""

from __future__ import annotations

import os
from pathlib import Path

from smart_mcts.datasets._core import fetch_file, get_data_home

#: logical name -> registry key
_SMARTS_FILES: dict[str, str] = {
    "functional_groups": "smartchemist/functional_groups.csv",
    "cyclic": "smartchemist/cyclic.csv",
    "biologicals": "smartchemist/biologicals.csv",
    "hierarchy": "smartchemist/smarts_with_hierarchy.csv",
}

#: The SMARTS collection is CC-BY-ND 4.0; this is the upstream license +
#: attribution text, always fetched alongside the patterns.
_LICENSE_FILE = "smartchemist/License_for_patterns_here"

_INDEX_CACHE_NAME = "smartchemist/annotator_index.pkl"


def fetch_smartchemist_license(
    data_dir: str | os.PathLike | None = None,
    *,
    verbose: bool = False,
    force_update: bool = False,
) -> Path:
    """Download (once) the CC-BY-ND 4.0 license/attribution text for the SMARTS collection."""
    return fetch_file(_LICENSE_FILE, data_dir, verbose=verbose, force_update=force_update)


def fetch_smartchemist_smarts(
    data_dir: str | os.PathLike | None = None,
    *,
    verbose: bool = False,
    force_update: bool = False,
) -> dict[str, Path]:
    """Download (once) the SmartChemist CSVs and return their local paths.

    Keys: ``functional_groups``, ``cyclic``, ``biologicals``, ``hierarchy``, and
    ``license`` (the CC-BY-ND 4.0 text, always fetched with the patterns).
    """
    paths = {
        name: fetch_file(key, data_dir, verbose=verbose, force_update=force_update)
        for name, key in _SMARTS_FILES.items()
    }
    paths["license"] = fetch_smartchemist_license(data_dir, verbose=verbose, force_update=force_update)
    return paths


def load_smartchemist_annotator(
    data_dir: str | os.PathLike | None = None,
    *,
    cache: bool = True,
    rebuild: bool = False,
    verbose: bool = False,
    force_update: bool = False,
):
    """Build a ready :class:`~smart_mcts.smartchemist_annotator.SmartChemistAnnotator`.

    Fetches the SMARTS CSVs, compiles the ~41 k-pattern library once (~30 s), and
    caches the compiled index (pickle) next to the CSVs so later calls are fast.

    Parameters
    ----------
    data_dir
        Root data directory; see :func:`smart_mcts.datasets.get_data_home`.
    cache
        Read/write the compiled index at ``<data_home>/smartchemist/annotator_index.pkl``.
    rebuild
        Recompile and overwrite the cached index even if it exists.
    verbose
        Show download progress bars.
    force_update
        Re-download the CSVs even if cached; implies ``rebuild``.
    """
    # imported here to keep ``pooch`` off the import path of the core annotator
    from smart_mcts.smartchemist_annotator import SmartChemistAnnotator

    paths = fetch_smartchemist_smarts(data_dir, verbose=verbose, force_update=force_update)
    index_path = get_data_home(data_dir) / _INDEX_CACHE_NAME

    if cache and index_path.exists() and not (rebuild or force_update):
        return SmartChemistAnnotator.load(index_path)

    annotator = SmartChemistAnnotator.from_csvs(
        functional_groups=paths["functional_groups"],
        cyclic=paths["cyclic"],
        biologicals=paths["biologicals"],
        hierarchy_csv=paths["hierarchy"],
    )
    if cache:
        index_path.parent.mkdir(parents=True, exist_ok=True)
        annotator.save(index_path)
    return annotator
