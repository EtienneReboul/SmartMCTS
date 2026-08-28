"""Loader for the MOSES benchmark set (``molecularsets/moses``)."""

from __future__ import annotations

import os

import pandas as pd

from smart_mcts.datasets._core import fetch_file

_MOSES_KEY = "moses/dataset_v1.csv"
_SUBSETS = ("train", "test", "test_scaffolds")


def load_moses(
    subset: str | None = None,
    data_dir: str | os.PathLike | None = None,
    *,
    as_frame: bool = False,
    verbose: bool = False,
    force_update: bool = False,
) -> pd.DataFrame | list[str]:
    """Load the MOSES drug-like molecule set.

    Downloads ``dataset_v1.csv`` (~84 MB, ~1.94 M molecules) on first use and
    caches it under the data home. The file has two columns, ``SMILES`` and
    ``SPLIT`` (``train`` / ``test`` / ``test_scaffolds``).

    Parameters
    ----------
    subset
        One of ``"train"``, ``"test"``, ``"test_scaffolds"`` to keep only that
        split, or ``None`` for the whole set.
    data_dir
        Root data directory. ``None`` uses ``$SMART_MCTS_DATA`` or the OS cache
        directory. See :func:`smart_mcts.datasets.get_data_home`.
    as_frame
        If ``True`` return the ``DataFrame`` (``SMILES``, ``SPLIT``); otherwise
        return the ``SMILES`` column as a list of strings.
    verbose
        Show a download progress bar.
    force_update
        Re-download even if the cached file is present and valid.
    """
    if subset is not None and subset not in _SUBSETS:
        raise ValueError(f"subset must be None or one of {_SUBSETS}, got {subset!r}")

    path = fetch_file(_MOSES_KEY, data_dir, verbose=verbose, force_update=force_update)
    df = pd.read_csv(path)

    if subset is not None:
        df = df.loc[df["SPLIT"] == subset].reset_index(drop=True)

    return df if as_frame else df["SMILES"].tolist()
