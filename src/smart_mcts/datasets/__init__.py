"""Dataset downloads for Smart MCTS, handled with :mod:`pooch`.

Data lives outside the repository by default (the OS cache directory, or
``$SMART_MCTS_DATA``). See :mod:`smart_mcts.datasets._core` for the scheme.

    >>> from smart_mcts.datasets import load_moses, load_smartchemist_annotator
    >>> smiles = load_moses("train")                 # doctest: +SKIP
    >>> annotator = load_smartchemist_annotator()    # doctest: +SKIP
"""

from smart_mcts.datasets._core import (
    DATA_HOME_ENV_VAR,
    REGISTRY,
    fetch_file,
    get_data_home,
)
from smart_mcts.datasets.moses import load_moses
from smart_mcts.datasets.smartchemist import (
    fetch_smartchemist_license,
    fetch_smartchemist_smarts,
    load_smartchemist_annotator,
)

__all__ = [
    "DATA_HOME_ENV_VAR",
    "REGISTRY",
    "fetch_file",
    "get_data_home",
    "load_moses",
    "fetch_smartchemist_license",
    "fetch_smartchemist_smarts",
    "load_smartchemist_annotator",
]
