"""Data-home resolution and the ``pooch`` download registry.

The scheme mirrors ``skfp.datasets`` from `scikit-fingerprints
<https://github.com/scikit-fingerprints/scikit-fingerprints>`_:

* every loader takes ``data_dir=None`` and resolves it to a *data home* that
  lives **outside the repository** by default;
* each dataset gets its own sub-directory under the data home;
* a file that is already present is not downloaded again unless
  ``force_update=True``;
* every download is checksum-verified (SHA-256) against a pinned registry, so a
  truncated or tampered file is rejected.

Resolution order for the data home (first hit wins):

1. an explicit ``data_dir`` argument;
2. the ``$SMART_MCTS_DATA`` environment variable;
3. the per-user OS cache directory, ``pooch.os_cache("smart_mcts")``
   (``~/Library/Caches/smart_mcts`` on macOS, ``~/.cache/smart_mcts`` on Linux,
   ``%LOCALAPPDATA%\\smart_mcts\\Cache`` on Windows).
"""

from __future__ import annotations

import os
from pathlib import Path

import pooch

#: Environment variable that overrides the default data-home location.
DATA_HOME_ENV_VAR = "SMART_MCTS_DATA"

# Upstream revisions are pinned so the registry hashes stay valid. Bump these
# (and the matching hashes below) to pull newer data.
_SMARTCHEMIST_REV = "55268b4be9a66b9ca583ca5e82507689677c533b"
_MOSES_REV = "dd7ed6ab38e23afd3ef5371d67939a1760bd8599"

_SMARTCHEMIST_BASE = f"https://raw.githubusercontent.com/torbengutermuth/SmartChemist/{_SMARTCHEMIST_REV}/smarts/"
# MOSES ships its CSVs through Git LFS; the ``media`` host serves the real bytes.
_MOSES_BASE = f"https://media.githubusercontent.com/media/molecularsets/moses/{_MOSES_REV}/data/"

#: ``pooch`` registry: relative path under the data home -> ``"sha256:<hex>"``.
REGISTRY: dict[str, str] = {
    "smartchemist/functional_groups.csv": "sha256:d88f10346dc7ff18891501ce1c200eb4e64f9e3d2c1d10c940e2f1c373a49b03",
    "smartchemist/cyclic.csv": "sha256:db9b67aa285a8e026fcb12ffeb412941ad125e3a1600236dbe1e795dd1727d17",
    "smartchemist/biologicals.csv": "sha256:8a6040ba516e294afa5725ddeee92a616349c36fd0c0c252d41d0b138e9853fd",
    "smartchemist/smarts_with_hierarchy.csv": "sha256:16455de171ac891fc4e177ab78e614eb7e2fcecafe00d28ce5bdd529c22c0ba1",
    # CC-BY-ND 4.0 license + attribution text for the SMARTS pattern collection.
    "smartchemist/License_for_patterns_here": "sha256:5e9a75ced22706937aa22109e7ec6e1cbb5fc68138dac8bedc61e488e09cb82e",
    "moses/dataset_v1.csv": "sha256:bb47a94d347afd476d3828b5e26dceeabc42a2d8cf92a791d00349f22fea0d8b",
}

#: Per-file download URLs (base URLs differ between the two upstreams).
URLS: dict[str, str] = {
    "smartchemist/functional_groups.csv": _SMARTCHEMIST_BASE + "functional_groups.csv",
    "smartchemist/cyclic.csv": _SMARTCHEMIST_BASE + "cyclic.csv",
    "smartchemist/biologicals.csv": _SMARTCHEMIST_BASE + "biologicals.csv",
    "smartchemist/smarts_with_hierarchy.csv": _SMARTCHEMIST_BASE + "smarts_with_hierarchy.csv",
    "smartchemist/License_for_patterns_here": _SMARTCHEMIST_BASE + "License_for_patterns_here",
    "moses/dataset_v1.csv": _MOSES_BASE + "dataset_v1.csv",
}


def get_data_home(data_dir: str | os.PathLike | None = None) -> Path:
    """Return the data-home directory, creating it if needed.

    See the module docstring for the resolution order. The directory (and any
    missing parents) is created before returning.
    """
    if data_dir is not None:
        home = Path(data_dir).expanduser()
    elif os.environ.get(DATA_HOME_ENV_VAR):
        home = Path(os.environ[DATA_HOME_ENV_VAR]).expanduser()
    else:
        home = Path(pooch.os_cache("smart_mcts"))

    home.mkdir(parents=True, exist_ok=True)
    return home


def _make_pooch(data_home: Path) -> pooch.Pooch:
    """Build a :class:`pooch.Pooch` rooted at ``data_home``."""
    return pooch.create(
        path=data_home,
        base_url="",
        registry=REGISTRY,
        urls=URLS,
        env=DATA_HOME_ENV_VAR,
    )


def fetch_file(
    key: str,
    data_dir: str | os.PathLike | None = None,
    *,
    verbose: bool = False,
    force_update: bool = False,
) -> Path:
    """Download (once) and return the local path of a registered data file.

    Parameters
    ----------
    key
        A key of :data:`REGISTRY`, e.g. ``"moses/dataset_v1.csv"``.
    data_dir
        Root data directory; see :func:`get_data_home`.
    verbose
        Show a download progress bar.
    force_update
        Re-download even if the file is already present and valid.
    """
    data_home = get_data_home(data_dir)
    pup = _make_pooch(data_home)

    if force_update:
        (data_home / key).unlink(missing_ok=True)

    # HTTPDownloader is pooch's own default downloader; its inline types are
    # stricter than the class satisfies, hence the ignore.
    downloader = pooch.HTTPDownloader(progressbar=verbose)
    return Path(pup.fetch(key, downloader=downloader))  # ty: ignore[invalid-argument-type]
