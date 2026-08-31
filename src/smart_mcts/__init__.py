"""Top-level package for Smart MCTS."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("smart-mcts")
except PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.0.0"

#: Attribution notice required by the SmartChemist Pattern Collection license
#: (CC-BY-ND 4.0). Shown by ``smart-mcts --version``; see the project README.
SMARTCHEMIST_NOTICE = (
    "This software uses the SMARTChemist Pattern Collection developed at the "
    "University of Hamburg, Germany. For details see https://chemist.smarts.plus."
)

__all__ = ["__version__", "SMARTCHEMIST_NOTICE"]
