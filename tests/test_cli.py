"""Tests for the ``smart-mcts`` command-line entry point."""

from typer.testing import CliRunner

from smart_mcts import SMARTCHEMIST_NOTICE, __version__
from smart_mcts.cli import app

runner = CliRunner()


def test_version_flag_shows_version_and_attribution():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert f"smart-mcts {__version__}" in result.output
    assert "SMARTChemist Pattern Collection" in result.output
    assert "chemist.smarts.plus" in result.output


def test_attribution_notice_matches_license_wording():
    # exact sentence the SmartChemist CC-BY-ND 4.0 license asks software to show
    assert SMARTCHEMIST_NOTICE == (
        "This software uses the SMARTChemist Pattern Collection developed at the "
        "University of Hamburg, Germany. For details see https://chemist.smarts.plus."
    )
