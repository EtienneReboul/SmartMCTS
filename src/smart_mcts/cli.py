"""Console script for smart_mcts."""

import typer
from rich.console import Console

from smart_mcts import SMARTCHEMIST_NOTICE, __version__, utils

app = typer.Typer()
console = Console()


def _version_callback(value: bool) -> None:
    if not value:
        return
    console.print(f"smart-mcts {__version__}")
    console.print(SMARTCHEMIST_NOTICE)
    raise typer.Exit


@app.command()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-V",
        help="Show the version and third-party attribution, then exit.",
        callback=_version_callback,
        is_eager=True,
    ),
) -> None:
    """Console script for smart_mcts."""
    console.print("Replace this message by putting your code into smart_mcts.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    utils.do_something_useful()


if __name__ == "__main__":
    app()
