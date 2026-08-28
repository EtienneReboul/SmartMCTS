"""Console script for smart_mcts."""

import typer
from rich.console import Console

from smart_mcts import utils

app = typer.Typer()
console = Console()


@app.command()
def main() -> None:
    """Console script for smart_mcts."""
    console.print("Replace this message by putting your code into smart_mcts.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    utils.do_something_useful()


if __name__ == "__main__":
    app()
