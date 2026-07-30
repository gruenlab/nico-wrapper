"""Umbrella command-line interface for nico-wrapper."""

from __future__ import annotations

import typer

from nico_wrapper.covariation.cli import app as covariation_app
from nico_wrapper.niche.cli import app as niche_app
from nico_wrapper.preprocess.cli import app as preprocess_app
from nico_wrapper.transfer.cli import app as transfer_app

app = typer.Typer(
    name="nico-wrapper",
    help="Run the NiCo preprocessing, transfer, niche, and covariation workflows.",
    no_args_is_help=True,
)

app.add_typer(preprocess_app, name="preprocess")
app.add_typer(transfer_app, name="transfer")
app.add_typer(niche_app, name="niche")
app.add_typer(covariation_app, name="covariation")


if __name__ == "__main__":
    app()
