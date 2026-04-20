from pathlib import Path
from typing import Optional

import click

from ._tracer_config import SUPPORTED_TRACERS


@click.command(name="oasis3-pet-pipeline")
@click.argument(
    "bids_dir",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
)
@click.argument(
    "output_dir",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
)
@click.option(
    "--tracer",
    type=click.Choice(SUPPORTED_TRACERS, case_sensitive=False),
    required=True,
    help="PET tracer to process (AV1451, AV45, or PIB).",
)
@click.option(
    "--inventory-dir",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
    default=None,
    help=(
        "When provided, saves per-tracer inventory and "
        "no-usable-session CSVs to this directory."
    ),
)
def cli(
    bids_dir: str,
    output_dir: str,
    tracer: str,
    inventory_dir: Optional[str] = None,
) -> None:
    """OASIS-3 PET Run Selection Pipeline.

    Scans the BIDS dataset at BIDS_DIR, classifies every PET run for the
    specified tracer, logs sessions with no usable scan, then coregisters
    and averages the usable late-phase frames into OUTPUT_DIR.
    """
    from ._pipeline import run_pipeline

    run_pipeline(
        bids_dir=Path(bids_dir),
        output_dir=Path(output_dir),
        tracer=tracer,
        inventory_dir=Path(inventory_dir) if inventory_dir else None,
    )


if __name__ == "__main__":
    cli()
