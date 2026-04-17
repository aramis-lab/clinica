from pathlib import Path
from typing import Optional

import click


@click.command(name="av1451-pipeline")
@click.argument(
    "data_dir",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
)
@click.argument(
    "output_dir",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
)
@click.option(
    "--inventory-dir",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
    default=None,
    help=(
        "When provided, saves oasis3_av1451_inventory.csv and "
        "oasis3_av1451_no_usable_session.csv to this directory."
    ),
)
def cli(
    data_dir: str,
    output_dir: str,
    inventory_dir: Optional[str] = None,
) -> None:
    """AV1451 Tau-PET Run Selection Pipeline for OASIS-3.

    Classifies every PET run in DATA_DIR, logs sessions with no usable scan,
    then coregisters and averages the usable late-phase frames into OUTPUT_DIR.
    """
    from ._pipeline import run_pipeline

    run_pipeline(
        data_dir=Path(data_dir),
        output_dir=Path(output_dir),
        inventory_dir=Path(inventory_dir) if inventory_dir else None,
    )


if __name__ == "__main__":
    cli()
