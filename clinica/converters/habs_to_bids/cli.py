from pathlib import Path

import click


@click.command(name="habs-to-bids")
@click.argument(
    "sourcedata",
    type=click.Path(exists=True, resolve_path=True, path_type=Path),
)
@click.argument(
    "rawdata",
    type=click.Path(writable=True, resolve_path=True, path_type=Path),
)
def cli(sourcedata: str, rawdata: str) -> None:
    """HABS to BIDS converter."""
    from ._converter import convert

    convert(sourcedata, rawdata)


if __name__ == "__main__":
    cli()
