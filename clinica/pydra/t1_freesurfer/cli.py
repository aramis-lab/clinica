import click

from clinica.pipelines.engine import clinica_pipeline


@clinica_pipeline
@click.command(name="pydra-t1-freesurfer", hidden=True)
def cli() -> None:
    """Cross-sectional pre-processing of T1w volumes with FreeSurfer."""
    pass


if __name__ == "__main__":
    cli()
