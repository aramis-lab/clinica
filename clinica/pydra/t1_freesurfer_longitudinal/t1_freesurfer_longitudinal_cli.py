import click

from clinica.pipelines.engine import clinica_pipeline


@clinica_pipeline
@click.command(name="pydra-t1-freesurfer-longitudinal", hidden=True)
def cli() -> None:
    """Longitudinal pre-processing of T1w images with FreeSurfer."""
    pass


if __name__ == "__main__":
    cli()
