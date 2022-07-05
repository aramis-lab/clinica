import click

import clinica.pydra.t1_linear.t1_linear as pt1
from clinica.pipelines import cli_param
from clinica.pydra.engine import Pipeline

pipeline_name = "pydra-t1-linear"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
) -> None:
    """Affine registration of Flair images to the MNI standard space (Pydra engine)."""

    p1 = Pipeline(pipeline_name, bids_directory, caps_directory)
    core_workflow = pt1.build_core_workflow()
    p1.build_workflow(core_workflow)
    p1.run()


if __name__ == "__main__":
    cli()
