import click

import clinica.pydra.engine_utils as pu
import clinica.pydra.t1_linear.t1_linear as pt1
from clinica.pipelines import cli_param

pipeline_name = "pydra-t1-linear"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
) -> None:
    """Affine registration of Flair images to the MNI standard space (Pydra engine)."""

    pipeline = pt1.build_core_workflow(
        "t1-linear-pydra", bids_directory, caps_directory
    )
    pu.run(pipeline)


if __name__ == "__main__":
    cli()
