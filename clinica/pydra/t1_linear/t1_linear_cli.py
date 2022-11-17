import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_linear.t1_linear as pydra_t1_linear
from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "pydra-t1-linear"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
) -> None:
    """Affine registration of T1w images to the MNI standard space (Pydra engine)."""

    t1_linear_pipeline = pydra_t1_linear.build_core_workflow(
        name="t1-linear-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters={},
    )
    pydra_utils.run(t1_linear_pipeline)


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
