from typing import Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "dwi-connectome"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.option_group.pipeline_specific_options
@cli_param.option_group.option(
    "-nt",
    "--n_tracks",
    type=int,
    default=1e6,
    show_default=True,
    help="Set the desired number of streamlines to generate the tractography and connectome.",
)
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
def cli(
    caps_directory: str,
    n_tracks: int = 1e6,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Connectome-based processing of DWI datasets.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DWI_Connectome/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .dwi_connectome_pipeline import DwiConnectome

    parameters = {"n_tracks": n_tracks}

    pipeline = DwiConnectome(
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
    )

    exec_pipeline = (
        pipeline.run(plugin="MultiProc", plugin_args={"n_procs": n_procs})
        if n_procs
        else pipeline.run()
    )

    if isinstance(exec_pipeline, Graph):
        print_end_pipeline(
            pipeline_name, pipeline.base_dir, pipeline.base_dir_was_specified
        )


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
