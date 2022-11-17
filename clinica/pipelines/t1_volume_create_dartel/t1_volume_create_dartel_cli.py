from typing import List, Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "t1-volume-create-dartel"


@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.dartel_tissues
def cli(
    bids_directory: str,
    caps_directory: str,
    group_label: str,
    dartel_tissues: List[int] = (1, 2, 3),
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Inter-subject registration using Dartel (creating a new Dartel template).

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects. For this pipeline, it is associated to the DARTEL template that you had created when running the t1-volume pipeline.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .t1_volume_create_dartel_pipeline import T1VolumeCreateDartel

    parameters = {"group_label": group_label, "dartel_tissues": dartel_tissues}

    pipeline = T1VolumeCreateDartel(
        bids_directory=bids_directory,
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
