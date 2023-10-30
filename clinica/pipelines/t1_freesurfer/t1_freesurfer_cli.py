from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "t1-freesurfer"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option_group.pipeline_specific_options
@cli_param.option_group.option(
    "-raa",
    "--recon_all_args",
    default="-qcache",
    show_default=True,
    help=(
        "Additional flags for recon-all command line "
        "Please note that = is compulsory after --recon_all_args/-raa flag "
        "(this is not the case for other flags)."
    ),
)
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.overwrite_outputs
@cli_param.option.yes
@cli_param.option.atlas_path
@option.global_option_group
@option.n_procs
@click.pass_context
def cli(
    ctx: click.Context,
    bids_directory: str,
    caps_directory: str,
    recon_all_args: str,
    working_directory: Optional[str] = None,
    subjects_sessions_tsv: Optional[str] = None,
    n_procs: Optional[int] = None,
    overwrite_outputs: bool = False,
    yes: bool = False,
    atlas_path: Optional[str] = None,
) -> None:
    """Cross-sectional pre-processing of T1w images with FreeSurfer.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from ..t1_freesurfer_atlas import t1_freesurfer_atlas_cli
    from .t1_freesurfer_pipeline import T1FreeSurfer

    parameters = {"recon_all_args": recon_all_args, "skip_question": yes}

    pipeline = T1FreeSurfer(
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
        overwrite_caps=overwrite_outputs,
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

    if atlas_path is not None:
        ctx.invoke(
            t1_freesurfer_atlas_cli.cli,
            caps_directory=caps_directory,
            atlas_path=atlas_path,
            n_procs=n_procs,
        )


if __name__ == "__main__":
    cli()
