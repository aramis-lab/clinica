from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "t1-freesurfer-longitudinal"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.overwrite_outputs
@cli_param.option.atlas_path
@option.global_option_group
@option.n_procs
@click.pass_context
def cli(
    ctx: click.Context,
    caps_directory: str,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    overwrite_outputs: bool = False,
    atlas_path: Optional[str] = None,
) -> None:
    """Longitudinal pre-processing of T1w images with FreeSurfer.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer_Longitudinal/
    """
    import datetime
    import os

    from clinica.utils.longitudinal import get_participants_long_id
    from clinica.utils.participant import get_subject_session_list
    from clinica.utils.stream import cprint

    from . import t1_freesurfer_longitudinal_correction_cli, t1_freesurfer_template_cli
    from .longitudinal_utils import save_part_sess_long_ids_to_tsv

    cprint(
        "The t1-freesurfer-longitudinal pipeline is divided into 2 parts:\n"
        "\tt1-freesurfer-unbiased-template pipeline: Creation of unbiased template\n"
        "\tt1-freesurfer-longitudinal-correction pipeline: Longitudinal correction."
    )

    if not subjects_sessions_tsv:
        l_sess, l_part = get_subject_session_list(caps_directory, None, False, False)
        l_long = get_participants_long_id(l_part, l_sess)
        now = datetime.datetime.now().strftime("%H%M%S")
        subjects_sessions_tsv = now + "_participants.tsv"
        save_part_sess_long_ids_to_tsv(
            l_part, l_sess, l_long, os.getcwd(), subjects_sessions_tsv
        )

    cprint("Part 1/2: Running t1-freesurfer-unbiased-template pipeline.")
    ctx.invoke(
        t1_freesurfer_template_cli.cli,
        caps_directory=caps_directory,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
        overwrite_outputs=overwrite_outputs,
    )

    cprint("Part 2/2 Running t1-freesurfer-longitudinal-correction pipeline.")
    ctx.invoke(
        t1_freesurfer_longitudinal_correction_cli.cli,
        caps_directory=caps_directory,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
        overwrite_outputs=overwrite_outputs,
        atlas_path=atlas_path,
    )


if __name__ == "__main__":
    cli()
