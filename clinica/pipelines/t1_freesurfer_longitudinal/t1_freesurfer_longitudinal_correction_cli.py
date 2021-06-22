from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "t1-freesurfer-longitudinal-correction"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.option_group.standard_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option.overwrite_outputs
def cli(
    caps_directory: str,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    overwrite_outputs: bool = False,
) -> None:
    """Longitudinal pre-processing correction of T1w images with FreeSurfer.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer_Longitudinal/"
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .t1_freesurfer_longitudinal_correction_pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
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


if __name__ == "__main__":
    cli()
