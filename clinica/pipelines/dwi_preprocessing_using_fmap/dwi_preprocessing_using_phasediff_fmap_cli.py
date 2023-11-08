from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "dwi-preprocessing-using-phasediff-fmap"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option_group.pipeline_specific_options
@cli_param.option.low_bval
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@option.global_option_group
@option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.use_cuda
@cli_param.option.initrand
def cli(
    bids_directory: str,
    caps_directory: str,
    low_bval: int = 5,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    use_cuda: bool = False,
    initrand: bool = False,
) -> None:
    """Preprocessing of raw DWI datasets using a phase difference image.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DWI_Preprocessing/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    parameters = {
        "low_bval": low_bval,
        "use_cuda": use_cuda,
        "initrand": initrand,
    }

    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
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


if __name__ == "__main__":
    cli()
