from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "dwi-preprocessing-using-t1"


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
@cli_param.option.delete_cache
@cli_param.option.random_seed
@cli_param.option.caps_name
def cli(
    bids_directory: str,
    caps_directory: str,
    low_bval: int = 5,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    use_cuda: bool = False,
    initrand: bool = False,
    delete_cache: bool = False,
    random_seed: Optional[int] = None,
    caps_name: Optional[str] = None,
) -> None:
    """Preprocessing of raw DWI datasets using a T1w image.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DWI_Preprocessing/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .pipeline import DwiPreprocessingUsingT1

    parameters = {
        "low_bval": low_bval,
        "use_cuda": use_cuda,
        "initrand": initrand,
        "delete_cache": delete_cache,
        "random_seed": random_seed,
        "double_precision": True,  # TODO: expose to user ?
    }

    pipeline = DwiPreprocessingUsingT1(
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
        caps_name=caps_name,
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
