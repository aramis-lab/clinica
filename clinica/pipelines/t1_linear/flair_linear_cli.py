from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "flair-linear"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option_group.pipeline_specific_options
@cli_param.option.uncropped
@cli_param.option.random_seed
@cli_param.option.antspy
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.caps_name
@option.global_option_group
@option.n_procs
def cli(
    bids_directory: str,
    caps_directory: str,
    uncropped_image: bool = False,
    random_seed: Optional[int] = None,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    use_antspy: bool = False,
    caps_name: Optional[str] = None,
) -> None:
    """Affine registration of Flair images to the MNI standard space.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Linear/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .anat_linear_pipeline import AnatLinear

    parameters = {
        "uncropped_image": uncropped_image,
        "random_seed": random_seed,
    }

    # Most of the time, you will want to instantiate your pipeline with a
    # BIDS and CAPS directory as inputs:
    pipeline = AnatLinear(
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
        use_antspy=use_antspy,
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
