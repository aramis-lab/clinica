from typing import Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_volume.parcellation.pipeline as pydra_parcellation
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-t1-volume-parcellation"


@clinica_pipeline
@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.modulate
def cli(
    caps_directory: str,
    group_label: str,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    modulate: bool = True,
) -> None:
    """Computation of mean GM concentration for a set of regions (Pydra engine).

    GROUP_LABEL is an user-defined identifier to target a specific group of subjects.
    For this pipeline, it is associated to the DARTEL template that you had created
    when running the t1-volume pipeline.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """
    parameters = {
        "group_label": group_label,
        "modulate": modulate,
    }
    pipeline = pydra_parcellation.build_core_workflow(
        name=pipeline_name,
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
    )
    pydra_utils.run(pipeline, n_procs=n_procs)


if __name__ == "__main__":
    cli()
