from optparse import Option
from typing import List, Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_volume.create_dartel.pipeline as pydra_create_dartel
from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "pydra-create-dartel"


@click.command(name=pipeline_name, hidden=True)
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
    """Pydra Inter-subject registration using Dartel (creating a new Dartel template).

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects. For this pipeline, it is associated to the DARTEL template that you had created when running the t1-volume pipeline.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """

    parameters = {"group_label": group_label, "dartel_tissues": dartel_tissues}
    t1_volume_create_dartel_pipeline = pydra_create_dartel.t1volume_create_dartel(
        name="t1-volume-tissue-segmentation-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(t1_volume_create_dartel_pipeline)


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
