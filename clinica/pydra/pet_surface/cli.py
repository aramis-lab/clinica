from typing import Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

NAME = "pydra-pet-surface"


@clinica_pipeline
@click.command(name=NAME)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.argument.pvc_psf_tsv
@cli_param.option_group.common_pipelines_options
@cli_param.option.reconstruction_method
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option.yes
def cli(
    bids_directory: str,
    caps_directory: str,
    acq_label: str,
    suvr_reference_region: str,
    pvc_psf_tsv: str,
    reconstruction_method: Optional[str] = None,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    yes: bool = False,
) -> None:
    ...
