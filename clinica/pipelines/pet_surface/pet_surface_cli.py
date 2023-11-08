from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pet-surface"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.argument.pvc_psf_tsv
@cli_param.option_group.common_pipelines_options
@cli_param.option.reconstruction_method
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.yes
@option.global_option_group
@option.n_procs
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
    """Surface-based processing of PET images.

       ACQ_LABEL corresponds to the label given to the PET acquisition, specifying the tracer used.
    Frequently used values are '18FFDG' or '18FAV45'.

       The reference region must be specified to perform intensity normalization.
    Accepted values include: 'pons', 'cerebellumPons', 'pons2', 'cerebellumPons2'.

       PVC_PSF_TSV is the TSV file containing the psf_x, psf_y and psf_z of the PSF for each PET image.

    Prerequisite: You need to have performed the t1-freesurfer pipeline on your T1-weighted MR images.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Surface/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .pet_surface_pipeline import PetSurface

    parameters = {
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "reconstruction_method": reconstruction_method,
        "pvc_psf_tsv": pvc_psf_tsv,
        "longitudinal": False,
        "skip_question": yes,
    }

    pipeline = PetSurface(
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
