from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "pet-surface-longitudinal"


@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.argument.pvc_psf_tsv
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
def cli(
    bids_directory: str,
    caps_directory: str,
    acq_label: str,
    suvr_reference_region: str,
    pvc_psf_tsv: str,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Longitudinal surface-based processing of PET images.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Surface_Longitudinal/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .pet_surface_pipeline import PetSurface

    parameters = {
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "pvc_psf_tsv": pvc_psf_tsv,
        "longitudinal": True,
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
