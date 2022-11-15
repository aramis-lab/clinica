from typing import Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "pet-linear"


@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.option_group.pipeline_specific_options
@cli_param.option_group.option(
    "-ui",
    "--uncropped_image",
    is_flag=True,
    help="Do not crop the image with template (cropped image are suggested for using with DL models)",
)
@cli_param.option_group.option(
    "--save_pet_in_t1w_space",
    is_flag=True,
    help="Save the PET image in the T1w space computed in the intermediate step of the pipeline",
)
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
def cli(
    bids_directory: str,
    caps_directory: str,
    acq_label: str,
    suvr_reference_region: str,
    uncropped_image: bool = False,
    save_pet_in_t1w_space: bool = False,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Affine registration of PET images to the MNI standard space.

       ACQ_LABEL corresponds the label given to the PET acquisition, specifying the tracer used.
    Frequently used values are '18FFDG' or '18FAV45'.

       The reference region must be specified to perform intensity normalization.
    Accepted values include: 'pons', 'cerebellumPons', 'pons2', 'cerebellumPons2'.

    Prerequisite: You need to have performed the t1-linear pipeline on your T1-weighted MR images.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Linear/"
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .pet_linear_pipeline import PETLinear

    parameters = {
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "uncropped_image": uncropped_image,
        "save_PETinT1w": save_pet_in_t1w_space,
    }

    pipeline = PETLinear(
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


run_cli.add_command(cli)


if __name__ == "__main__":
    cli()
