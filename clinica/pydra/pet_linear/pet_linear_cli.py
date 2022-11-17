from typing import List, Optional, Tuple

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.pet_linear.pipeline as pydra_pet_linear
from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "pydra-pet-linear"


@click.command(name=pipeline_name, hidden=False)
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
    """Affine registration of PET images to the MNI standard space (Pydra engine).

    Prerequisite: You need to have performed the t1-linear pipeline on your T1-weighted MR images.
    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Linear/"

    Parameters
    ----------
    bids_directory : str
        Path to the input folder containing a BIDS compliant dataset.

    caps_directory : str
        Path to a CAPS folder containing the results of the required pipelines (T1-linear).
        This will also be used for writing the outputs.

    acq_label : str
        Label given to the PET acquisition, specifying the tracer used.
        Frequently used values are '18FFDG' or '18FAV45'.

    suvr_reference_region : str
        ???
        The reference region must be specified to perform intensity normalization.
        Accepted values include: 'pons', 'cerebellumPons', 'pons2', 'cerebellumPons2'.

    uncropped_image : bool, optional
        Whether to also return the uncropped image or not. Default=False.

    save_pet_in_t1w_space : bool, optional
        ???. Default=False.

    subjects_sessions_tsv : str, optional
        ???. Default=None.

    working_directory : str, optional
        ???. Default=None.

    n_procs : int, optional
        Number of processes to use. Default=None.
    """
    parameters = {
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "uncropped_image": uncropped_image,
        "save_PETinT1w": save_pet_in_t1w_space,
    }
    pipeline = pydra_pet_linear.build_core_workflow(
        name="pet-linear-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        # tsv_file=subjects_sessions_tsv,
        # base_dir=working_directory,
        parameters=parameters,
    )
    pydra_utils.run(pipeline)


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
