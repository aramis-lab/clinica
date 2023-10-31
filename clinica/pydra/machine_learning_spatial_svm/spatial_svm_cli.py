from typing import Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.machine_learning_spatial_svm.pipeline as pydra_machine_learning_spatial_svm
from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-machine-learning-prepare-spatial-svm"


@clinica_pipeline
@click.command(pipeline_name)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.orig_input_data_ml
@cli_param.option_group.custom_pipeline_options(
    "Pipeline options if you use inputs from pet-volume pipeline"
)
@cli_param.option.acq_label
@cli_param.option.suvr_reference_region
@cli_param.option.use_pvc_data
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@option.global_option_group
@option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option_group.option(
    "-fwhm",
    "--full_width_half_maximum",
    default=4.0,
    help=(
        "Amount of regularization (in mm). In practice, we found the default value "
        "(--full_width_half_maximum %(default)s) to be optimal. We therefore "
        "do not recommend to change it unless you have a specific reason to do so."
    ),
)
def cli(
    caps_directory: str,
    group_label: str,
    orig_input_data_ml: str,
    acq_label: Optional[str] = None,
    suvr_reference_region: Optional[str] = None,
    use_pvc_data: bool = False,
    full_width_half_maximum: float = 4.0,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Prepare input data for SVM with spatial and anatomical regularization.

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects.

       The third positional argument can be t1-volume to use tissue maps or pet-volume to use standardized uptake value ratio (SUVR) maps.

    Prerequisites: You need to execute the t1-volume pipeline to run the pipeline on T1-weighted MRI data, and the t1-volume + pet-volume pipelines to apply it to PET data.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/MachineLearning_PrepareSVM/"
    """

    from clinica.utils.exceptions import ClinicaException

    parameters = {
        # Clinica compulsory arguments
        "group_label": group_label,
        "orig_input_data_ml": orig_input_data_ml,
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": acq_label,
        "use_pvc_data": use_pvc_data,
        "suvr_reference_region": suvr_reference_region,
        # Advanced arguments
        "fwhm": full_width_half_maximum,
    }

    pipeline = pydra_machine_learning_spatial_svm.build_core_workflow(
        name="machine-learning-spatial-svm-pydra",
        input_dir=None,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(pipeline, n_procs)


if __name__ == "__main__":
    cli()
