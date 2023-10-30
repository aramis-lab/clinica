from typing import Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.statistics_volume.pipeline as pydra_statistics_volume
from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-statistics-volume"


@clinica_pipeline
@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.orig_input_data_volume
@cli_param.argument.subject_visits_with_covariates_tsv
@cli_param.argument.contrast
@cli_param.option_group.pipeline_specific_options
@cli_param.option_group.option(
    "-dartel",
    "--group_label_dartel",
    default="*",
    show_default=True,
    help="Name of the DARTEL template that Clinica needs to use to grab the input files.",
)
@cli_param.option_group.option(
    "-fwhm",
    "--full_width_at_half_maximum",
    default=8,
    show_default=True,
    help="Full Width at Half Maximum (FWHM) of the smoothing used in your input files.",
)
@cli_param.option_group.custom_pipeline_options(
    "Pipeline options if you use inputs from pet-volume pipeline"
)
@cli_param.option.acq_label
@cli_param.option.suvr_reference_region
@cli_param.option.use_pvc_data
@cli_param.option_group.custom_pipeline_options(
    "Pipeline options if you selected custom-pipeline"
)
@cli_param.option_group.option(
    "-cf",
    "--custom_file",
    help=(
        "Custom file string. Specify filename using * when the subject or session name "
        "appears e.g. '*_trc-18FFDG_pet_space-Ixi549Space_pet.nii.gz' will grab "
        "the corresponding file in all the subjects/sessions. "
        "This flag must be specified with the --measure_label flag). "
        "See Wiki for an example."
    ),
)
@cli_param.option_group.option(
    "-ml",
    "--measure_label",
    help=(
        "Name of the feature type, it will be saved on the CAPS "
        "_measure-MEASURE_LABEL key-value association. "
        "This flag must be specified with the --custom_file flag). "
        "See Wiki for an example."
    ),
)
@cli_param.option_group.common_pipelines_options
@cli_param.option.working_directory
@option.global_option_group
@option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option_group.option(
    "-ct",
    "--cluster_threshold",
    default=0.001,
    show_default=True,
    help="Threshold to define a cluster in the process of cluster-wise correction.",
)
def cli(
    caps_directory: str,
    group_label: str,
    orig_input_data_volume: str,
    subject_visits_with_covariates_tsv: str,
    contrast: str,
    group_label_dartel: str = "*",
    full_width_at_half_maximum: int = 8,
    acq_label: Optional[str] = None,
    suvr_reference_region: Optional[str] = None,
    use_pvc_data: bool = False,
    custom_file: Optional[str] = None,
    measure_label: Optional[str] = None,
    cluster_threshold: float = 1e-3,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Volume-based mass-univariate analysis with SPM.

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects.

       The type of volume-based feature can be defined by using the third argument: t1-volume to use gray matter maps, pet-volume to use PET data or custom-pipeline for you own data in CAPS directory.

       SUBJECT_VISITS_WITH_COVARIATES_TSV is a TSV file containing a list of subjects with their sessions and all the covariates and factors of the model.

       CONTRAST is a string defining the contrast matrix or the variable of interest for the GLM.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Volume/
    """
    from clinica.utils.exceptions import ClinicaException

    parameters = {
        # Clinica compulsory arguments
        "group_label": group_label,
        "orig_input_data_volume": orig_input_data_volume,
        "contrast": contrast,
        # Optional arguments
        "group_label_dartel": group_label_dartel,
        "full_width_at_half_maximum": full_width_at_half_maximum,
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": acq_label,
        "use_pvc_data": use_pvc_data,
        "suvr_reference_region": suvr_reference_region,
        # Optional arguments for custom pipeline
        "measure_label": measure_label,
        "custom_file": custom_file,
        # Advanced arguments
        "cluster_threshold": cluster_threshold,
        "tsv_file": subject_visits_with_covariates_tsv,
        "base_dir": working_directory,
    }

    statistics_volume_pipeline = pydra_statistics_volume.build_core_workflow(
        name="pydra-statistics-volume",
        input_dir=None,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(statistics_volume_pipeline, n_procs=n_procs)


if __name__ == "__main__":
    cli()
