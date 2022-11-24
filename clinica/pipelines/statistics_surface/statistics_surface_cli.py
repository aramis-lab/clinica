from typing import List, Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "statistics-surface"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@click.argument(
    "orig_input_data",
    type=click.Choice(["t1-freesurfer", "pet-surface", "custom-pipeline"]),
)
@click.argument(
    "glm_type",
    type=click.Choice(["group_comparison", "correlation"]),
)
@click.argument(
    "subject_visits_with_covariates_tsv",
    type=click.Path(exists=True, resolve_path=True),
)
@cli_param.argument.contrast
@cli_param.option_group.pipeline_specific_options
@cli_param.option_group.option(
    "-c",
    "--covariates",
    multiple=True,
    help=(
        "List of covariates. Each covariate must match the column name of the TSV file. "
        "By default, no covariate is taken."
    ),
)
@cli_param.option_group.option(
    "-fwhm",
    "--full_width_at_half_maximum",
    default=20,
    show_default=True,
    help="FWHM for the surface smoothing.",
)
@cli_param.option_group.custom_pipeline_options(
    "Pipeline options if you use inputs from pet-surface pipeline"
)
@cli_param.option.acq_label
@cli_param.option.suvr_reference_region
@cli_param.option_group.custom_pipeline_options(
    "Pipeline options if you selected custom-pipeline"
)
@cli_param.option_group.option(
    "-cf",
    "--custom_file",
    help=(
        "Pattern of file inside CAPS directory using @subject, @session, "
        "@fwhm, @hemi. This flag must be specified with the --measure_label flag). "
        "See Wiki for an example."
    ),
)
@cli_param.option_group.option(
    "-ml",
    "--measure_label",
    help=(
        "Name of the feature type, it will be saved on the CAPS "
        "_measure-FEATURE_LABEL key-value association. "
        "This flag must be specified with the --custom_file flag). "
        "See Wiki for an example."
    ),
)
@cli_param.option_group.common_pipelines_options
@cli_param.option.working_directory
@cli_param.option.n_procs
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
    orig_input_data: str,
    glm_type: str,
    subject_visits_with_covariates_tsv: str,
    contrast: str,
    covariates: Optional[List[str]] = None,
    full_width_at_half_maximum: int = 20,
    acq_label: Optional[str] = None,
    suvr_reference_region: Optional[str] = None,
    custom_file: Optional[str] = None,
    measure_label: Optional[str] = None,
    cluster_threshold: float = 1e-3,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Surface-based mass-univariate analysis with SurfStat.

    GROUP_LABEL is a user-defined identifier to target a specific group of
    subjects.

    The type of surface-based feature can be defined by using the third
    argument: t1-freesurfer for cortical thickness, pet-surface for projected
    PET data or custom-pipeline for you own data in CAPS directory.

    The type of analysis of the model is defined by the argument
    'group_comparison' or 'correlation'.

    SUBJECT_VISITS_WITH_COVARIATES_TSV is a TSV file containing a list of
    subjects with their sessions and all the covariates and factors in your
    model.

    CONTRAST is a string defining the contrast matrix or the variable of
    interest in the GLM, e.g. 'group' or 'age'

    Prerequisite: You need to have performed the t1-freesurfer pipeline on
    your T1-weighted MR images or pet-surface pipeline for measurements of
    activity map from PET.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Surface/
    """
    from networkx import Graph

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.ux import print_end_pipeline

    from .statistics_surface_pipeline import StatisticsSurface
    from .statistics_surface_utils import (
        get_pet_surface_custom_file,
        get_t1_freesurfer_custom_file,
    )

    # PET-Surface pipeline
    if orig_input_data == "pet-surface":
        if not acq_label:
            raise ClinicaException(
                "You selected pet-surface pipeline without setting --acq_label flag. "
                "Clinica will now exit."
            )
        if not suvr_reference_region:
            raise ClinicaException(
                "You selected pet-surface pipeline without setting --suvr_reference_region flag. "
                "Clinica will now exit."
            )

    # FreeSurfer cortical thickness
    if orig_input_data == "t1-freesurfer":
        custom_file = get_t1_freesurfer_custom_file()
        measure_label = "ct"
    # PET cortical projection
    elif orig_input_data == "pet-surface":
        custom_file = get_pet_surface_custom_file(acq_label, suvr_reference_region)
        measure_label = acq_label
    else:
        if not all([custom_file, measure_label]):
            raise ClinicaException(
                "You must set --measure_label and --custom_file flags."
            )

    parameters = {
        # Clinica compulsory arguments
        "group_label": group_label,
        "orig_input_data": orig_input_data,
        "glm_type": glm_type,
        "contrast": contrast,
        # Optional arguments
        "covariates": covariates,
        "full_width_at_half_maximum": full_width_at_half_maximum,
        # Optional arguments for inputs from pet-surface pipeline
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        # Optional arguments for custom pipeline
        "custom_file": custom_file,
        "measure_label": measure_label,
        # Advanced arguments (i.e. tricky parameters)
        "cluster_threshold": cluster_threshold,
    }

    pipeline = StatisticsSurface(
        caps_directory=caps_directory,
        tsv_file=subject_visits_with_covariates_tsv,
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
