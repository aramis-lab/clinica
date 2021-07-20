from typing import List, Optional

import click

from clinica.pipelines import cli_param

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
@cli_param.option_group.pipeline_options
@cli_param.option_group.option(
    "-c",
    "--covariate",
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
@cli_param.option_group.standard_options
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option_group.advanced_options
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

    See "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Surface/
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
        if acq_label is None:
            raise ClinicaException(
                "You selected pet-surface pipeline without setting --acq_label flag. "
                "Clinica will now exit."
            )
        if suvr_reference_region is None:
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
        if (custom_file is None) or (measure_label is None):
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


if __name__ == "__main__":
    cli()
