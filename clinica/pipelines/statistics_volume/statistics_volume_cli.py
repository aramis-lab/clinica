from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "statistics-volume"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.orig_input_data
@cli_param.argument.subject_visits_with_covariates_tsv
@cli_param.argument.contrast
@click.option(
    "-dartel",
    "--group_label_dartel",
    default="*",
    show_default=True,
    help="Name of the DARTEL template that Clinica needs to use to grab the input files.",
)
@click.option(
    "-fwhm",
    "--full_width_at_half_maximum",
    default=8,
    show_default=True,
    help="Full Width at Half Maximum (FWHM) of the smoothing used in your input files.",
)
@cli_param.option.acq_label
@cli_param.option.suvr_reference_region
@cli_param.option.use_pvc_data
@click.option(
    "-cf",
    "--custom_file",
    help=(
        "Custom file string. Specify filename using * when the subject or session name "
        "appears e.g. '*_task-rest_acq-fdg_pet_space-Ixi549Space_pet.nii.gz' will grab "
        "the corresponding file in all the subjects/sessions. "
        "This flag must be specified with the --measure_label flag). "
        "See Wiki for an example."
    ),
)
@click.option(
    "-ml",
    "--measure_label",
    help=(
        "Name of the feature type, it will be saved on the CAPS "
        "_measure-MEASURE_LABEL key-value association. "
        "This flag must be specified with the --custom_file flag). "
        "See Wiki for an example."
    ),
)
@click.option(
    "-ct",
    "--cluster_threshold",
    default=0.001,
    show_default=True,
    help="Threshold to define a cluster in the process of cluster-wise correction.",
)
@cli_param.option.working_directory
@cli_param.option.n_procs
def cli(
    caps_directory: str,
    group_label: str,
    orig_input_data: str,
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
    """Volume-based mass-univariate analysis with SPM."

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Volume/
    """
    from networkx import Graph

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.ux import print_end_pipeline

    from .statistics_volume_pipeline import StatisticsVolume

    # PET-Volume pipeline
    if orig_input_data == "pet-volume":
        if acq_label is None:
            raise ClinicaException(
                "You selected pet-volume pipeline without setting --acq_label flag. "
                "Clinica will now exit."
            )
        if suvr_reference_region is None:
            raise ClinicaException(
                "You selected pet-volume pipeline without setting --suvr_reference_region flag. "
                "Clinica will now exit."
            )

    # Custom pipeline
    if orig_input_data == "custom-pipeline":
        if (custom_file is None) or (measure_label is None):
            raise ClinicaException(
                "You must set --measure_label and --custom_file flags."
            )

    parameters = {
        # Clinica compulsory arguments
        "group_label": group_label,
        "orig_input_data": orig_input_data,
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
    }

    pipeline = StatisticsVolume(
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
