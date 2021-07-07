from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "machinelearning-prepare-spatial-svm"


@click.command(pipeline_name)
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.orig_input_data_ml
@cli_param.option_group.pipeline_specific_options
@cli_param.option.acq_label
@cli_param.option.suvr_reference_region
@cli_param.option.use_pvc_data
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
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

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/MachineLearning_PrepareSVM/"
    """
    from networkx import Graph

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.ux import print_end_pipeline

    from .spatial_svm_pipeline import SpatialSVM

    if orig_input_data_ml == "pet-volume":
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

    pipeline = SpatialSVM(
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
