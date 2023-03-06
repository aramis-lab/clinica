import typing as ty
from pathlib import PurePath

from pydra.mark import annotate, task


@task
@annotate(
    {"return": {"first_group_idx": list, "second_group_idx": list, "class_names": list}}
)
def get_group_1_and_2_task(tsv, contrast):
    """Pydra task for computing indexes of each group.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.get_group_1_and_2`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        get_group_1_and_2,
    )

    return get_group_1_and_2(tsv, contrast)


@task
@annotate({"return": {"current_model": str, "covariates": list}})
def model_creation_task(
    tsv, contrast, idx_group1, idx_group2, file_list, template_file
):
    """Pydra task for creating the .m file for the instantiation of the 2-sample t-test model in SPM

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.model_creation`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        model_creation,
    )

    return model_creation(
        tsv, contrast, idx_group1, idx_group2, file_list, template_file
    )


@task
@annotate({"return": {"output_mat_file": str}})
def run_m_script_task(m_file):
    """Pydra task for running an .m script

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.run_m_script`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    return run_m_script(m_file)


@task
@annotate({"return": {"current_model_estimation": str}})
def estimate_task(mat_file, template_file):
    """Pydra task for estimation script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.estimate`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import estimate

    return estimate(mat_file, template_file)


@task
@annotate({"return": {"current_model_result": str}})
def results_task(mat_file, template_file, method, threshold):
    """Pydra task for result script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.results`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import results

    return results(mat_file, template_file, method, threshold)


@task
@annotate({"return": {"current_model_estimation": str}})
def contrast_task(mat_file, template_file, covariates, class_names):
    """Pydra task for contrast script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.contrast`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import contrast

    return contrast(mat_file, template_file, covariates, class_names)


@task
@annotate(
    {
        "return": {
            "spmT_0001": str,
            "spmT_0002": str,
            "spm_figures": list,
            "variance_of_error": str,
            "resels_per_voxels": str,
            "mask": str,
            "regression_coeff": list,
            "contrasts": list[str],
        }
    }
)
def read_output_task(spm_mat, class_names, covariates, group_label, fwhm, measure):
    """Pydra task for reading outputs of the statistics volume pipeline.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.read_output`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import read_output

    return read_output(spm_mat, class_names, covariates, group_label, fwhm, measure)
