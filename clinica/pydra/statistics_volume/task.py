import typing as ty

from pydra.mark import annotate, task


@task
@annotate({"return": {"idx_group1": list, "idx_group2": list, "class_names": list}})
def get_group_1_and_2_task(tsv: str, contrast: str) -> dict:
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
def write_matlab_model_task(
    tsv: str,
    contrast: str,
    idx_group1: list,
    idx_group2: list,
    file_list: list,
    template_file: str,
):
    """Pydra task for creating the .m file for the instantiation of the 2-sample t-test model in SPM

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.model_creation`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        write_matlab_model,
    )

    return write_matlab_model(
        tsv, contrast, idx_group1, idx_group2, file_list, template_file
    )


@task
@annotate({"return": {"output_mat_file": str}})
def run_m_script_task(m_file: str) -> str:
    """Pydra task for running an .m script

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.run_m_script`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    return run_m_script(m_file)


@task
@annotate({"return": {"current_model_estimation": str}})
def clean_template_file_task(mat_file: str, template_file: str) -> str:
    """Pydra task for estimation script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.estimate`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        clean_template_file,
    )

    return clean_template_file(mat_file, template_file)


@task
@annotate({"return": {"current_model_result": str}})
def clean_spm_result_file_task(
    mat_file: str, template_file: str, method: str, threshold: float
) -> str:
    """Pydra task for result script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.results`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        clean_spm_result_file,
    )

    return clean_spm_result_file(mat_file, template_file, method, threshold)


@task
@annotate({"return": {"current_model_estimation": str}})
def clean_spm_contrast_file_task(
    mat_file: str, template_file: str, covariates: list, class_names: list
):
    """Pydra task for contrast script.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.contrast`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        clean_spm_contrast_file,
    )

    return clean_spm_contrast_file(mat_file, template_file, covariates, class_names)


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
def copy_and_rename_spm_output_files_task(
    spm_mat: str,
    class_names: list,
    covariates: list,
    group_label: str,
    fwhm: int,
    measure: str,
):
    """Pydra task for reading outputs of the statistics volume pipeline.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume.statistics_volume_utils.read_output`.
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        copy_and_rename_spm_output_files,
    )

    return copy_and_rename_spm_output_files(
        spm_mat, class_names, covariates, group_label, fwhm, measure
    )
