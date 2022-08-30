from os import PathLike
from pathlib import Path
from typing import Dict, Optional

from ._inputs import (
    _build_thickness_array,
    _get_average_surface,
    _get_t1_freesurfer_custom_file_template,
    _read_and_check_tsv_file,
)
from ._model import create_glm_model


def clinica_surfstat(
    input_dir: PathLike,
    output_dir: PathLike,
    tsv_file: PathLike,
    design_matrix: str,
    contrast: str,
    glm_type: str,
    group_label: str,
    freesurfer_home: PathLike,
    surface_file: Optional[PathLike],
    feature_label: str,
    fwhm: Optional[int] = 20,
    threshold_uncorrected_pvalue: Optional[float] = 0.001,
    threshold_corrected_pvalue: Optional[float] = 0.05,
    cluster_threshold: Optional[float] = 0.001,
) -> None:
    """This function mimics the previous function `clinica_surfstat`
    written in MATLAB and relying on the MATLAB package SurfStat.

    This implementation is written in pure Python and rely on the
    package brainstat for GLM modeling.

    Results, both plots and matrices, will be written in the location
    specified through `output_dir`.

    The names of the output files within `output_dir` follow the
    conventions:

        <ROOT>_<SUFFIX>.<EXTENSION>

    EXTENSION can be:

        - "mat": for storing matrices (this is mainly for backward
          compatibility with the previous MATLAB implementation of
          this function).
        - "json": for storing matrices in a more Pythonic way.
        - "png" for surface figures.

    SUFFIX can be:

        - "coefficients": relative to the model's beta coefficients.
        - "TStatistics": relative to the T-statistics.
        - "uncorrectedPValue": relative to the uncorrected P-values.
        - "correctedPValues": relative to the corrected P-values.
        - "FDR": Relative to the False Discovery Rate.

    ROOT can be:

        - For group comparison GLM with an interaction term:

            interaction-<contrast_name>_measure-<feature_label>_fwhm-<fwhm>

        - For group comparison GLM without an interaction term:

            group-<group_label>_<contrast_name>_measure-<feature_label>_fwhm-<fwhm>

        - For correlation GLM:

            group-<group_label>_correlation-<contrast_name>-<contrast_sign>_measure-<feature_label>_fwhm-<fwhm>

    Parameters
    ----------
    input_dir : PathLike
        Path to the input folder.

    output_dir : PathLike
        Path to the output folder for storing results.

    tsv_file : PathLike
        Path to the TSV file `subjects.tsv` which contains the
        necessary metadata to run the statistical analysis.

        .. warning::
            The column names need to be accurate because they
            are used to defined contrast and model terms.
            Please double check for typos.

    design_matrix : str
        The design matrix specified in string format.
        For example "1 + Label"

    contrast : str
        The contrast to be used in the GLM, specified in string format.

        .. warning::
            The contrast needs to be in the design matrix.

    glm_type : {"group_comparison", "correlation"}
        Type of GLM to run:
            - "group_comparison": Performs group comparison.
              For example "AD - ND".
            - "correlation": Performs correlation analysis.

    group_label : str
        The label for the group. This is used in the output file names
        (see main description of the function).

    freesurfer_home : PathLike
        The path to the home folder of Freesurfer.
        This is required to get the fsaverage templates.

    surface_file : PathLike, optional
        The path to the surface file to analyze.
        Typically the cortical thickness.
        If `None`, the surface file will be the t1 freesurfer template.

    feature_label : str
        The label used for the measure. This is used in the output file
        names (see main description of the function).

    fwhm : int, optional
        The smoothing FWHM. This is used in the output file names.
        Default=20.

    threshold_uncorrected_pvalue : float, optional
        The threshold to be used with uncorrected P-values. Default=0.001.

    threshold_corrected_pvalue : float, optional
        The threshold to be used with corrected P-values. Default=0.05.

    cluster_threshold : float, optional
        The threshold to be used to declare clusters as significant. Default=0.05.
    """
    # Load subjects data
    df_subjects = _read_and_check_tsv_file(tsv_file)
    if surface_file is None:
        surface_file = _get_t1_freesurfer_custom_file_template(input_dir)
    thickness = _build_thickness_array(input_dir, surface_file, df_subjects, fwhm)

    # Load average surface template
    fsaverage_path = freesurfer_home / Path("subjects/fsaverage/surf")
    average_surface, average_mesh = _get_average_surface(fsaverage_path)

    # Build and run GLM model
    glm_model = create_glm_model(
        glm_type,
        design_matrix,
        df_subjects,
        contrast,
        feature_label,
        group_label=group_label,
        fwhm=fwhm,
        threshold_uncorrected_pvalue=threshold_uncorrected_pvalue,
        threshold_corrected_pvalue=threshold_corrected_pvalue,
        cluster_threshold=cluster_threshold,
    )
    glm_model.fit(thickness, average_surface)
    glm_model.save_results(output_dir, ["json", "mat"])
    glm_model.plot_results(output_dir, ["nilearn_plot_surf_stat_map"], average_mesh)
