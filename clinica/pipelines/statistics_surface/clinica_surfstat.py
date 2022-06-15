import numpy as np
from os import PathLike
from pathlib import Path
from typing import Dict

from ._inputs import (
    _read_and_check_tsv_file,
    _get_t1_freesurfer_custom_file_template,
    _build_thickness_array,
    _get_average_surface,
)
from ._model import GLMFactory

DEFAULT_FWHM = 20


def clinica_surfstat(
    input_dir: PathLike,
    output_dir: PathLike,
    tsv_file: PathLike,
    design_matrix: str,
    contrast: str,
    glm_type: str,
    group_label: str,
    freesurfer_home: PathLike,
    surface_file: PathLike,
    feature_label: str,
    parameters: Dict,
):
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
    input_dir : Input folder.

    output_dir : Output folder for storing results.

    tsv_file : Path to the TSV file `subjects.tsv` which contains the
        necessary metadata to run the statistical analysis.

        .. warning::
            The column names need to be accurate because they
            are used to defined contrast and model terms.
            Please double check for typos.

    design_matrix : Design matrix in string format.
        For example "1+Label"

    contrast : The contrast to be used in the GLM.

        .. warning::
            The contrast needs to be in the design matrix.

    glm_type : {"group_comparison", "correlation"}
        Type of GLM to run:
            - "group_comparison": Performs group comparison.
              For example "AD - ND".
            - "correlation": Performs correlation analysis.

    group_label : Label for the group.
        This is used in the output file names (see main description
        of the function).

    freesurfer_home : Path to the home folder of Freesurfer.
        This is required to get the fsaverage templates.

    surface_file : Path to the surface file to analyze.
        Typically the cortical thickness.

    feature_label : Label used for the measure.
        This is used in the output file names (see main description
        of the function).

    parameters : Dictionary of additional parameters
        - "sizeoffwhm": Smoothing. This is used in the output file names.
          Default=20.
        - "thresholduncorrectedpvalue": Threshold to be used with uncorrected
          P-values. Default=0.001.
        - "thresholdcorrectedpvalue": Threshold to be used with corrected
          P-values. Default=0.05.
        - "clusterthreshold": Threshold to be used to declare clusters as
          significant. Default=0.001.
    """
    # Load subjects data
    df_subjects = _read_and_check_tsv_file(tsv_file)
    if "sizeoffwhm" in parameters:
        fwhm = parameters["sizeoffwhm"]
    else:
        fwhm = DEFAULT_FWHM
        parameters["sizeoffwhm"] = fwhm
    surface_file = _get_t1_freesurfer_custom_file_template(input_dir)
    thickness = _build_thickness_array(input_dir, surface_file, df_subjects, fwhm)

    # Load average surface template
    fsaverage_path = freesurfer_home / Path("subjects/fsaverage/surf")
    average_surface, average_mesh = _get_average_surface(fsaverage_path)

    # Build and run GLM model
    glm_factory = GLMFactory(feature_label)
    glm_model = glm_factory.create_model(
        glm_type,
        design_matrix,
        df_subjects,
        contrast,
        group_label=group_label,
        **parameters,
    )
    glm_model.fit(thickness, average_surface)
    glm_model.save_results(output_dir, ["json", "mat"])
    glm_model.plot_results(output_dir, ["nilearn_plot_surf_stat_map"], average_mesh)
