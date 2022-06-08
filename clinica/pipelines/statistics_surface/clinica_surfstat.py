from os import PathLike
from pathlib import Path
from typing import Dict

import numpy as np
from brainstat.stats.SLM import SLM
from clinica.utils.stream import cprint

from ._contrasts import _get_contrasts_and_filenames
from ._inputs import (
    _build_thickness_array,
    _extract_parameters,
    _get_average_surface,
    _get_t1_freesurfer_custom_file_template,
    _read_and_check_tsv_file,
)
from ._model import _build_model
from ._outputs import _plot_results, _print_clusters, _save_results


def _compute_results(
    model: SLM,
    mask: np.ndarray,
    threshold_uncorrected_pvalue: float,
    threshold_corrected_pvalue: float,
) -> Dict:
    """Take a fitted SLM model and store the results in a dictionary.

    Parameters
    ----------
    model : Fitted SLM model.
    mask : Mask used in the analysis.
    threshold_uncorrected_pvalue : Threshold used with non corrected P-values.
    threshold_corrected_pvalue : Threshold used with corrected P-values.

    Returns
    -------
    results : Dictionary with the results.
    """
    from scipy.stats import t

    results = dict()
    results["coefficients"] = np.nan_to_num(model.coef)
    results["TStatistics"] = np.nan_to_num(model.t)
    results["uncorrectedPValue"] = dict()
    results["uncorrectedPValue"]["P"] = 1 - t.cdf(results["TStatistics"], model.df)
    results["uncorrectedPValue"]["mask"] = mask
    results["uncorrectedPValue"]["thresh"] = threshold_uncorrected_pvalue
    results["FDR"] = model._fdr()
    results["correctedPValue"] = dict()
    results["correctedPValue"]["P"] = model.P["pval"]["P"]
    results["correctedPValue"]["C"] = model.P["pval"]["C"]
    results["correctedPValue"]["mask"] = mask
    results["correctedPValue"]["thresh"] = threshold_corrected_pvalue
    return results


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

    design_matrix : Design matrix in string format. For example "1+Label"
    contrast : The contrast to be used in the GLM.

        .. warning::
            The contrast needs to be in the design matrix.

    glm_type : {"group_comparison", "correlation"}
        Type of GLM to run:
            - "group_comparison": Performs group comparison.
              For example "AD - ND".
            - "correlation": Performs correlation analysis.

    group_label : Label for the group.
    freesurfer_home : Path to the home folder of Freesurfer.
        This is required to get the fsaverage templates.
    surface_file : Path to the surface file.
    """
    (
        fwhm,
        threshold_uncorrected_pvalue,
        threshold_corrected_pvalue,
        cluster_threshold,
    ) = _extract_parameters(parameters)
    fsaverage_path = freesurfer_home / Path("subjects/fsaverage/surf")
    cprint(
        msg=f"fsaverage path : {fsaverage_path}",
        lvl="info",
    )
    df_subjects = _read_and_check_tsv_file(tsv_file)
    surface_file = _get_t1_freesurfer_custom_file_template(input_dir)
    thickness = _build_thickness_array(input_dir, surface_file, df_subjects, fwhm)
    mask = thickness[0, :] > 0
    average_surface, average_mesh = _get_average_surface(fsaverage_path)
    cprint(
        msg=f"The GLM model is: {design_matrix} and the GLM type is: {glm_type}",
        lvl="info",
    )
    contrasts, filenames = _get_contrasts_and_filenames(glm_type, contrast, df_subjects)
    naming_parameters = {
        "fwhm": fwhm,
        "group_label": group_label,
        "feature_label": feature_label,
    }
    model = _build_model(design_matrix, df_subjects)
    for contrast_name, model_contrast in contrasts.items():
        filename_root = Path(output_dir) / Path(
            filenames[contrast_name].safe_substitute(
                contrast_name=contrast_name, **naming_parameters
            )
        )
        slm_model = SLM(
            model,
            contrast=model_contrast,
            surf=average_surface,
            mask=mask,
            two_tailed=True,
            correction=["fdr", "rft"],
            cluster_threshold=cluster_threshold,
        )
        cprint(
            msg=f"Fitting the SLM model with contrast {contrast_name}...",
            lvl="info",
        )
        slm_model.fit(thickness)
        results = _compute_results(
            slm_model, mask, threshold_uncorrected_pvalue, threshold_corrected_pvalue
        )
        _save_results(results, filename_root, out_formats="all")
        try:
            _plot_results(results, filename_root, average_mesh)
        except:  # noqa
            cprint(
                msg="Plotting failed...",
                lvl="error",
            )
            pass
        _print_clusters(slm_model, threshold_corrected_pvalue)
