from enum import Enum
from typing import Optional

import pandas as pd

from ._base import GLM
from ._correlation import CorrelationGLM
from ._group import GroupGLM, GroupGLMWithInteraction

__all__ = ["GLMModelType", "create_glm_model"]


class GLMModelType(str, Enum):
    """Supported types of GLM models."""

    CORRELATION = "correlation"
    GROUP_COMPARISON = "group_comparison"


def create_glm_model(
    glm_type: GLMModelType,
    design: str,
    df: pd.DataFrame,
    contrast: str,
    feature_label: str,
    group_label: Optional[str] = "group",
    fwhm: Optional[int] = 20,
    threshold_uncorrected_pvalue: Optional[float] = 0.001,
    threshold_corrected_pvalue: Optional[float] = 0.05,
    cluster_threshold: Optional[float] = 0.001,
) -> GLM:
    """Factory method for building a GLM model instance corresponding to the
    provided type and design matrix.

    Parameters
    ----------
    glm_type : GLMModelType
        The type of GLM to be created.

    design : str
        The design matrix specified in string format.
        If this contains a "*", it will be interpreted as an interaction effect.

    df : pd.DataFrame
        The subjects DataFrame.

    contrast : str
        The contrast specified in string format.

    feature_label : str
        The label used for building output filenames.

    group_label : str, optional
        The label to use for group GLM models. Default="group".

    fwhm : int, optional
        The smoothing FWHM. This is used in the output file names.
        Default=20.

    threshold_uncorrected_pvalue : float, optional
        The threshold to be used with uncorrected P-values. Default=0.001.

    threshold_corrected_pvalue : float, optional
        The threshold to be used with corrected P-values. Default=0.05.

    cluster_threshold : float, optional
        The threshold to be used to declare clusters as significant. Default=0.001.

    Returns
    -------
    model : GLM
        An instance of the `GLM` class.

    Raises
    ------
    ValueError
        If the glm_type is not supported.
    """
    from clinica.utils.stream import cprint

    cprint(
        msg=f"The GLM model is: {design} and the GLM type is: {glm_type}",
        lvl="info",
    )
    params = {
        "group_label": group_label,
        "fwhm": fwhm,
        "threshold_uncorrected_pvalue": threshold_uncorrected_pvalue,
        "threshold_corrected_pvalue": threshold_corrected_pvalue,
        "cluster_threshold": cluster_threshold,
    }
    if glm_type == GLMModelType.CORRELATION:
        return CorrelationGLM(design, df, feature_label, contrast, **params)
    if glm_type == GLMModelType.GROUP_COMPARISON:
        if "*" in design:
            return GroupGLMWithInteraction(
                design, df, feature_label, contrast, **params
            )
        return GroupGLM(design, df, feature_label, contrast, **params)
