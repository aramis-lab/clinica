"""This file contains utilities for building contrasts for
different types of GLMs.
"""


import numpy as np
import pandas as pd
from string import Template
from typing import Dict, Tuple

from ._inputs import _check_contrast
from ._utils import _is_categorical


def _get_contrasts_and_filenames(
        glm_type: str,
        contrast: str,
        df: pd.DataFrame,
):
    """Transforms the contrast in string format into matrix format.

    Parameters
    ----------
    glm_type : Type of GLM to run (either 'group_comparision' or 'correlation')
    contrast : Desired contrast.
    df : Subject DataFrame.

    Returns
    -------
    contrasts : Dictionary <contrast_name>:<contrast_vector>
    filenames : Dictionary <contrast_name>:<root_filename>
    """
    (
        abs_contrast, contrast_sign, with_interaction
    ) = _check_contrast(contrast, df, glm_type)
    if glm_type == "group_comparison":
        if not with_interaction:
            return _get_group_contrast_without_interaction(abs_contrast, df)
        else:
            return _get_group_contrast_with_interaction(abs_contrast, df)
    elif glm_type == "correlation":
        return _get_correlation_contrast(abs_contrast, df, contrast_sign)
    else:
        raise ValueError(
            "Check out if you define the glmtype flag correctly, "
            "or define your own general linear model, e,g MGLM."
        )


def _get_group_contrast_with_interaction(
        contrast: str,
        df: pd.DataFrame,
) -> Tuple[Dict, Dict]:
    """Build contrasts and filename roots for group GLMs with interaction.

    Parameters
    ----------
    contrast : Contrast in string format.
    df : Subjects DataFrame.

    Returns
    -------
    contrasts : Dictionary <contrast_name>:<contrast_vector>
    filenames : Dictionary <contrast_name>:<root_filename>
    """
    contrasts = dict()
    filenames = dict()
    contrast_elements = [_.strip() for _ in contrast.split("*")]
    categorical = [_is_categorical(df, _) for _ in contrast_elements]
    if len(contrast_elements) != 2 or sum(categorical) != 1:
        raise ValueError(
            "The contrast must be an interaction between one continuous "
            "variable and one categorical variable. Your contrast contains "
            f"the following variables : {contrast_elements}"
        )
    idx = 0 if categorical[0] else 1
    categorical_contrast = contrast_elements[idx]
    continue_contrast = contrast_elements[(idx + 1) % 2]
    group_values = np.unique(df[categorical_contrast])
    built_contrast = df[continue_contrast] * (
        (df[categorical_contrast] == group_values[0]).astype(int)
    ) - df[continue_contrast] * (
        (df[categorical_contrast] == group_values[1]).astype(int)
    )
    contrasts[contrast] = built_contrast
    filenames[contrast] = (
        Template("interaction-${contrast_name}_measure-${feature_label}_fwhm-${fwhm}")
    )
    return contrasts, filenames


def _get_group_contrast_without_interaction(
        contrast: str,
        df: pd.DataFrame,
):
    """Build contrasts and filename roots for group GLMs without interaction.

    Parameters
    ----------
    contrast : Contrast in string format.
    df : Subjects DataFrame.

    Returns
    -------
    contrasts : Dictionary <contrast_name>:<contrast_vector>
    filenames : Dictionary <contrast_name>:<root_filename>
    """
    contrasts = dict()
    filenames = dict()
    if not _is_categorical(df, contrast):
        raise ValueError(
            "Contrast should refer to a categorical variable for group comparison. "
            "Please select 'correlation' for 'glm_type' otherwise."
        )
    group_values = np.unique(df[contrast])
    for contrast_type, (i, j) in zip(["positive", "negative"], [(0, 1), (1, 0)]):
        contrast_name = f"{group_values[i]}-lt-{group_values[j]}"
        contrasts[contrast_name] = (
            (df[contrast] == group_values[i]).astype(int) -
            (df[contrast] == group_values[j]).astype(int)
        )
        filenames[contrast_name] = (
            Template("group-${group_label}_${contrast_name}_measure-${feature_label}_fwhm-${fwhm}")
        )
    return contrasts, filenames


def _get_correlation_contrast(
        contrast: str,
        df: pd.DataFrame,
        contrast_sign: str,
):
    """Build contrasts and filename roots for correlation GLMs.

    Parameters
    ----------
    contrast : Contrast in string format.
    df : Subjects DataFrame.
    contrast_sign : Sign of the contrast (either 'positive' or 'negative').

    Returns
    -------
    contrasts : Dictionary <contrast_name>:<contrast_vector>
    filenames : Dictionary <contrast_name>:<root_filename>
    """
    contrasts = dict()
    filenames = dict()
    built_contrast = df[contrast]
    if contrast_sign == "negative":
        built_contrast *= -1
    contrasts[contrast] = built_contrast
    filenames[contrast] = Template(
        "group-${group_label}_correlation-${contrast_name}-"
        f"{contrast_sign}_"
        "measure-${feature_label}_fwhm-${fwhm}"
    )
    return contrasts, filenames

