import warnings
from typing import List, Optional

import pandas as pd

from ._base import GLM
from ._contrast import GroupContrast, GroupContrastWithInteraction


class GroupGLM(GLM):
    """Class implementing group GLM models.

    Attributes
    ----------
    See documentation for `GLM` class.

    group_label : str, optional
        The Label to use for group GLM models. Default="group".
    """

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        group_label: Optional[str] = "group",
        fwhm: Optional[int] = 20,
        threshold_uncorrected_pvalue: Optional[float] = 0.001,
        threshold_corrected_pvalue: Optional[float] = 0.05,
        cluster_threshold: Optional[float] = 0.001,
    ):
        super().__init__(
            design,
            df,
            feature_label,
            contrast,
            fwhm,
            threshold_uncorrected_pvalue,
            threshold_corrected_pvalue,
            cluster_threshold,
        )
        self.with_interaction = False
        self.group_label = group_label

    def _build_contrasts(
        self, contrast: str, subjects_df: pd.DataFrame
    ) -> List[GroupContrast]:
        return [
            GroupContrast.from_string(contrast, subjects_df, sign)
            for sign in ("positive", "negative")
        ]

    def _get_output_filename(self, contrast: GroupContrast) -> str:
        """Build the filename root part from class attributes and provided contrast.

        Parameters
        ----------
        contrast : GroupContrast
            The contrast to use for building the filename.
        """
        return f"group-{self.group_label}_{contrast.name}_measure-{self.feature_label}_fwhm-{self.fwhm}"


class GroupGLMWithInteraction(GroupGLM):
    """This class implements a GLM model for group comparison with
    interaction effects.

    Attributes
    ----------
    See attributes of parent class `GroupGLM`.
    """

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        group_label: Optional[str] = "group",
        fwhm: Optional[int] = 20,
        threshold_uncorrected_pvalue: Optional[float] = 0.001,
        threshold_corrected_pvalue: Optional[float] = 0.05,
        cluster_threshold: Optional[float] = 0.001,
    ):
        super().__init__(
            design,
            df,
            feature_label,
            contrast,
            group_label,
            fwhm,
            threshold_uncorrected_pvalue,
            threshold_corrected_pvalue,
            cluster_threshold,
        )
        self.with_interaction: bool = True
        warnings.warn(
            "You included interaction as covariate in your model, "
            "please carefully check the format of your tsv files."
        )

    def _build_contrasts(
        self, contrast: str, subjects_df: pd.DataFrame
    ) -> List[GroupContrastWithInteraction]:
        return [GroupContrastWithInteraction.from_string(contrast, subjects_df)]

    def _get_output_filename(self, contrast: GroupContrastWithInteraction) -> str:
        """Build the filename from class attributes and provided contrast.

        Parameters
        ----------
        contrast : GroupContrastWithInteraction
            The contrast to use for building the filename.
        """
        return (
            f"interaction-{contrast.name}_measure-{self.feature_label}_fwhm-{self.fwhm}"
        )
