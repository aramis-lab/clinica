from typing import List, Optional

import pandas as pd

from ._base import GLM
from ._contrast import CorrelationContrast


class CorrelationGLM(GLM):
    """Class implementing the correlation type GLM model.

    Attributes
    ----------
    See documentation for `GLM` class.

    group_label : str, optional
        The label to use for group GLM models. Default=None.
    """

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        group_label: Optional[str],
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
        self.with_interaction: bool = False
        self.group_label = group_label

    def _build_contrasts(
        self, contrast: str, subjects_df: pd.DataFrame
    ) -> List[CorrelationContrast]:
        return [CorrelationContrast.from_string(contrast, subjects_df)]

    def _get_output_filename(self, contrast: CorrelationContrast) -> str:
        """Build the filename from class attributes and provided contrast.

        Parameters
        ----------
        contrast : CorrelationContrast
            The contrast to use for building the filename.
        """
        return (
            f"group-{self.group_label}_correlation-{contrast.absolute_name}"
            f"-{contrast.sign}_measure-{self.feature_label}_fwhm-{self.fwhm}"
        )
