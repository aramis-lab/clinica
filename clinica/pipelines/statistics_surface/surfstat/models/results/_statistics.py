from dataclasses import dataclass

import numpy as np
import pandas as pd
from brainstat.stats.SLM import SLM

from ._base import Results

__all__ = ["StatisticsResults"]


@dataclass
class PValueResults(Results):
    """This class implements a container for raw (uncorrected)
    P-value results obtained with a GLM model.

    Attributes
    ----------
    pvalues : np.ndarray
        Array of uncorrected P-values.

    mask : np.ndarray
        The binary mask.

    threshold : float
        The threshold used.
    """

    pvalues: np.ndarray
    mask: np.ndarray
    threshold: float

    @property
    def thresh(self):
        """For compatibility with previous Matlab implementation."""
        return self.threshold

    @property
    def P(self):
        """For compatibility with previous Matlab implementation."""
        return self.pvalues

    @classmethod
    def from_t_statistics(
        cls,
        tstats: np.ndarray,
        df: pd.DataFrame,
        mask: np.ndarray,
        threshold: float,
    ):
        """Instantiate the class from an array of T-statistics.

        Parameters
        ----------
        tstats : np.ndarray
            Array of T-statistics.

        df : pd.DataFrame
            The subjects DataFrame.

        mask : np.ndarray
            The binary mask.

        threshold : float
            The threshold to be used.
        """
        from scipy.stats import t

        return cls(1 - t.cdf(tstats, df), mask, threshold)


@dataclass
class CorrectedPValueResults(PValueResults):
    """This class implements a container for corrected P-value
    results obtained with a GLM model.

    Attributes
    ----------
    cluster_pvalues : np.ndarray
        The cluster P-values.
    """

    cluster_pvalues: np.ndarray

    @property
    def C(self):
        """For compatibility with previous Matlab implementation."""
        return self.cluster_pvalues


@dataclass
class StatisticsResults(Results):
    """This class implements a container for results obtained with
    the GLM model classes. It holds information relative to a GLM
    run with one specific contrast.

    Attributes
    ----------
    coefficients : np.ndarray
        The beta coefficients of the fitted GLM model.

    tstats : np.ndarray
        The corresponding T-statistics.

    uncorrected_p_value : PValueResults
        The corresponding uncorrected p values, stored in a `PValueResults` instance.

    fdr : np.ndarray
        The corresponding False Discovery Rate.

    corrected_p_value : CorrectedPValueResults
        The corresponding corrected p values, stored in a `CorrectedPValueResults` instance.
    """

    coefficients: np.ndarray
    tstats: np.ndarray
    uncorrected_p_values: PValueResults
    fdr: np.ndarray
    corrected_p_values: CorrectedPValueResults

    @property
    def TStatistics(self):
        """Needed for compatibility with previous implementation in Matlab."""
        return self.tstats

    @property
    def uncorrectedPValue(self):
        """Needed for compatibility with previous implementation in Matlab."""
        return self.uncorrected_p_values

    @property
    def correctedPValue(self):
        """Needed for compatibility with previous implementation in Matlab."""
        return self.corrected_p_values

    @property
    def FDR(self):
        """Needed for compatibility with previous implementation in Matlab."""
        return self.fdr

    @classmethod
    def from_slm_model(
        cls,
        model: SLM,
        mask: np.ndarray,
        threshold_uncorrected_p_value: float,
        threshold_corrected_p_value: float,
    ):
        """Instantiate from a SLM model.

        Parameters
        ----------
        model : brainstat.stats.SLM
            SLM model instance to use.

        mask : np.ndarray
            The binary mask to use.

        threshold_uncorrected_p_value : float
            The threshold to use with uncorrected P-values.

        threshold_corrected_p_value : float
            The threshold to use with corrected P-values.
        """
        idx = np.argwhere(np.isnan(model.t))
        corrected_pvals = model.P["pval"]["P"]
        corrected_pvals[idx] = 1.0
        tstats = np.nan_to_num(model.t)
        uncorrected_p_values = PValueResults.from_t_statistics(
            tstats,
            model.df,
            mask,
            threshold_uncorrected_p_value,
        )
        corrected_p_values = CorrectedPValueResults(
            corrected_pvals,
            mask,
            threshold_corrected_p_value,
            model.P["pval"]["C"],
        )
        return cls(
            np.nan_to_num(model.coef),
            tstats,
            uncorrected_p_values,
            model.Q,
            corrected_p_values,
        )
