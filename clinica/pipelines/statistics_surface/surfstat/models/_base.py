import abc
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect
from nilearn.surface import Mesh

from clinica.utils.stream import cprint

from ._contrast import Contrast
from ._utils import build_model, print_clusters
from .results import (
    StatisticsResults,
    StatisticsResultsPlotter,
    StatisticsResultsSerializer,
)

__all__ = ["GLM"]


class GLM:
    """This class implements the functionalities common to all GLM models
    used in the Clinica SurfaceStatistics pipeline.

    Attributes
    ----------
    design : str
        The design matrix specified in string format.
        If this contains a "*", it will be interpreted as an interaction effect.

    df : pd.DataFrame
        The subjects DataFrame.

    feature_label : str
        The label used for building output filenames.

    contrast : str
        The contrast specified in string format.

    fwhm : int, optional
        The smoothing FWHM. This is used in the output file names.
        Default=20.

    threshold_uncorrected_pvalue : float, optional
        The threshold to be used with uncorrected P-values. Default=0.001.

    threshold_corrected_pvalue : float, optional
        The threshold to be used with corrected P-values. Default=0.05.

    cluster_threshold : float, optional
        The threshold to be used to declare clusters as significant. Default=0.001.
    """

    default_threshold_uncorrected_pvalue = 0.001
    default_threshold_corrected_pvalue = 0.05
    default_cluster_threshold = 0.001

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        fwhm: Optional[int] = 20,
        threshold_uncorrected_pvalue: Optional[float] = None,
        threshold_corrected_pvalue: Optional[float] = None,
        cluster_threshold: Optional[float] = None,
    ):
        self._two_tailed: bool = False
        self._correction = ["fdr", "rft"]
        self.feature_label: str = feature_label
        self.fwhm: Optional[int] = fwhm
        self.threshold_uncorrected_pvalue: float = (
            threshold_uncorrected_pvalue or self.default_threshold_uncorrected_pvalue
        )
        self.threshold_corrected_pvalue: float = (
            threshold_corrected_pvalue or self.default_threshold_corrected_pvalue
        )
        self.cluster_threshold: float = (
            cluster_threshold or self.default_cluster_threshold
        )
        self.results_ = None
        self.slm_models_ = None
        self.filenames: dict = {}
        self.model: FixedEffect = build_model(design, df)
        self.contrasts: List[Contrast] = self._build_contrasts(contrast, df)

    @abc.abstractmethod
    def _build_contrasts(
        self, contrast: str, subjects_df: pd.DataFrame
    ) -> List[Contrast]:
        raise NotImplementedError

    @property
    def contrast_names(self) -> List[str]:
        if self.contrasts is not None:
            return [contrast.name for contrast in self.contrasts]
        return list()

    def get_contrast_by_name(self, contrast_name: str) -> Contrast:
        if self.contrasts:
            contrast = [c for c in self.contrasts if c.name == contrast_name]
            if len(contrast) == 1:
                return contrast[0]
        raise ValueError(f"Unknown contrast {contrast_name}")

    @property
    def results(self):
        if self._is_fitted():
            return self.results_

    def get_output_filename(self, contrast_name: str) -> str:
        """Returns the output file name root for the provided contrast.

        .. note::
            This method needs to be implemented in subclasses.

        Parameters
        ----------
        contrast_name : str
            Contrast for which to get the output filename.
        """
        contrast = self.get_contrast_by_name(contrast_name)
        return self._get_output_filename(contrast)

    @abc.abstractmethod
    def _get_output_filename(self, contrast: Contrast) -> str:
        raise NotImplementedError

    def _is_fitted(self) -> bool:
        return self.results_ is not None

    def fit(
        self, data: np.ndarray, surface: Dict, mask: Optional[np.ndarray] = None
    ) -> None:
        """Fit the GLM model instance.

        Parameters
        ----------
        data : np.ndarray
            The data on which to fit the GLM model.

        surface : dict
            The Brainstat surface on which to fit the GLM model.

        mask : np.ndarray, optional
            The mask to be used to mask the data. Default=None.
        """
        if mask is None:
            mask = data[0, :] > 0
        self.results_ = dict()
        self.slm_models_ = dict()
        for contrast in self.contrasts:
            slm_model = SLM(
                self.model,
                contrast=contrast.built_contrast,
                surf=surface,
                mask=mask,
                two_tailed=self._two_tailed,
                correction=self._correction,
                cluster_threshold=self.cluster_threshold,
            )
            cprint(
                msg=f"Fitting the GLM model with contrast {contrast.name}...",
                lvl="info",
            )
            slm_model.fit(data)
            print_clusters(slm_model, self.threshold_corrected_pvalue)
            self.results_[contrast.name] = StatisticsResults.from_slm_model(
                slm_model,
                mask,
                self.threshold_uncorrected_pvalue,
                self.threshold_corrected_pvalue,
            )
            self.slm_models_[contrast.name] = slm_model

    def save_results(self, output_dir: Path, method: Union[str, List[str]]) -> None:
        """Save results to the provided output directory.

        Parameters
        ----------
        output_dir : PathLike
            The output directory in which to write the results.

        method : str or List[str]
            The method(s) to write the results.
        """
        if not self._is_fitted():
            raise ValueError(
                "GLM model needs to be fitted before accessing the results."
            )
        if isinstance(method, str):
            method = [method]
        for contrast_name, result in self.results_.items():
            result_serializer = StatisticsResultsSerializer(
                output_dir / self.get_output_filename(contrast_name)
            )
            for meth in method:
                result_serializer.save(result, meth)

    def plot_results(
        self,
        output_dir: Path,
        method: Union[str, List[str]],
        mesh: Mesh,
    ) -> None:
        """Plot results to the provided directory.

        Parameters
        ----------
        output_dir : PathLike
            The output directory in which to write the plot files.

        method : str or List[str]
            The method(s) to make the plots.

        mesh : nilearn.surface.Mesh
            The mesh on which to plot the result data.
        """
        if not self._is_fitted():
            raise ValueError(
                "GLM model needs to be fitted before accessing the results."
            )
        if isinstance(method, str):
            method = [method]
        for contrast_name, result in self.results_.items():
            plotter = StatisticsResultsPlotter(
                output_dir / self.get_output_filename(contrast_name), mesh
            )
            for meth in method:
                plotter.plot(result, meth)
