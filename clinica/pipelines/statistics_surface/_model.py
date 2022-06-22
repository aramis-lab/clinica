import abc
import warnings
from functools import reduce
from os import PathLike
from pathlib import Path
from string import Template
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect
from nilearn.surface import Mesh

from clinica.utils.stream import cprint

DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE = 0.001
DEFAULT_THRESHOLD_CORRECTED_P_VALUE = 0.05
DEFAULT_CLUSTER_THRESHOLD = 0.001
MISSING_TERM_ERROR_MSG = Template(
    "Term ${term} from the design matrix is not in the columns of the "
    "provided TSV file. Please make sure that there is no typo."
)


def _print_clusters(model: SLM, threshold: float) -> None:
    """This function prints the results related to total number
    of clusters, as well as the significative clusters.

    Parameters
    ----------
    model : Fitted SLM model.
    threshold : Cluster defining threshold.
    """
    cprint("#" * 40)
    cprint("After correction (Cluster-wise Correction for Multiple Comparisons): ")
    df = model.P["clus"][1]
    cprint(df)
    cprint(f"Clusters found: {len(df)}")
    cprint(
        f"Significative clusters (after correction): {len(df[df['P'] <= threshold])}"
    )


def _is_categorical(df: pd.DataFrame, column: str) -> bool:
    """This function returns whether the provided DataFrame's column
    is categorical or not.

    .. note::
        There might be more clever ways to do that.

    Parameters
    ----------
    df : DataFrame to analyze.
    column : Name of the column to check.

    Returns
    -------
    True if the column contains categorical values, False otherwise.
    """
    if column not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=column))
    return not df[column].dtype.name.startswith("float")


def _build_model(design_matrix: str, df: pd.DataFrame) -> FixedEffect:
    """Build a brainstat model from the design matrix in
    string format.

    This function assumes that the design matrix is formatted
    in the following way:

        1 + factor_1 + factor_2 + ...

    Or:

        factor_1 + factor_2 + ...

    in the latter case the intercept will be added automatically.

    Parameters
    ----------
    design_matrix : Design matrix specified as a string.
    df : Subjects DataFrame.

    Returns
    -------
    model : BrainStats model.
    """
    if len(design_matrix) == 0:
        raise ValueError("Design matrix cannot be empty.")
    if "+" in design_matrix:
        terms = [_.strip() for _ in design_matrix.split("+")]
    else:
        terms = [design_matrix.strip()]
    model = []
    for term in terms:
        # Intercept is automatically included in brainstat
        if term == "1":
            continue
        # Handles the interaction effects
        if "*" in term:
            sub_terms = [_.strip() for _ in term.split("*")]
            model_term = reduce(
                lambda x, y: x * y, [_build_model_term(_, df) for _ in sub_terms]
            )
        else:
            model_term = _build_model_term(term, df)
        model.append(model_term)
    if len(model) == 1:
        return model[0]
    return reduce(lambda x, y: x + y, model)


def _build_model_term(term: str, df: pd.DataFrame) -> FixedEffect:
    """Builds a BrainStats model term from the subjects
    DataFrame and a column name.

    Parameters
    ----------
    term : Name of the column of the DataFrame to be used.
    df : Subjects DataFrame.

    Returns
    -------
    BrainStats FixedEffect.
    """
    if term not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=term))
    return FixedEffect(df[term])


class GLM:
    """This class implements the functionalities common to all GLM models
    used in the Clinica SurfaceStatistics pipeline.
    """

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        **kwargs,
    ):
        self._two_tailed = True  # Could be exposed to users?
        self._correction = ["fdr", "rft"]  # Could be exposed to users?
        self.df = df
        self.feature_label = feature_label
        self.fwhm = kwargs.pop("sizeoffwhm")
        self.threshold_uncorrected_pvalue = kwargs.pop(
            "thresholduncorrectedpvalue",
            DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE,
        )
        self.threshold_corrected_pvalue = kwargs.pop(
            "thresholdcorrectedpvalue",
            DEFAULT_THRESHOLD_CORRECTED_P_VALUE,
        )
        self.cluster_threshold = kwargs.pop(
            "clusterthreshold",
            DEFAULT_CLUSTER_THRESHOLD,
        )
        self.results_ = None
        self.slm_models_ = None
        self.contrasts = dict()
        self.filenames = dict()
        self.model = _build_model(design, df)
        self.build_contrasts(contrast)

    @property
    def contrast_names(self) -> List[str]:
        if self.contrasts is not None:
            return list(self.contrasts.keys())
        return list()

    @property
    def results(self):
        if self._is_fitted():
            return self.results_

    @abc.abstractmethod
    def build_contrasts(self, contrast: str):
        """Build the contrasts from the provided contrast in string format.

        .. note::
            This method needs to be implemented in subclasses.

        Parameters
        ----------
        contrast: Contrast in string format.
        """
        pass

    @abc.abstractmethod
    def filename_root(self, contrast: str):
        """Returns the output file name root for the provided contrast.

        .. note::
            This method needs to be implemented in subclasses.

        Parameters
        ----------
        contrast: Contrast for which to get the output filename.
        """
        pass

    def _is_fitted(self):
        return self.results_ is not None

    def fit(self, data: np.ndarray, surface: Dict, mask: Optional[np.ndarray] = None):
        if mask is None:
            mask = data[0, :] > 0
        self.results_ = dict()
        self.slm_models_ = dict()
        for contrast_name, contrast in self.contrasts.items():
            slm_model = SLM(
                self.model,
                contrast=contrast,
                surf=surface,
                mask=mask,
                two_tailed=self._two_tailed,
                correction=self._correction,
                cluster_threshold=self.cluster_threshold,
            )
            cprint(
                msg=f"Fitting the GLM model with contrast {contrast_name}...",
                lvl="info",
            )
            slm_model.fit(data)
            _print_clusters(slm_model, self.threshold_corrected_pvalue)
            self.results_[contrast_name] = StatisticsResults.from_slm_model(
                slm_model,
                mask,
                self.threshold_uncorrected_pvalue,
                self.threshold_corrected_pvalue,
            )
            self.slm_models_[contrast_name] = slm_model

    def save_results(self, output_dir: PathLike, method: Union[str, List[str]]):
        """Save results to the provided output directory."""
        if not self._is_fitted():
            raise ValueError(
                "GLM model needs to be fitted before accessing the results."
            )
        if isinstance(method, str):
            method = [method]
        for contrast, result in self.results_.items():
            result_serializer = StatisticsResultsSerializer(
                Path(output_dir) / Path(self.filename_root(contrast))
            )
            for meth in method:
                result_serializer.save(result, meth)

    def plot_results(
        self,
        output_dir: PathLike,
        method: Union[str, List[str]],
        mesh: Mesh,
    ):
        """Plot results to the provided directory."""
        if not self._is_fitted():
            raise ValueError(
                "GLM model needs to be fitted before accessing the results."
            )
        if isinstance(method, str):
            method = [method]
        for contrast, result in self.results_.items():
            plotter = StatisticsResultsPlotter(
                Path(output_dir) / Path(self.filename_root(contrast)), mesh
            )
            for meth in method:
                plotter.plot(result, meth)


class CorrelationGLM(GLM):
    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        **kwargs,
    ):
        self.with_interaction = False
        self.absolute_contrast_name = None
        self.contrast_sign = None
        self.group_label = kwargs.pop("group_label")
        super().__init__(design, df, feature_label, contrast, **kwargs)

    def build_contrasts(self, contrast: str):
        absolute_contrast_name = contrast
        contrast_sign = "positive"
        if contrast.startswith("-"):
            absolute_contrast_name = contrast[1:].lstrip()
            contrast_sign = "negative"
        built_contrast = self.df[absolute_contrast_name]
        if contrast_sign == "negative":
            built_contrast *= -1
        self.contrasts[contrast] = built_contrast
        self.absolute_contrast_name = absolute_contrast_name
        self.contrast_sign = contrast_sign

    def filename_root(self, contrast: str):
        if contrast not in self.contrasts:
            raise ValueError(f"Unknown contrast {contrast}.")
        return (
            f"group-{self.group_label}_correlation-{self.absolute_contrast_name}"
            f"-{self.contrast_sign}_measure-{self.feature_label}_fwhm-{self.fwhm}"
        )


class GroupGLM(GLM):
    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        **kwargs,
    ):
        self.with_interaction = False
        self.group_label = kwargs.pop("group_label", "group")
        super().__init__(design, df, feature_label, contrast, **kwargs)

    def build_contrasts(self, contrast: str):
        if not _is_categorical(self.df, contrast):
            raise ValueError(
                "Contrast should refer to a categorical variable for group comparison. "
                "Please select 'correlation' for 'glm_type' otherwise."
            )
        group_values = np.unique(self.df[contrast])
        for contrast_type, (i, j) in zip(["positive", "negative"], [(0, 1), (1, 0)]):
            contrast_name = f"{group_values[i]}-lt-{group_values[j]}"
            self.contrasts[contrast_name] = (
                self.df[contrast] == group_values[i]
            ).astype(int) - (self.df[contrast] == group_values[j]).astype(int)

    def filename_root(self, contrast: str):
        if contrast not in self.contrasts:
            raise ValueError(f"Unknown contrast {contrast}.")
        return f"group-{self.group_label}_{contrast}_measure-{self.feature_label}_fwhm-{self.fwhm}"


class GroupGLMWithInteraction(GroupGLM):
    """This class implements a GLM model for group comparison with
    interaction effects.

    Attributes
    ----------
    See attributes of parent class `GroupGLM`.
    """

    def __init__(
        self, design: str, df: pd.DataFrame, feature_label: str, contrast: str, **kwargs
    ):
        super().__init__(design, df, feature_label, contrast, **kwargs)
        self.with_interaction = True

    def build_contrasts(self, contrast: str):
        contrast_elements = [_.strip() for _ in contrast.split("*")]
        categorical = [_is_categorical(self.df, _) for _ in contrast_elements]
        if len(contrast_elements) != 2 or sum(categorical) != 1:
            raise ValueError(
                "The contrast must be an interaction between one continuous "
                "variable and one categorical variable. Your contrast contains "
                f"the following variables : {contrast_elements}"
            )
        idx = 0 if categorical[0] else 1
        categorical_contrast = contrast_elements[idx]
        continue_contrast = contrast_elements[(idx + 1) % 2]
        group_values = np.unique(self.df[categorical_contrast])
        built_contrast = self.df[continue_contrast] * (
            (self.df[categorical_contrast] == group_values[0]).astype(int)
        ) - self.df[continue_contrast] * (
            (self.df[categorical_contrast] == group_values[1]).astype(int)
        )
        self.contrasts[contrast] = built_contrast

    def filename_root(self, contrast: str):
        if contrast not in self.contrasts:
            raise ValueError(f"Unknown contrast {contrast}.")
        return f"interaction-{contrast}_measure-{self.feature_label}_fwhm-{self.fwhm}"


class GLMFactory:
    """Factory class for building GLM models.

    Attributes
    ----------
    feature_label: Label used for building output filenames.
    """

    def __init__(self, feature_label: str) -> None:
        self.feature_label = feature_label

    def create_model(
        self,
        glm_type: str,
        design: str,
        df: pd.DataFrame,
        contrast: str,
        **kwargs,
    ) -> GLM:
        """Factory method for building a GLM model instance corresponding to the
        provided type and design matrix.

        Parameters
        ----------
        glm_type: Type of GLM to be created. Either "correlation" or "group_comparison".
        design: Design matrix in string format. If this contains a "*", it will be
            interpreted as an interaction effect.
        df: Subjects DataFrame.
        contrast: Contrast in string format.
        """
        cprint(
            msg=f"The GLM model is: {design} and the GLM type is: {glm_type}",
            lvl="info",
        )
        if glm_type == "correlation":
            return CorrelationGLM(design, df, self.feature_label, contrast, **kwargs)
        elif glm_type == "group_comparison":
            if "*" in design:
                warnings.warn(
                    "You included interaction as covariate in your model, "
                    "please carefully check the format of your tsv files."
                )
                return GroupGLMWithInteraction(
                    design, df, self.feature_label, contrast, **kwargs
                )
            return GroupGLM(design, df, self.feature_label, contrast, **kwargs)
        raise ValueError(
            f"GLM factory received an unknown GLM type: {glm_type}."
            f"Only 'correlation' and 'group_comparison' are supported."
        )


class PValueResults:
    """This class implements a container for raw (uncorrected)
    P-value results obtained with a GLM model.

    Attributes
    ----------
    pvalues: Array of uncorrected P-values.
    mask: Binary mask.
    thresh: Threshold.
    """

    def __init__(
        self,
        pvalues: np.ndarray,
        mask: np.ndarray,
        threshold: float,
    ) -> None:
        self.P = pvalues
        self.mask = mask
        self.thresh = threshold

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
        tstats: Array of T-statistics.
        df: Subjects DataFrame.
        mask: Binary mask.
        threshold: Threshold.
        """
        from scipy.stats import t

        return cls(1 - t.cdf(tstats, df), mask, threshold)

    def to_dict(self, jsonable: bool = True):
        """Returns the `PValueResults` instance in dict format.

        Parameters
        ----------
        jsonable: If `True`, the output can be directly written in JSON.
            Otherwise, it might contain non-serializable objects.
        """
        import inspect

        json_dict = dict()
        for attribute in inspect.getmembers(self):
            name, value = attribute
            if not name.startswith("_"):
                if not inspect.ismethod(value):
                    if isinstance(value, np.ndarray) and jsonable:
                        json_dict[name] = value.tolist()
                    else:
                        json_dict[name] = value
        return json_dict


class CorrectedPValueResults(PValueResults):
    """This class implements a container for corrected P-value
    results obtained with a GLM model.

    Attributes
    ----------

    """

    def __init__(
        self,
        pvalues: np.ndarray,
        cluster_pvalues: np.ndarray,
        mask: np.ndarray,
        threshold: float,
    ) -> None:
        super().__init__(pvalues, mask, threshold)
        self.C = cluster_pvalues


class StatisticsResults:
    """This class implements a container for results obtained with
    the GLM model classes. It holds information relative to a GLM
    run with one specific contrast.

    Attributes
    ----------
    coefficients: The beta coefficients of the fitted GLM model.
    tstats: The corresponding T-statistics.
    uncorrected_p_value: The corresponding uncorrected p values,
        stored in a `PValueResults` instance.
    FDR: The corresponding False Discovery Rate.
    corrected_p_value: The corresponding corrected p values,
        stored in a `CorrectedPValueResults` instance.
    """

    def __init__(
        self,
        coefficients: np.ndarray,
        tstats: np.ndarray,
        uncorrected_p_values: PValueResults,
        fdr: np.ndarray,
        corrected_p_values: CorrectedPValueResults,
    ) -> None:
        self.coefficients = coefficients
        self.TStatistics = tstats
        self.uncorrectedPValues = uncorrected_p_values
        self.FDR = fdr
        self.correctedPValues = corrected_p_values

    @classmethod
    def from_slm_model(
        cls,
        model: SLM,
        mask: np.ndarray,
        threshold_uncorrected_p_value: float,
        threshold_corrected_p_value: float,
    ):
        tstats = np.nan_to_num(model.t)
        uncorrected_p_values = PValueResults.from_t_statistics(
            tstats,
            model.df,
            mask,
            threshold_uncorrected_p_value,
        )
        corrected_p_values = CorrectedPValueResults(
            model.P["pval"]["P"],
            model.P["pval"]["C"],
            mask,
            threshold_corrected_p_value,
        )
        return cls(
            np.nan_to_num(model.coef),
            tstats,
            uncorrected_p_values,
            model._fdr(),
            corrected_p_values,
        )

    def to_dict(self, jsonable: bool = True):
        """Returns the results in dict format.

        Parameters
        ----------
        jsonable: If `True`, the output can be directly written in JSON.
            Otherwise, it might contain non-serializable objects.
        """
        import inspect

        json_dict = dict()
        for attribute in inspect.getmembers(self):
            name, value = attribute
            if not name.startswith("_"):
                if not inspect.ismethod(value):
                    if hasattr(value, "to_dict"):
                        json_dict[name] = value.to_dict(jsonable=jsonable)
                    elif isinstance(value, np.ndarray) and jsonable:
                        json_dict[name] = value.tolist()
                    else:
                        json_dict[name] = value
        return json_dict


class StatisticsResultsPlotter:
    def __init__(self, output_file: PathLike, mesh: Mesh):
        self.output_file = output_file
        self.mesh = mesh
        self.plotting_extension = ".png"
        self.no_plot = {"coefficients"}  # Elements which should not be plotted

    def plot(self, result: StatisticsResults, method: str):
        plotter = self._get_plotter(method)
        return plotter(result)

    def _get_plotter(self, method: str):
        if method == "nilearn_plot_surf_stat_map":
            return self._plot_stat_maps
        else:
            raise NotImplementedError(f"Plotting method {method} is not implemented.")

    def _plot_stat_maps(self, result: StatisticsResults):
        from nilearn.plotting import plot_surf_stat_map

        for name, res in result.to_dict(jsonable=False).items():
            if name not in self.no_plot:
                texture = res
                threshold = None
                plot_filename = (
                    str(self.output_file) + "_" + name + self.plotting_extension
                )
                if isinstance(res, dict):
                    texture = res["P"]
                    threshold = res["thresh"]
                cprint(msg=f"Saving plot to {plot_filename}", lvl="info")
                plot_surf_stat_map(
                    self.mesh,
                    texture,
                    threshold=threshold,
                    output_file=plot_filename,
                    title=name,
                )


class StatisticsResultsSerializer:
    """This class is responsible for writing instances of `StatisticsResults`
    to disk through different methods.

    Attributes
    ----------
    output_file: Path and filename root to be used.
    """

    def __init__(self, output_file: PathLike):
        self.output_file = output_file
        self.json_extension = "_results.json"
        self.json_indent = 4
        self.mat_extension = ".mat"

    def save(self, result: StatisticsResults, method: str) -> None:
        """Save provided `StatisticsResults` to disk with provided method.

        Parameters
        ----------
        result: StatisticsResults to be saved.
        method: Saving method to use.
        """
        writer = self._get_writer(method)
        return writer(result)

    def _get_writer(self, method: str):
        if method.lower() == "json":
            return self._write_to_json
        elif method.lower() == "mat":
            return self._write_to_mat
        else:
            raise NotImplementedError(
                f"Serializing method {method} is not implemented."
            )

    def _write_to_json(self, result: StatisticsResults):
        """Write the provided `StatisticsResults` to JSON format.

        Parameters
        ----------
        results : Results to write to disk in JSON format.
        """
        import json
        import os

        out_json_file = Path(str(self.output_file) + self.json_extension)
        if not os.path.exists(out_json_file.parents[0]):
            os.makedirs(out_json_file.parents[0])
        cprint(
            msg=f"Writing results to JSON in {out_json_file}...",
            lvl="info",
        )
        with open(out_json_file, "w") as fp:
            json.dump(result.to_dict(jsonable=True), fp, indent=self.json_indent)

    def _write_to_mat(self, result: StatisticsResults):
        """Write the provided `StatisticsResults` to MAT format.

        Parameters
        ----------
        results: Results to write to disk in MAT format.
        """
        from scipy.io import savemat

        # These labels are used for compatibility with the previous
        # MATLAB implementation of the Statistics Surface Pipeline
        # of Clinica.
        struct_labels = {
            "coefficients": "coef",
            "TStatistics": "tvaluewithmask",
            "uncorrectedPValues": "uncorrectedpvaluesstruct",
            "correctedPValues": "correctedpvaluesstruct",
            "FDR": "FDR",
        }
        for name, res in result.to_dict(jsonable=False).items():
            mat_filename = str(self.output_file) + "_" + name + self.mat_extension
            cprint(
                msg=f"Writing {name} results to MAT in  {mat_filename}",
                lvl="info",
            )
            savemat(mat_filename, {struct_labels[name]: res})
