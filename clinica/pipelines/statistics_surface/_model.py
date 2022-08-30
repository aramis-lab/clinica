import abc
import warnings
from dataclasses import dataclass
from functools import reduce
from os import PathLike
from pathlib import Path
from string import Template
from typing import Callable, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect
from nilearn.surface import Mesh

from clinica.utils.stream import cprint

MISSING_TERM_ERROR_MSG = Template(
    "Term ${term} from the design matrix is not in the columns of the "
    "provided TSV file. Please make sure that there is no typo."
)


def _print_clusters(model: SLM, threshold: float) -> None:
    """This function prints the results related to total number
    of clusters, as well as the significative clusters.

    Parameters
    ----------
    model : brainstat.stats.SLM
        Fitted SLM model.

    threshold : float
        Cluster defining threshold.
    """
    cprint("#" * 40)
    cprint("After correction (Cluster-wise Correction for Multiple Comparisons): ")
    df = model.P["clus"][0]
    cprint(df)
    cprint(f"Clusters found: {len(df)}")
    cprint(
        f"Significative clusters (after correction): {len(df[df['P'] <= threshold])}"
    )


def _check_column_in_df(df: pd.DataFrame, column: str) -> None:
    """Checks if the provided column name is in the provided DataFrame.
    Raises a ValueError if not.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to analyze.

    column : str
        Name of the column to check.
    """
    if column not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=column))


def _categorical_column(df: pd.DataFrame, column: str) -> bool:
    """Returns `True` if the column is categorical and `False` otherwise.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to analyze.

    column : str
        The name of the column to check.

    Returns
    -------
    bool :
        `True` if the column contains categorical values, `False` otherwise.
    """
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
    design_matrix : str
        Design matrix specified as a string.

    df : pd.DataFrame
        Subjects DataFrame.

    Returns
    -------
    model : FixedEffect
        BrainStats model.
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
    return reduce(lambda x, y: x + y, model)


def _build_model_term(
    term: str,
    df: pd.DataFrame,
    add_intercept: Optional[bool] = True,
) -> FixedEffect:
    """Builds a BrainStats model term from the subjects
    DataFrame and a column name.

    Parameters
    ----------
    term : str
        The name of the column of the DataFrame to be used.

    df : pd.DataFrame
        The subjects DataFrame.

    add_intercept : bool
        If `True`, adds an intercept term.

    Returns
    -------
    FixedEffect :
        BrainStats model term.
    """
    return FixedEffect(df[term], add_intercept=add_intercept)


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

    def __init__(
        self,
        design: str,
        df: pd.DataFrame,
        feature_label: str,
        contrast: str,
        fwhm: Optional[int] = 20,
        threshold_uncorrected_pvalue: Optional[float] = 0.001,
        threshold_corrected_pvalue: Optional[float] = 0.05,
        cluster_threshold: Optional[float] = 0.001,
    ):
        self._two_tailed = False
        self._correction = ["fdr", "rft"]
        self.df = df
        self.feature_label = feature_label
        self.fwhm = fwhm
        self.threshold_uncorrected_pvalue = threshold_uncorrected_pvalue
        self.threshold_corrected_pvalue = threshold_corrected_pvalue
        self.cluster_threshold = cluster_threshold
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
        contrast : str
            Contrast in string format.
        """
        pass

    @abc.abstractmethod
    def filename_root(self, contrast: str):
        """Returns the output file name root for the provided contrast.

        .. note::
            This method needs to be implemented in subclasses.

        Parameters
        ----------
        contrast : str
            Contrast for which to get the output filename.
        """
        pass

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

    def save_results(self, output_dir: PathLike, method: Union[str, List[str]]) -> None:
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
        for contrast, result in self.results_.items():
            plotter = StatisticsResultsPlotter(
                Path(output_dir) / Path(self.filename_root(contrast)), mesh
            )
            for meth in method:
                plotter.plot(result, meth)


class CorrelationGLM(GLM):
    """Class implementing the correlation type GLM model.

    Attributes
    ----------
    See documentation for `GLM` class.

    group_label : str, optinal
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
        self.with_interaction = False
        self.absolute_contrast_name = None
        self.contrast_sign = None
        self.group_label = group_label
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

    def build_contrasts(self, contrast: str):
        """Build the contrast from the string specification.

        Parameters
        ----------
        contrast : str
            The contrast to build.
        """
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
        """Build the filename root part from class attributes and provided contrast.

        Parameters
        ----------
        contrast : str
            The contrast to use for building the filename.
        """
        if contrast not in self.contrasts:
            raise ValueError(f"Unknown contrast {contrast}.")
        return (
            f"group-{self.group_label}_correlation-{self.absolute_contrast_name}"
            f"-{self.contrast_sign}_measure-{self.feature_label}_fwhm-{self.fwhm}"
        )


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
        self.with_interaction = False
        self.group_label = group_label
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

    def build_contrasts(self, contrast: str):
        """Build the contrast from the string specification.

        Parameters
        ----------
        contrast : str
            The contrast to build.
        """
        _check_column_in_df(self.df, contrast)
        if not _categorical_column(self.df, contrast):
            raise ValueError(
                "Contrast should refer to a categorical variable for group comparison. "
                "Please select 'correlation' for 'glm_type' otherwise."
            )
        group_values = np.unique(self.df[contrast])
        for contrast_type, (i, j) in zip(["positive", "negative"], [(0, 1), (1, 0)]):
            contrast_name = f"{group_values[j]}-lt-{group_values[i]}"
            self.contrasts[contrast_name] = (
                self.df[contrast] == group_values[i]
            ).astype(int) - (self.df[contrast] == group_values[j]).astype(int)

    def filename_root(self, contrast: str):
        """Build the filename root part from class attributes and provided contrast.

        Parameters
        ----------
        contrast : str
            The contrast to use for building the filename.
        """
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
        self.with_interaction = True
        warnings.warn(
            "You included interaction as covariate in your model, "
            "please carefully check the format of your tsv files."
        )

    def build_contrasts(self, contrast: str):
        """Build the contrast from the string specification.

        Parameters
        ----------
        contrast : str
            The contrast to build.
        """
        contrast_elements = [_.strip() for _ in contrast.split("*")]
        for contrast_element in contrast_elements:
            _check_column_in_df(self.df, contrast_element)
        categorical = [_categorical_column(self.df, _) for _ in contrast_elements]
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
        built_contrast = self.df[continue_contrast].where(
            self.df[categorical_contrast] == group_values[0], 0
        ) - self.df[continue_contrast].where(
            self.df[categorical_contrast] == group_values[1], 0
        )
        self.contrasts[contrast] = built_contrast

    def filename_root(self, contrast: str):
        """Build the filename root part from class attributes and provided contrast.

        Parameters
        ----------
        contrast : str
            The contrast to use for building the filename.
        """
        if contrast not in self.contrasts:
            raise ValueError(f"Unknown contrast {contrast}.")
        return f"interaction-{contrast}_measure-{self.feature_label}_fwhm-{self.fwhm}"


def create_glm_model(
    glm_type: str,
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
    glm_type : str
        The type of GLM to be created. Either "correlation" or "group_comparison".

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
    if glm_type == "correlation":
        return CorrelationGLM(design, df, feature_label, contrast, **params)
    elif glm_type == "group_comparison":
        if "*" in design:
            return GroupGLMWithInteraction(
                design, df, feature_label, contrast, **params
            )
        return GroupGLM(design, df, feature_label, contrast, **params)
    raise ValueError(
        f"create_glm_model received an unknown GLM type: {glm_type}."
        f"Only 'correlation' and 'group_comparison' are supported."
    )


def _convert_arrays_to_lists(data: dict) -> dict:
    """If the input dictionary contains numpy arrays, this function will
    cast them to lists and return the same dictionary with the lists instead
    of the numpy arrays.

    Parameters
    ----------
    data : dict
        The dictionary to clean.

    Returns
    -------
    new_data : dict
        The dictionary with arrays casted to lists.
    """
    new_data = dict()
    for k, v in data.items():
        if isinstance(v, dict):
            new_data[k] = _convert_arrays_to_lists(v)
        elif isinstance(v, np.ndarray):
            new_data[k] = v.tolist()
        else:
            new_data[k] = v
    return new_data


class Results:
    """Common class for GLM results."""

    def to_dict(self) -> dict:
        """Returns the `Results` instance in dict format.

        Private attributes and all methods are not returned.

        This function does not perform any casting.

        Returns
        -------
        data : dict
            Resulting dictionary.
        """
        import inspect

        data = dict()
        for attribute in inspect.getmembers(self):
            name, value = attribute
            if not name.startswith("_"):
                if not inspect.ismethod(value):
                    if hasattr(value, "to_dict"):
                        data[name] = value.to_dict()
                    else:
                        data[name] = value
        return data

    def to_json(self, indent: Optional[int] = 4) -> str:
        """Returns the json of the `Results` instance.

        Parameters
        ----------
        indent : int, optional
            Indent to use. Default=4.

        Returns
        -------
        str :
            The JSON dumps of the results.
        """
        import json

        return json.dumps(_convert_arrays_to_lists(self.to_dict()), indent=indent)


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
        """Instanciate from a SLM model.

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
            model.P["pval"]["C"],
            mask,
            threshold_corrected_p_value,
        )
        return cls(
            np.nan_to_num(model.coef),
            tstats,
            uncorrected_p_values,
            model.Q,
            corrected_p_values,
        )


class StatisticsResultsPlotter:
    """Class responsible to plotting results of GLM fit.

    Attributes
    ----------
    output_file : PathLike
        Path to the output file.

    mesh : nilearn.surface.Mesh
        The mesh to be used for plotting results.
    """

    def __init__(self, output_file: PathLike, mesh: Mesh):
        self.output_file = output_file
        self.mesh = mesh
        self.plotting_extension = ".png"
        self.no_plot = {"coefficients"}  # Elements which should not be plotted

    def plot(self, result: StatisticsResults, method: str) -> None:
        """Plot the results.

        Parameters
        ----------
        result : StatisticsResults
            The results to be plotted.

        method : str
            The plotting method to use.
        """
        plotter = self._get_plotter(method)
        plotter(result)

    def _get_plotter(self, method: str) -> Callable[[StatisticsResults], None]:
        """Returns the plotting method from its name.

        Parameters
        ----------
        method : str
            Name of the plotting method to use.

        Returns
        -------
        Callable :
            Plotting method.
        """
        if method == "nilearn_plot_surf_stat_map":
            return self._plot_stat_maps
        else:
            raise NotImplementedError(f"Plotting method {method} is not implemented.")

    def _plot_stat_maps(self, result: StatisticsResults) -> None:
        """Wrapper around the `nilearn.plotting.plot_surf_stat_map` method.

        Parameters
        ----------
        result : StatisticsResults
            The results to plot.
        """
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
    output_file : PathLike
        Path and filename root to be used.
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
        result : StatisticsResults
            Results to be saved.

        method : str
            Name of the saving method to use.
        """
        writer = self._get_writer(method)
        writer(result)

    def _get_writer(self, method: str) -> Callable[[StatisticsResults], None]:
        """Returns a writter method from its name.

        Parameters
        ----------
        method : str
            The name of the writting method to use.

        Returns
        -------
        Callable :
            The writting method.
        """
        if method.lower() == "json":
            return self._write_to_json
        elif method.lower() == "mat":
            return self._write_to_mat
        else:
            raise NotImplementedError(
                f"Serializing method {method} is not implemented."
            )

    def _write_to_json(self, results: StatisticsResults) -> None:
        """Write the provided `StatisticsResults` to JSON format.

        Parameters
        ----------
        results : StatisticsResults
            The results to write to disk in JSON format.
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
            json.dump(results.to_json(indent=self.json_indent), fp)

    def _write_to_mat(self, results: StatisticsResults) -> None:
        """Write the provided `StatisticsResults` to MAT format.

        Parameters
        ----------
        results : StatisticsResults
            The results to write to disk in MAT format.
        """
        from scipy.io import savemat

        # These labels are used for compatibility with the previous
        # MATLAB implementation of the Statistics Surface Pipeline
        # of Clinica.
        struct_labels = {
            "coefficients": "coef",
            "TStatistics": "tvaluewithmask",
            "uncorrectedPValue": "uncorrectedpvaluesstruct",
            "correctedPValue": "correctedpvaluesstruct",
            "FDR": "FDR",
        }
        for name, res in results.to_dict().items():
            if name in struct_labels:
                mat_filename = str(self.output_file) + "_" + name + self.mat_extension
                cprint(
                    msg=f"Writing {name} results to MAT in  {mat_filename}",
                    lvl="info",
                )
                savemat(mat_filename, {struct_labels[name]: res})
