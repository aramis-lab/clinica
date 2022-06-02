"""This file contains functions to write results of the
StatisticsSurface Pipeline to disk.
"""

import warnings
import numpy as np
from os import PathLike
from typing import Dict, List, Union


def _save_results_to_json(
        results: Dict,
        filename_root: PathLike,
        verbose: bool = True
):
    """Write the provided results to JSON format.

    Parameters
    ----------
    results : Results to write to disk.
    filename_root : Filename root common to all output files.
    verbose : Verbose mode.
    """
    import json
    out_json_file = str(filename_root) + "_results.json"
    if verbose:
        print(f"--> Writing results to JSON in {out_json_file}...")
    jsonable = dict()
    for name, struct in results.items():
        if isinstance(struct, np.ndarray):
            jsonable[name] = struct.tolist()
        elif isinstance(struct, dict):
            jsonable[name] = dict()
            for k, v in struct.items():
                if isinstance(v, np.ndarray):
                    jsonable[name][k] = v.tolist()
                else:
                    jsonable[name][k] = v
        else:
            jsonable[name] = struct
    with open(out_json_file, "w") as fp:
        json.dump(jsonable, fp, indent=4)


def _save_results_to_mat(
        results: Dict,
        filename_root: PathLike,
        verbose: bool = True,
):
    """Write the provided results to mat format.

    .. note::

        There are much more convenient ways to write
        arrays to disk in Python. The main reason this
        function exists is for providing backward
        compatibility with previous versions of the
        Clinica's StatisticsSurface pipeline.

    Parameters
    ----------
    results : Results to write to disk.
    filename_root : Filename root common to all output files.
    verbose : Verbose mode.
    """
    if verbose:
        print("--> Writing results to mat files...")
    # These labels are used for compatibility with the previous
    # MATLAB implementation of the Statistics Surface Pipeline
    # of Clinica.
    STRUCT_LABELS = {
        "coefficients": "coef",
        "TStatistics": "tvaluewithmask",
        "uncorrectedPValue": "uncorrectedpvaluesstruct",
        "correctedPValue": "correctedpvaluesstruct",
        "FDR": "FDR",
    }
    for name, result in results.items():
        _save_to_mat(
            result,
            str(filename_root) + "_" + name,
            STRUCT_LABELS[name],
            verbose=verbose,
        )


def _save_results_to_bids(
        results: Dict,
        filename_root: PathLike,
        verbose: bool = True
):
    """Write provided results to BIDS format.

    .. warning::
        This is not implemented yet.

    """
    warnings.warn("Writing results to BIDS is not implemented yet.")


WRITERS = {
    "json": _save_results_to_json,
    "mat": _save_results_to_mat,
    "bids": _save_results_to_bids,
}


def _save_results(
        results: Dict,
        filename_root: PathLike,
        out_formats: Union[str, List] = "all",
        verbose: bool = True
):
    """Write the provided results to all requested output formats.

    Parameters
    ----------
    results : Results to write to disk.
    filename_root : Filename root common to all output files.
    out_formats : Either a list of output formats (among "mat", "json", and "bids"), or "all".
    verbose : Verbose mode.
    """
    if out_formats == "all":
        out_formats = list(WRITERS.keys())
    for output_format in out_formats:
        if output_format not in WRITERS:
            warnings.warn(
                f"Could not write to {output_format} because writer doesn't exist."
            )
        WRITERS[output_format](results, filename_root, verbose=verbose)


def _print_clusters(model, threshold: float):
    """This function prints the results related to total number
    of clusters, as well as the significative clusters.

    Parameters
    ----------
    model : Fitted SLM model.
    threshold : Cluster defining threshold.
    """
    print("#" * 40)
    print("After correction (Clusterwise Correction for Multiple Comparisons): ")
    df = model.P['clus'][1]
    print(df)
    print(f"Clusters found: {len(df)}")
    print(f"Significative clusters (after correction): {len(df[df['P'] <= threshold])}")


def _plot_stat_map(
        mesh: np.ndarray,
        texture: np.ndarray,
        filename: str,
        threshold: float = None,
        title: str = None,
        verbose: bool = True
):
    """Plot a given texture over the provided mesh using Nilearn's
    plot_surf_stat_map function.

    Parameters
    ----------
    mesh : Mesh of the surface to plot.
    texture : Texture of the surface to plot.
    threshold : Threshold to be used for plotting.
    title : Title to display on the plot.
    verbose : Verbose mode.
    """
    try:
        from nilearn.plotting import plot_surf_stat_map
    except ImportError:
        raise ImportError(
            "Nilearn is required to plot the surfaces obtained with this pipeline."
        )
    plot_filename = filename + ".png"
    if verbose:
        print(f"--> Saving plot to {plot_filename}")
    plot_surf_stat_map(
        mesh, texture, threshold=threshold, output_file=plot_filename, title=title,
    )


def _plot_results(
        results: Dict,
        filename_root: PathLike,
        mesh: np.ndarray,
        verbose: bool = True
):
    """This function will plot all possible surfaces in the
    provided results dictionary.

    Parameters
    ----------
    results : Dictionary of results.
    filename_root : The common root for the output filenames.
    mesh : The mesh to be used for plotting.
    verbose : Verbose mode.
    """
    RESULTS_NO_PLOT = set(["coefficients"])
    for name, result in results.items():
        if name not in RESULTS_NO_PLOT:
            if isinstance(result, dict):
                texture = result['P']
            else:
                texture = result
            _plot_stat_map(
                mesh, texture, str(filename_root) + name,
                threshold=None, title=name, verbose=verbose,
            )


def _save_to_mat(struct: Dict, filename: str, key: str, verbose: bool = True):
    """Write a given struct to a mat file.

    Parameters
    ----------
    struct : Dictinary to write.
    filename : File name for writing (without the '.mat' extension).
    key : Key to be used to refer to the provided struct.
    verbose : Verbose mode.
    """
    from scipy.io import savemat
    mat_filename = filename + ".mat"
    if verbose:
        print(f"--> Saving matrix to {mat_filename}")
    savemat(mat_filename, {key: struct})

