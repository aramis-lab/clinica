"""This file contains functions for loading data from disk
and performing some checks on them.
"""

import warnings
import numpy as np
import pandas as pd
from os import PathLike
from pathlib import Path
from typing import Dict, Tuple
from string import Template


DEFAULT_FWHM = 20
DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE = 0.001
DEFAULT_THRESHOLD_CORRECTED_P_VALUE = 0.05
DEFAULT_CLUSTER_THRESHOLD = 0.001
TSV_FIRST_COLUMN = "participant_id"
TSV_SECOND_COLUMN = "session_id"


def _extract_parameters(parameters: Dict) -> Tuple[float, float, float, float]:
    """This function extracts the parameters from the dictionary passed
    to the main function. If required parameters are not provided, the
    default values will be used.

    Parameters
    ----------
    parameters : Dictionary of parameters passed to the main function.

    Returns
    -------
    fwhm : Smoothing
    threshold_uncorrected_pvalue : Threshold to be used with uncorrected P-values.
    threshold_corrected_pvalue : Threshold to be used with corrected P-values.
    cluster_threshold : Threshold to be used to declare clusters as significant.
    """
    fwhm = DEFAULT_FWHM
    if "sizeoffwhm" in parameters:
        fwhm = parameters["sizeoffwhm"]
    threshold_uncorrected_pvalue = DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE
    if "thresholduncorrectedpvalue" in parameters:
        threshold_uncorrected_pvalue = parameters["thresholduncorrectedpvalue"]
    threshold_corrected_pvalue = DEFAULT_THRESHOLD_CORRECTED_P_VALUE
    if "thresholdcorrectedpvalue" in parameters:
        threshold_corrected_pvalue = parameters["thresholdcorrectedpvalue"]
    cluster_threshold = DEFAULT_CLUSTER_THRESHOLD
    if "clusterthreshold" in parameters:
        cluster_threshold = parameters["clusterthreshold"]
    return fwhm, threshold_uncorrected_pvalue, threshold_corrected_pvalue, cluster_threshold


def _read_and_check_tsv_file(tsv_file: PathLike) -> pd.DataFrame:
    """This function reads the TSV file provided and performs some basic checks.

    Parameters
    ----------
    tsv_file : TSV file to open.

    Returns
    -------
    tsv_data : DataFrame obtained from the file.
    """
    if not Path(tsv_file).exists():
        raise FileNotFoundError(f"File {tsv_file} does not exist.")
    tsv_data = pd.read_csv(tsv_file, sep="\t")
    if len(tsv_data.columns) < 2:
        raise ValueError(
            f"The TSV data in {tsv_file} should have at least 2 columns."
        )
    if tsv_data.columns[0] != TSV_FIRST_COLUMN:
        raise ValueError(
            f"The first column in {tsv_file} should always be {TSV_FIRST_COLUMN}."
        )
    if tsv_data.columns[1] != TSV_SECOND_COLUMN:
        raise ValueError(
            f"The second column in {tsv_file} should always be {TSV_SECOND_COLUMN}."
        )
    return tsv_data


def _get_t1_freesurfer_custom_file_template(base_dir: PathLike) -> Template:
    """Returns a Template for the path to the desired surface file."""
    return Template(
        str(base_dir) + (
        "/${subject}/${session}/t1/freesurfer_cross_sectional/${subject}_${session}"
        "/surf/${hemi}.thickness.fwhm${fwhm}.fsaverage.mgh"
        )
    )


def _build_thickness_array(
        input_dir: PathLike,
        surface_file: Template,
        df: pd.DataFrame,
        fwhm: float,
) -> np.ndarray:
    """This function builds the cortical thickness array.

    Parameters
    ----------
    input_dir : Input directory.
    surface_file : Template for the path to the surface file of interest.
    df : Subjects DataFrame
    fwhm : Smoothing parameter only used to retrieve the right surface file.

    Returns
    -------
    thickness : Cortical thickness. Hemispheres and subjects are stacked.
    """
    from nibabel.freesurfer.mghformat import load
    thickness = []
    for idx, row in df.iterrows():
        subject = row[TSV_FIRST_COLUMN]
        session = row[TSV_SECOND_COLUMN]
        parts = (
            load(
                Path(input_dir) / Path(surface_file.safe_substitute(
                    subject=subject, session=session, fwhm=fwhm, hemi=hemi
                ))
            ).get_fdata() for hemi in ['lh', 'rh']
        )
        combined = np.vstack(parts)
        thickness.append(combined.flatten())
    thickness = np.vstack(thickness)
    if thickness.shape[0] != len(df):
        raise ValueError(
            f"Unexpected shape for thickness array : {thickness.shape}. "
            f"Expected {len(df_subjects)} rows."
        )
    return thickness


def _get_average_surface(fsaverage_path: PathLike) -> dict:
    """This function extracts the average surface and the average mesh
    from the path to the fsaverage templates.

    .. note::

        Note that the average surface is returned as a dictionary
        with 'coord' and 'tri' as keys, while the average mesh is
        returned as a Nilearn Mesh object (basically a NamedTuple
        with 'coordinates' and 'faces' attributes). The surface
        isn't returned as a Nilearn Surface object for compatibility
        with BrainStats.

    .. warning::

        There is an issue with faces having  a value of 0 as index.
        This is most likely a bug in BrainStat as MATLAB indexing
        starts at 1 while Python starts at zero.

    Parameters
    ----------
    fsaverage_path : Path to the fsaverage templates.

    Returns
    -------
    average_surface : Average surface as a dictionary for BrainStat compatibility.
    average_mesh : Average mesh as a Nilearn Mesh object.
    """
    from nilearn.surface import Mesh, load_surf_mesh
    meshes = [
        load_surf_mesh(str(fsaverage_path / Path(f"{hemi}.pial")))
        for hemi in ['lh', 'rh']
    ]
    coordinates = np.vstack([mesh.coordinates for mesh in meshes])
    faces = np.vstack([
        meshes[0].faces,
        meshes[1].faces + meshes[0].coordinates.shape[0]
    ])
    average_mesh = Mesh(
        coordinates=coordinates,
        faces=faces,
    )
    ##################
    ## UGLY HACK !!! Need investigation
    ##################
    # Uncomment the following line if getting an error
    # with negative values in bincount in Brainstat.
    # Not sure, but might be a bug in BrainStat...
    #
    faces += 1
    #################
    average_surface = {
        "coord": coordinates,
        "tri": faces,
    }
    return average_surface, average_mesh


def _check_contrast(
        contrast: str,
        df: pd.DataFrame,
        glm_type: str,
) -> Tuple[str, str, bool]:
    """This function performs some basic checks on the provided contrast.

    Parameters
    ----------
    contrast : Contrast in string format.
    df : Subject DataFrame. It must contain the contrast elements as columns.
    glm_type : The type of GLM run. Can be 'group_comparison' or 'correlation'.

    Returns
    -------
    absolute_contrast : Contrast without the negative sign if any.
    contrat_sign : Either 'positive' or 'negative' depending on the provided contrast.
    with_interaction : Boolean indicating whether the contrast has interaction terms or not.
    """
    absolute_contrast = contrast
    with_interaction = False
    contrast_sign = "positive"
    if contrast.startswith("-"):
        absolute_contrast = contrast[1:].lstrip()
        contrast_sign = "negative"
    if "*" in contrast:
        with_interaction = True
        warnings.warn(
            "You included interaction as covariate in your model, "
            "please carefully check the format of your tsv files."
        )
    else:
        if absolute_contrast not in df.columns:
            raise ValueError(
                f"Column {absolute_contrast} does not exist in provided TSV file."
            )
        if glm_type == "group_comparison":
            unique_labels = np.unique(df[absolute_contrast])
            if len(unique_labels) != 2:
                raise ValueError(
                    "For group comparison, there should be just 2 different groups!"
                )
    return absolute_contrast, contrast_sign, with_interaction

