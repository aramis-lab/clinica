"""This file contains functions for loading data from disk
and performing some checks on them.
"""
from os import PathLike
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from nilearn.surface import Mesh

TSV_FIRST_COLUMN = "participant_id"
TSV_SECOND_COLUMN = "session_id"


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
        raise ValueError(f"The TSV data in {tsv_file} should have at least 2 columns.")
    if tsv_data.columns[0] != TSV_FIRST_COLUMN:
        raise ValueError(
            f"The first column in {tsv_file} should always be {TSV_FIRST_COLUMN}."
        )
    if tsv_data.columns[1] != TSV_SECOND_COLUMN:
        raise ValueError(
            f"The second column in {tsv_file} should always be {TSV_SECOND_COLUMN}."
        )
    return tsv_data


def _get_t1_freesurfer_custom_file_template(base_dir: PathLike) -> str:
    """Returns a Template for the path to the desired surface file."""
    return str(base_dir) + (
        "/%(subject)s/%(session)s/t1/freesurfer_cross_sectional/%(subject)s_%(session)s"
        "/surf/%(hemi)s.thickness.fwhm%(fwhm)s.fsaverage.mgh"
    )


def _build_thickness_array(
    input_dir: PathLike,
    surface_file: str,
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
        parts = []
        for hemi in ["lh", "rh"]:
            query = {"subject": subject, "session": session, "fwhm": fwhm, "hemi": hemi}
            parts.append(
                load(str(Path(input_dir) / Path(surface_file % query))).get_fdata()
            )
        combined = np.vstack(parts)
        thickness.append(combined.flatten())
    thickness = np.vstack(thickness)
    if thickness.shape[0] != len(df):
        raise ValueError(
            f"Unexpected shape for thickness array : {thickness.shape}. "
            f"Expected {len(df)} rows."
        )
    return thickness


def _get_average_surface(fsaverage_path: PathLike) -> Tuple[Dict, Mesh]:
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
    import copy

    from nilearn.surface import Mesh, load_surf_mesh

    meshes = [
        load_surf_mesh(str(fsaverage_path / Path(f"{hemi}.pial")))
        for hemi in ["lh", "rh"]
    ]
    coordinates = np.vstack([mesh.coordinates for mesh in meshes])
    faces = np.vstack(
        [meshes[0].faces, meshes[1].faces + meshes[0].coordinates.shape[0]]
    )
    average_mesh = Mesh(
        coordinates=coordinates,
        faces=copy.deepcopy(faces),
    )
    ##################
    # UGLY HACK !!! Need investigation
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
