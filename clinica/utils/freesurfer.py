"""This module contains FreeSurfer utilities."""
from enum import Enum, auto
from pathlib import PurePath
from typing import List, Optional, Tuple

import pandas as pd


class ColumnType(Enum):
    PARCELLATION = auto()  # The actual value is meaningless.
    SEGMENTATION = auto()


def _get_prefix(subject_id: str) -> str:
    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)

    return "_".join([i for i in image_id if i])


def _read_stats_file(stats_filename: PurePath, column_type: ColumnType) -> pd.DataFrame:
    """Read the provided statistics file and extract the relevant data.

    Raises
    ------
    ValueError
        If `column_type` is not supported.

    FileNotFoundError
        If `stats_filename` is not a valid file name.
    """
    if column_type == ColumnType.PARCELLATION:
        # Columns in ?h.BA.stats, ?h.aparc.stats or ?h.aparc.a2009s.stats file
        columns = [
            "StructName",
            "NumVert",
            "SurfArea",
            "GrayVol",
            "ThickAvg",
            "ThickStd",
            "MeanCurv",
            "GausCurv",
            "FoldInd",
            "CurvInd",
        ]
    elif column_type == ColumnType.SEGMENTATION:
        # Columns in aseg.stats or wmparc.stats file
        columns = [
            "Index",
            "SegId",
            "NVoxels",
            "Volume_mm3",
            "StructName",
            "normMean",
            "normStdDev",
            "normMin",
            "normMax",
            "normRange",
        ]
    else:
        raise ValueError(
            f"Unknown column type {column_type}. Use either parcellation or segmentation."
        )

    try:
        return pd.read_csv(
            stats_filename,
            names=columns,
            comment="#",
            header=None,
            delimiter="\s+",
            dtype=str,
        )
    except FileNotFoundError:
        raise FileNotFoundError(f"Stats file {stats_filename} could not be found.")


def _get_stats_filename_for_atlas(stats_folder: PurePath, atlas: str) -> dict:
    """Gets the stats file names for each hemisphere for requested atlas.

    Parameters
    ----------
    stats_folder : PurePath
        Path to the statistics folder.

    atlas : str
        Name of the atlas to get statistics for.

    Returns
    -------
    dict :
        Dictionary mapping hemispheres to statistics file names.
        Indexes are 'lh' for 'left hemisphere', and 'rh' for 'right hemisphere'.
    """
    atlas_dict = {"desikan": "aparc", "destrieux": "aparc.a2009s", "ba": "BA_exvivo"}
    atlas_filename = atlas_dict.get(atlas, atlas)

    return {
        hemi: stats_folder / f"{hemi}.{atlas_filename}.stats" for hemi in ("lh", "rh")
    }


def _generate_tsv_for_parcellation(
    stats_folder: PurePath,
    output_dir: PurePath,
    prefix: str,
    atlases: List[str],
) -> None:
    """Generate TSV files for parcellation files.

    This function will generate .tsv from:
        - the table in the .stats file
        - the secondary (commented out) information common to
          both 'left' and 'right' .stats files
    """
    info_dict = {
        InfoType.VOLUME: "GrayVol",
        InfoType.THICKNESS: "ThickAvg",
        InfoType.AREA: "SurfArea",
        InfoType.MEANCURV: "MeanCurv",
    }
    for atlas in atlases:
        hemi_to_stats_filename = _get_stats_filename_for_atlas(stats_folder, atlas)
        hemi_to_stats_df = {
            hemi: _read_stats_file(filename, ColumnType.PARCELLATION)
            for hemi, filename in hemi_to_stats_filename.items()
        }
        for info, col_name in info_dict.items():
            # Secondary information are common to left and right hemispheres
            stats_df = get_secondary_stats(hemi_to_stats_filename["lh"], info)
            for hemi in ("lh", "rh"):
                stats_df = pd.concat(
                    [
                        stats_df,
                        hemi_to_stats_df[hemi][["StructName", col_name]].rename(
                            columns={
                                "StructName": "label_name",
                                col_name: "label_value",
                            }
                        ),
                    ]
                )
            output_filename = (
                output_dir / f"{prefix}_parcellation-{atlas}_{str(info)}.tsv"
            )
            stats_df.to_csv(output_filename, sep="\t", index=False, encoding="utf-8")


def _generate_tsv_for_segmentation(
    stats_folder: PurePath,
    output_dir: PurePath,
    prefix: str,
) -> None:
    """Generate TSV files for segmentation."""
    for filename, suffix in zip(
        ["aseg.stats", "wmparc.stats"],
        ["segmentationVolumes.tsv", "parcellation-wm_volume.tsv"],
    ):
        stats_filename = stats_folder / filename
        primary_stats_df = _read_stats_file(stats_filename, ColumnType.SEGMENTATION)
        primary_stats_df.rename(
            columns={"StructName": "label_name", "Volume_mm3": "label_value"},
            inplace=True,
        )
        secondary_stats_df = get_secondary_stats(stats_filename, InfoType.VOLUME)
        full_stats_df = pd.concat(
            [
                primary_stats_df[["label_name", "label_value"]],
                secondary_stats_df,
            ]
        )
        full_stats_df.to_csv(
            output_dir / f"{prefix}_{suffix}", sep="\t", index=False, encoding="utf-8"
        )


def generate_regional_measures(
    segmentation_path: str,
    subject_id: str,
    atlases: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
) -> None:
    """Read stats files and generate regional measures TSV files.

    The stats files are located in
    <segmentation_path>/<subject_id>/stats/*.stats
    and the generated TSV files will be located in the folder
    <segmentation_path>/regional_measures

    Notes
    -----
    The .stats files contain:
        - a table with statistical information (e.g., structure volume)
        - secondary statistical information with all lines starting
          with the sentence '# Measure'.

    The .tsv files return the relevant statistical information from both sources.

    Parameters
    ----------
    segmentation_path : str
        Path to the FreeSurfer segmentation.

    subject_id : str
        Subject ID in the form `sub-CLNC01_ses-M00`, `sub-CLNC01_long-M00M18` or
        `sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18`.

    atlases : List[str], optional
        The atlas names to be used.
        By default, the following atlases will be used:
            - "desikan"
            - "destrieux"
            - "ba"

    output_dir : str, optional
        Path to the folder where the .tsv stats files will be stored.
        If no output folder is provided, files will be written to
        [path_segmentation]/regional_measures.
    """
    import os
    from pathlib import Path

    atlases = atlases or ["desikan", "destrieux", "ba"]
    if not isinstance(atlases, list):
        raise ValueError(
            f"atlases should be a list of strings. {type(atlases)} was provided instead."
        )
    prefix = _get_prefix(subject_id)
    seg_path = Path(os.path.expanduser(segmentation_path))
    stats_folder = seg_path / subject_id / "stats"

    if not stats_folder.is_dir():
        raise FileNotFoundError(
            f"Image {prefix.replace('_', ' | ')} does not contain FreeSurfer segmentation"
        )
    output_dir = Path(output_dir) if output_dir else seg_path / "regional_measures"
    os.makedirs(output_dir, exist_ok=True)
    _generate_tsv_for_parcellation(stats_folder, output_dir, prefix, atlases)
    _generate_tsv_for_segmentation(stats_folder, output_dir, prefix)


def extract_image_id_from_longitudinal_segmentation(freesurfer_id: str):
    """Extract image ID from longitudinal segmentation folder.

    This function will extract participant, session and longitudinal ID from `freesurfer_id`.

    Parameters
    ----------
    freesurfer_id : str
        The freesurfer ID.

    Returns
    -------
    image_id: NamedTuple
        NamedTuple containing the following fields:
            - "participant_id": The extracted participant ID.
            - "session_id": The extracted session ID.
              Note that this field can be an empty string.
            - "long_id": The extracted longitudinal ID.
              Note that this field can be an empty string.

    Examples
    --------
    >>> from clinica.utils.freesurfer import extract_image_id_from_longitudinal_segmentation
    >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC001_ses-M000')
    image_id(participant_id='sub-CLNC001', session_id='ses-M000', long_id='')
    >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC001_long-M018')
    image_id(participant_id='sub-CLNC001', session_id='', long_id='long-M018')
    >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC001_ses-M000.long.sub-CLNC001_long-M000M018')
    image_id(participant_id='sub-CLNC001', session_id='ses-M000', long_id='long-M000M018')
    """
    from collections import namedtuple

    image_id = namedtuple("image_id", ["participant_id", "session_id", "long_id"])

    # Case 'sub-CLNC001_ses-M000.long.sub-CLNC001_long-M000M018'
    if ".long." in freesurfer_id:
        participant_id = freesurfer_id.split(".long.")[0].split("_")[0]
        session_id = freesurfer_id.split(".long.")[0].split("_")[1]
        long_id = freesurfer_id.split(".long.")[1].split("_")[1]
    # Case 'sub-CLNC001_long-M000M018'
    elif "long-" in freesurfer_id:
        participant_id = freesurfer_id.split("_")[0]
        session_id = ""
        long_id = freesurfer_id.split("_")[1]
    # Case 'sub-CLNC001_ses-M000'
    else:
        try:
            participant_id = freesurfer_id.split("_")[0]
            session_id = freesurfer_id.split("_")[1]
            long_id = ""
        except IndexError:
            raise ValueError(
                f"The provided Freesurfer ID {freesurfer_id} could not be parsed."
            )

    return image_id(participant_id, session_id, long_id)


class InfoType(Enum):
    VOLUME = auto()
    THICKNESS = auto()
    AREA = auto()
    MEANCURV = auto()

    def __str__(self):
        return self._name_.lower()


def _stats_line_filter(line: str, info: InfoType) -> bool:
    """Helper function to filter lines from stats files.

    The function returns True if the line starts with '# Measure',
    and ends with a marker specific to the type of info provided.
    More precisely, the marker is:

        - `mm^3` if info type is `volume`.
        - `mm` if  info type is `thickness`.
        - `mm^2` if info type is `area`.

    .. note::
        Info type `meancurv` is not supported.

    Parameters
    ----------
    line : str
        The line of the stat file to analyze.

    info : InfoType
        The type of information we are looking for.

    Returns
    -------
    bool :
        True if line should be kept, False otherwise.
    """
    end_line_marker = {
        InfoType.VOLUME: "mm^3",
        InfoType.THICKNESS: "mm",
        InfoType.AREA: "mm^2",
    }
    return line.startswith("# Measure") and line.endswith(end_line_marker[info])


def _extract_region_and_stat_value(
    line: str, info: InfoType
) -> Optional[Tuple[str, str]]:
    """Helper function to extract the region label and corresponding stat value
    from a stat line kept by `_stats_line_filter`.

    Parameters
    ----------
    line : str
        The line of the stat file to analyze.

    info : InfoType
        The type of information we are looking for.

    Returns
    -------
    region : str
        The region name for this line.

        .. note::
            This will be an empty string if the line does not contain
            specific keywords related to the type of information.

    value : str
        The corresponding statistics value.

        .. note::
            This will be an empty string if the line does not contain
            specific keywords related to the type of information.
    """
    info_to_keywords = {
        InfoType.VOLUME: ["volume", "Volume"],
        InfoType.AREA: ["area", "Area"],
        InfoType.THICKNESS: ["thickness", "Thickness"],
    }
    word_list = line.replace(",", "").split()
    if any(x in word_list for x in info_to_keywords[info]):
        return word_list[2], word_list[-2]
    return None


def get_secondary_stats(
    stats_filename: PurePath, info_type: InfoType
) -> Optional[pd.DataFrame]:
    """Read the 'secondary' statistical info from .stats file.

    Extract the information from .stats file that is commented out
    (lines starting with '# Measure' prefix) and does not appear in the
    table at the end of the document.

    Parameters
    ----------
    stats_filename : PurePath
        Path to the .stats file.

    info_type : InfoType
        Type of information to read from .stats file.

    Returns
    -------
    secondary_stats_dict : pd.DataFrame
        DataFrame containing the regions of the brain (in a column named 'label_name'),
        and the corresponding volume/thickness/area values depending on the
        input info type (in a column named 'label_value').
    """
    if info_type == InfoType.MEANCURV:
        return None

    with open(stats_filename, "r") as stats_file:
        stats = stats_file.read()
    stats_lines = [
        line for line in stats.splitlines() if _stats_line_filter(line, info_type)
    ]
    regions_stats = list(
        filter(
            None,
            (_extract_region_and_stat_value(line, info_type) for line in stats_lines),
        )
    )

    return pd.DataFrame.from_records(
        regions_stats, columns=["label_name", "label_value"]
    )


def check_flags(in_t1w: str, recon_all_args: str) -> str:
    """Check `recon_all_args` flags for `in_t1w` image.

    Currently, this function only adds '-cw256' if the FOV of `in_t1w` is greater than 256.

    Parameters
    ----------
    in_t1w : str
        Path to the T1W image.

    recon_all_args : str
        ??????.

    Returns
    -------
    str :
        ????.
    """
    import nibabel as nib

    f = nib.load(in_t1w)
    voxel_size = f.header.get_zooms()
    t1_size = f.header.get_data_shape()
    if any([v * t > 256 for v, t in zip(voxel_size, t1_size)]):
        return " ".join([recon_all_args, "-cw256"])
    return recon_all_args
