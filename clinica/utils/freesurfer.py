"""This module contains FreeSurfer utilities."""
from enum import Enum, auto
from pathlib import PurePath
from typing import List, Optional

import pandas as pd


class ColumnType(Enum):
    PARCELLATION = auto()  # The actual value is meaningless.
    SEGMENTATION = auto()


def _get_prefix(subject_id: str) -> str:
    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)

    return "_".join([i for i in image_id if i])


def _read_stats_file(stats_filename: PurePath, column_type: ColumnType) -> pd.DataFrame:
    """Read the provided statistics file and extract the relevant data."""
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

    return pd.read_csv(
        stats_filename,
        names=columns,
        comment="#",
        header=None,
        delimiter="\s+",
        dtype=str,
    )


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
    atlas_dict = {"desikan": "aparc", "destrieux": "aparc.a2009s", "ba": "BA_exvivo"}
    info_dict = {
        "volume": "GrayVol",
        "thickness": "ThickAvg",
        "area": "SurfArea",
        "meancurv": "MeanCurv",
    }
    for atlas in atlases:
        atlas_filename = atlas_dict.get(atlas, atlas)
        stats_filename_dict = {
            hemi: stats_folder / f"{hemi}.{atlas_filename}.stats"
            for hemi in ("lh", "rh")
        }
        df_dict = {
            hemi: _read_stats_file(stats_filename_dict[hemi], ColumnType.PARCELLATION)
            for hemi in ("lh", "rh")
        }
        for info in ("volume", "thickness", "area", "meancurv"):
            # Secondary information (common to 'left' and 'right')
            secondary_stats_dict = get_secondary_stats(stats_filename_dict["lh"], info)
            # Join primary and secondary information
            key_list = [
                hemi + "_" + df_dict[hemi]["StructName"] for hemi in ("lh", "rh")
            ]
            key_list += list(secondary_stats_dict.keys())
            col_name = info_dict[info]
            value_list = [df_dict[hemi][col_name] for hemi in ("lh", "rh")]
            value_list += list(secondary_stats_dict.values())
            output_filename = output_dir / f"{prefix}_parcellation-{atlas}_{info}.tsv"
            write_tsv_file(output_filename, key_list, value_list)


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
        df = _read_stats_file(stats_filename, ColumnType.SEGMENTATION)
        secondary_stats_dict = get_secondary_stats(stats_filename, "volume")
        key_list = list(df["StructName"]) + list(secondary_stats_dict.keys())
        value_list = list(df["Volume_mm3"]) + list(secondary_stats_dict.values())
        output_filename = output_dir / f"{prefix}_{suffix}"
        write_tsv_file(output_filename, key_list, value_list)


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
    import errno
    import os
    from pathlib import Path

    atlases = atlases or ["desikan", "destrieux", "ba"]
    prefix = _get_prefix(subject_id)
    seg_path = Path(os.path.expanduser(segmentation_path))
    stats_folder = seg_path / subject_id / "stats"

    if not stats_folder.is_dir():
        raise IOError(
            "Image %s does not contain FreeSurfer segmentation"
            % prefix.replace("_", " | ")
        )
    output_dir = Path(output_dir) if output_dir else seg_path / "regional_measures"
    try:
        os.makedirs(output_dir)
    except OSError as exception:
        # if dest_dir exists, go on, if its other error, raise
        if exception.errno != errno.EEXIST:
            raise
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


def get_secondary_stats(stats_filename: PurePath, info_type: InfoType) -> dict:
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
    secondary_stats_dict : dict
        The keys are regions of the brain, and the associated values
        are the corresponding volume/thickness/area depending on the
        input info type.
    """
    secondary_stats_dict = {}

    # currently no additional information is provided by .stats file for
    # the mean curvature
    if info_type == InfoType.MEANCURV:
        return {}

    # define how lines are supposed to end in the stats file, depending
    # on the type of information that is searched for
    endline_dict = {
        InfoType.VOLUME: "mm^3",
        InfoType.THICKNESS: "mm",
        InfoType.AREA: "mm^2",
    }

    # define keywords that are supposed to appear in commented lines
    # containing statistical information
    info_keyword_dict = {
        InfoType.VOLUME: ["volume", "Volume"],
        InfoType.AREA: ["area", "Area"],
        InfoType.THICKNESS: ["thickness", "Thickness"],
    }
    with open(stats_filename, "r") as stats_file:
        stats = stats_file.read()
    stats_line_list = stats.splitlines()
    for stats_line in stats_line_list:
        startswith_condition = stats_line.startswith("# Measure")
        endswith_condition = stats_line.endswith(endline_dict[info_type])
        if startswith_condition and endswith_condition:
            stats_line_word_list = stats_line.replace(",", "").split()
            # sanity check: make sure any sensible variation of
            # 'volume', 'thickness' or 'area' appears inside the line
            if any(x in stats_line_word_list for x in info_keyword_dict[info_type]):
                info_region = stats_line_word_list[2]
                info_value = stats_line_word_list[-2]
                secondary_stats_dict[info_region] = info_value

    return secondary_stats_dict


def write_tsv_file(
    out_filename: PurePath,
    name_list: List[str],
    scalar_list: List[float],
) -> None:
    """Write a .tsv file with list of keys and values.

    Parameters
    ----------
    out_filename : PurePath
        Path to the .tsv file to write to.
        Must have the same length as `name_list`.

    name_list : List[str]
        list of keys. Must have the same length as `out_filename`.

    scalar_list : List[float]
        list of values corresponding to the keys.
    """
    from clinica.utils.stream import cprint

    if len(name_list) != len(scalar_list):
        raise ValueError(
            "Cannot write to TSV because keys and values have different "
            f"lengths {len(name_list)} vs. {len(scalar_list)}."
        )
    try:
        data = pd.DataFrame({"label_name": name_list, "label_value": scalar_list})
        data.to_csv(out_filename, sep="\t", index=False, encoding="utf-8")
    except Exception as exception:
        cprint(msg=f"Impossible to save {out_filename} file.", lvl="warning")
        raise exception


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
