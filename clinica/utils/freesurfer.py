"""This module contains FreeSurfer utilities."""


def generate_regional_measures_alt(
    segmentation_path, subject_id, atlas, output_dir=None
):
    """
    Read stats files located in
    <segmentation_path>/<subject_id>/stats/*.stats
    and generate TSV files in <segmentation_path>/regional_measures
    folder.

    Note: the .stats files contain both 1) a table with statistical
    information (e.g., structure volume) and 2) 'secondary' statistical
    information with all lines starting with the sentence '# Measure'.
    The .tsv files return the relevant statistical information from both
    sources.

    Args:
        segmentation_path (string): Path to the FreeSurfer segmentation.
        subject_id (string): Subject ID in the form sub-CLNC01_ses-M00, sub-CLNC01_long-M00M18 or
            sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18
        output_dir (string): folder where the .tsv stats files will be
            stored. Will be [path_segmentation]/regional_measures if no
            dir is provided by the user
    """
    import errno
    import os

    import pandas

    from clinica.utils.freesurfer import (
        extract_image_id_from_longitudinal_segmentation,
        get_secondary_stats,
        write_tsv_file,
    )

    if segmentation_path != "":
        image_id = extract_image_id_from_longitudinal_segmentation(subject_id)
        prefix = image_id.participant_id
        if image_id.session_id:
            prefix = prefix + "_" + image_id.session_id
        if image_id.long_id:
            prefix = prefix + "_" + image_id.long_id

        stats_folder = os.path.join(
            os.path.expanduser(segmentation_path), subject_id, "stats"
        )
        if not os.path.isdir(stats_folder):
            raise IOError(
                "Image %s does not contain FreeSurfer segmentation"
                % prefix.replace("_", " | ")
            )

        if not output_dir:
            output_dir = os.path.join(segmentation_path, "regional_measures")
        try:
            os.makedirs(output_dir)
        except OSError as exception:
            # if dest_dir exists, go on, if its other error, raise
            if exception.errno != errno.EEXIST:
                raise

        # Generate TSV files for parcellation files
        #
        # Columns in ?h.BA.stats, ?h.aparc.stats or ?h.aparc.a2009s.stats file
        columns_parcellation = [
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
        hemi_dict = {"left": "lh", "right": "rh"}
        info_dict = {
            "volume": "GrayVol",
            "thickness": "ThickAvg",
            "area": "SurfArea",
            "meancurv": "MeanCurv",
        }
        for atlas in [atlas]:
            stats_filename_dict = dict()
            df_dict = dict()
            # read both left and right .stats files
            for hemi in ("left", "right"):
                stats_filename_dict[hemi] = os.path.join(
                    stats_folder,
                    "{0}.{1}.stats".format(hemi_dict[hemi], atlas),
                )
                df_dict[hemi] = pandas.read_csv(
                    stats_filename_dict[hemi],
                    names=columns_parcellation,
                    comment="#",
                    header=None,
                    delimiter="\s+",
                    dtype=str,
                )
            # generate .tsv from 1) the table in .stats file and 2) the
            # secondary (commented out) information common to both 'left'
            # and 'right' .stats file
            for info in ("volume", "thickness", "area", "meancurv"):
                # Secondary information (common to 'left' and 'right')
                secondary_stats_dict = get_secondary_stats(
                    stats_filename_dict["left"], info
                )
                # Join primary and secondary information
                key_list = (
                    list("lh_" + df_dict["left"]["StructName"])
                    + list("rh_" + df_dict["right"]["StructName"])
                    + list(secondary_stats_dict.keys())
                )
                col_name = info_dict[info]
                value_list = (
                    list(df_dict["left"][col_name])
                    + list(df_dict["right"][col_name])
                    + list(secondary_stats_dict.values())
                )
                # Write .tsv
                write_tsv_file(
                    os.path.join(
                        output_dir,
                        "{0}_parcellation-{1}_{2}.tsv".format(prefix, atlas, info),
                    ),
                    key_list,
                    info,
                    value_list,
                )

        # Generate TSV files for segmentation files
        #
        # Columns in aseg.stats or wmparc.stats file
        columns_segmentation = [
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

        # Parsing aseg.stats
        stats_filename = os.path.join(stats_folder, "aseg.stats")
        df = pandas.read_csv(
            stats_filename,
            comment="#",
            header=None,
            delimiter="\s+",
            dtype=str,
            names=columns_segmentation,
        )
        secondary_stats_dict = get_secondary_stats(stats_filename, "volume")
        key_list = list(df["StructName"]) + list(secondary_stats_dict.keys())
        value_list = list(df["Volume_mm3"]) + list(secondary_stats_dict.values())
        write_tsv_file(
            os.path.join(output_dir, prefix + "_segmentationVolumes.tsv"),
            key_list,
            "volume",
            value_list,
        )

        # Parsing wmparc.stats
        stats_filename = os.path.join(stats_folder, "wmparc.stats")
        df = pandas.read_csv(
            stats_filename,
            comment="#",
            header=None,
            delimiter="\s+",
            dtype=str,
            names=columns_segmentation,
        )
        secondary_stats_dict = get_secondary_stats(stats_filename, "volume")
        key_list = list(df["StructName"]) + list(secondary_stats_dict.keys())
        value_list = list(df["Volume_mm3"]) + list(secondary_stats_dict.values())
        write_tsv_file(
            os.path.join(output_dir, prefix + "_parcellation-wm_volume.tsv"),
            key_list,
            "volume",
            value_list,
        )


def extract_image_id_from_longitudinal_segmentation(freesurfer_id):
    """Extract image ID from longitudinal segmentation folder.

    This function will extract participant, session and longitudinal ID from `freesurfer_id`.

    Example:
        >>> from clinica.utils.freesurfer import extract_image_id_from_longitudinal_segmentation
        >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC01_ses-M00')
        image_id(participant_id='sub-CLNC01', session_id='ses-M00', long_id='')
        >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC01_long-M0018')
        image_id(participant_id='sub-CLNC01', session_id='', long_id='long-M0018')
        >>> extract_image_id_from_longitudinal_segmentation('sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18')
        image_id(participant_id='sub-CLNC01', session_id='ses-M00', long_id='long-M00M18')
    """
    from collections import namedtuple

    image_id = namedtuple("image_id", ["participant_id", "session_id", "long_id"])

    # Case 'sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18'
    if ".long." in freesurfer_id:
        participant_id = freesurfer_id.split(".long.")[0].split("_")[0]
        session_id = freesurfer_id.split(".long.")[0].split("_")[1]
        long_id = freesurfer_id.split(".long.")[1].split("_")[1]
    # Case 'sub-CLNC01_long-M00M18'
    elif "long-" in freesurfer_id:
        participant_id = freesurfer_id.split("_")[0]
        session_id = ""
        long_id = freesurfer_id.split("_")[1]
    # Case 'sub-CLNC01_ses-M00'
    else:
        participant_id = freesurfer_id.split("_")[0]
        session_id = freesurfer_id.split("_")[1]
        long_id = ""

    return image_id(participant_id, session_id, long_id)


def get_secondary_stats(stats_filename, info_type):
    """Read the 'secondary' statistical info from .stats file.

    Extract the information from .stats file that is commented out
    (lines starting with '# Measure' prefix) and does not appear in the
    table at the end of the document.

    Args:
        stats_filename (string): path to the .stats file
        info_type (string): 'volume', 'thickness' or 'area'

    Returns:
        secondary_stats_dict (dictionary of float): keys are regions of
            the brain, associated values are the corresponding
            volume/thickness/area depending on the input info type
    """
    # initialise structure containing the secondary statistical info
    secondary_stats_dict = dict()

    # currently no additional information is provided by .stats file for
    # the mean curvature
    if info_type != "meancurv":
        # define how lines are supposed to end in the stats file, depending
        # on the type of information that is searched for
        endline_dict = dict()
        endline_dict["volume"] = "mm^3"
        endline_dict["thickness"] = "mm"
        endline_dict["area"] = "mm^2"

        # define keywords that are supposed to appear in commented lines
        # containing statistical information
        info_keyword_dict = dict()
        info_keyword_dict["volume"] = ["volume", "Volume"]
        info_keyword_dict["area"] = ["area", "Area"]
        info_keyword_dict["thickness"] = ["thickness", "Thickness"]

        # read stats file line by line and only keep relevant lines
        with open(stats_filename, "r") as stats_file:
            stats = stats_file.read()
        stats_line_list = stats.splitlines()
        for stats_line in stats_line_list:
            startswith_condition = stats_line.startswith("# Measure")
            endswith_condition = stats_line.endswith(endline_dict[info_type])
            if (startswith_condition) and (endswith_condition):
                stats_line_word_list = stats_line.replace(",", "").split()
                # sanity check: make sure any sensible variation of
                # 'volume', 'thickness' or 'area' appears inside the line
                if any(x in stats_line_word_list for x in info_keyword_dict[info_type]):
                    # add info
                    info_region = stats_line_word_list[2]
                    info_value = stats_line_word_list[-2]
                    secondary_stats_dict[info_region] = info_value

    return secondary_stats_dict


def generate_regional_measures(segmentation_path, subject_id, output_dir=None):
    """General regional measures files from FreeSurfer segmentation.

    Read stats files located in <segmentation_path>/<subject_id>/stats/*.stats
    and generate TSV files in <segmentation_path>/regional_measures folder.

    Note: the .stats files contain both 1) a table with statistical
    information (e.g., structure volume) and 2) 'secondary' statistical
    information with all lines starting with the sentence '# Measure'.
    The .tsv files return the relevant statistical information from both
    sources.

    Args:
        segmentation_path (string): Path to the FreeSurfer segmentation.
        subject_id (string): Subject ID in the form sub-CLNC01_ses-M00, sub-CLNC01_long-M00M18 or
            sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18
        output_dir (string): folder where the .tsv stats files will be
            stored. Will be [path_segmentation]/regional_measures if no
            dir is provided by the user
    """
    import os

    import pandas

    from clinica.utils.freesurfer import write_tsv_file

    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)
    prefix = image_id.participant_id
    if image_id.session_id:
        prefix = f"{prefix}_{image_id.session_id}"
    if image_id.long_id:
        prefix = f"{prefix}_{image_id.long_id}"

    stats_folder = os.path.join(
        os.path.expanduser(segmentation_path), subject_id, "stats"
    )

    if not os.path.isdir(stats_folder):
        raise IOError(
            f"Image {prefix.replace('_', ' | ')} does not contain FreeSurfer segmentation"
        )

    if not output_dir:
        output_dir = os.path.join(segmentation_path, "regional_measures")
    os.makedirs(output_dir, exist_ok=True)

    # Generate TSV files for parcellation files
    #
    # Columns in ?h.BA.stats, ?h.aparc.stats or ?h.aparc.a2009s.stats file
    columns_parcellation = [
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
    hemi_dict = {"left": "lh", "right": "rh"}
    atlas_dict = {"desikan": "aparc", "destrieux": "aparc.a2009s", "ba": "BA_exvivo"}
    info_dict = {
        "volume": "GrayVol",
        "thickness": "ThickAvg",
        "area": "SurfArea",
        "meancurv": "MeanCurv",
    }
    for atlas in ("desikan", "destrieux", "ba"):
        stats_filename_dict = dict()
        df_dict = dict()
        # read both left and right .stats files
        for hemi in ("left", "right"):
            stats_filename_dict[hemi] = os.path.join(
                stats_folder, "{0}.{1}.stats".format(hemi_dict[hemi], atlas_dict[atlas])
            )
            df_dict[hemi] = pandas.read_csv(
                stats_filename_dict[hemi],
                names=columns_parcellation,
                comment="#",
                header=None,
                delimiter="\s+",
                dtype=str,
            )
        # generate .tsv from 1) the table in .stats file and 2) the
        # secondary (commented out) information common to both 'left'
        # and 'right' .stats file
        for info in ("volume", "thickness", "area", "meancurv"):
            # Secondary information (common to 'left' and 'right')
            secondary_stats_dict = get_secondary_stats(
                stats_filename_dict["left"], info
            )
            # Join primary and secondary information
            key_list = (
                list("lh_" + df_dict["left"]["StructName"])
                + list("rh_" + df_dict["right"]["StructName"])
                + list(secondary_stats_dict.keys())
            )
            col_name = info_dict[info]
            value_list = (
                list(df_dict["left"][col_name])
                + list(df_dict["right"][col_name])
                + list(secondary_stats_dict.values())
            )
            # Write .tsv
            write_tsv_file(
                os.path.join(
                    output_dir,
                    f"{prefix}_parcellation-{atlas}_{info}.tsv",
                ),
                key_list,
                info,
                value_list,
            )

    # Generate TSV files for segmentation files
    #
    # Columns in aseg.stats or wmparc.stats file
    columns_segmentation = [
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

    # Parsing aseg.stats
    stats_filename = os.path.join(stats_folder, "aseg.stats")
    df = pandas.read_csv(
        stats_filename,
        comment="#",
        header=None,
        delimiter="\s+",
        dtype=str,
        names=columns_segmentation,
    )
    secondary_stats_dict = get_secondary_stats(stats_filename, "volume")
    key_list = list(df["StructName"]) + list(secondary_stats_dict.keys())
    value_list = list(df["Volume_mm3"]) + list(secondary_stats_dict.values())
    write_tsv_file(
        os.path.join(output_dir, prefix + "_segmentationVolumes.tsv"),
        key_list,
        "volume",
        value_list,
    )

    # Parsing wmparc.stats
    stats_filename = os.path.join(stats_folder, "wmparc.stats")
    df = pandas.read_csv(
        stats_filename,
        comment="#",
        header=None,
        delimiter="\s+",
        dtype=str,
        names=columns_segmentation,
    )
    secondary_stats_dict = get_secondary_stats(stats_filename, "volume")
    key_list = list(df["StructName"]) + list(secondary_stats_dict.keys())
    value_list = list(df["Volume_mm3"]) + list(secondary_stats_dict.values())
    write_tsv_file(
        os.path.join(output_dir, prefix + "_parcellation-wm_volume.tsv"),
        key_list,
        "volume",
        value_list,
    )


def write_tsv_file(out_filename, name_list, scalar_name, scalar_list):
    """Write a .tsv file with list of keys and values.

    Args:
        out_filename (string): name of the .tsv file
        name_list (list of string): list of keys
        scalar_name (string): 'volume', 'thickness', 'area' or
            'meanCurv'. Not used for now. Might be used as part of a
            pandas data frame
        scalar_list (list of float): list of values corresponding to the keys
    """
    import warnings

    import pandas

    try:
        data = pandas.DataFrame({"label_name": name_list, "label_value": scalar_list})
        data.to_csv(out_filename, sep="\t", index=False, encoding="utf-8")
    except Exception as exception:
        warnings.warn("Impossible to save {0} file".format(out_filename))
        raise exception

    return out_filename


def check_flags(in_t1w, recon_all_args):
    """Check `recon_all_args` flags for `in_t1w` image.

    Currently, this function only adds '-cw256' if the FOV of `in_t1w` is greater than 256.
    """
    import nibabel as nib

    # from clinica.utils.stream import cprint

    f = nib.load(in_t1w)
    voxel_size = f.header.get_zooms()
    t1_size = f.header.get_data_shape()
    if (
        (voxel_size[0] * t1_size[0] > 256)
        or (voxel_size[1] * t1_size[1] > 256)
        or (voxel_size[2] * t1_size[2] > 256)
    ):
        # cprint(f"Setting MRI Convert to crop images to 256 FOV for {in_t1w} file.")
        optional_flag = " -cw256"
    else:
        # cprint(f"No need to add -cw256 flag for {in_t1w} file.")
        optional_flag = ""
    flags = "{0}".format(recon_all_args) + optional_flag

    return flags
