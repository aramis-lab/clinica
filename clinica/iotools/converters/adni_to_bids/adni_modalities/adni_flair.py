"""Module for converting FLAIR of ADNI."""
from os import PathLike
from typing import List, Optional


def convert_adni_flair(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: Optional[List[str]] = None,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert FLAIR images of ADNI into BIDS format.

    Parameters
    ----------
    source_dir : PathLike
        Path to the ADNI directory.

    csv_dir : PathLike
        Path to the clinical data directory.

    destination_dir : PathLike
        Path to the destination BIDS directory.

    conversion_dir : PathLike
        Path to the TSV files including the paths to original images.

    subjects : List of str, optional
        List of subjects.

    mod_to_update : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        load_clinical_csv,
        paths_to_bids,
    )
    from clinica.utils.stream import cprint

    if not subjects:
        adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")
        subjects = list(adni_merge.PTID.unique())

    cprint(
        f"Calculating paths of FLAIR images. Output will be stored in {conversion_dir}."
    )
    images = compute_flair_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of FLAIR images found. Exporting images into BIDS ...")
    # flair_paths_to_bids(images, dest_dir)
    paths_to_bids(
        images, destination_dir, "flair", mod_to_update=mod_to_update, n_procs=n_procs
    )
    cprint(msg="FLAIR conversion done.", lvl="debug")


def compute_flair_paths(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute the paths to the FLAIR images and store them in a TSV file.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns:
        images: a dataframe with all the paths to the FLAIR images that will be converted into BIDS
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        load_clinical_csv,
        visits_to_timepoints,
    )

    flair_col_df = [
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Series_ID",
        "Image_ID",
        "Field_Strength",
        "Scanner",
    ]
    flair_df = pd.DataFrame(columns=flair_col_df)
    flair_dfs_list = []

    # Loading needed .csv files
    adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")
    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_12_08_15")
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "AFL"]
    mri_list = load_clinical_csv(csv_dir, "MRILIST")

    # Selecting FLAIR DTI images that are not MPR
    mri_list = mri_list[mri_list.SEQUENCE.str.contains("flair", case=False, na=False)]
    unwanted_sequences = ["_MPR_"]
    mri_list = mri_list[
        mri_list.SEQUENCE.map(
            lambda x: not any(subs in x for subs in unwanted_sequences)
        )
    ]

    for subj in subjs_list:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subj, mri_list_subj, adnimerge_subj, "FLAIR")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            flair = flair_image(
                subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )

            if flair is not None:
                row_to_append = pd.DataFrame(
                    flair,
                    index=[
                        "i",
                    ],
                )
                flair_dfs_list.append(row_to_append)

    if flair_dfs_list:
        flair_df = pd.concat(flair_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # Eq_1 images
        ("141_S_0767", "m84"),
        ("067_S_5205", "bl"),
        ("127_S_4928", "m24"),
        ("024_S_4674", "m06"),
        ("123_S_2363", "m24"),
        ("053_S_4578", "m48"),
        ("128_S_4586", "m48"),
        ("053_S_4813", "m48"),
        ("053_S_5272", "m24"),
        ("013_S_1186", "m48"),
        ("031_S_2022", "bl"),
        ("031_S_2022", "m06"),
        ("031_S_2233", "bl"),
        ("031_S_2233", "m03"),
        ("029_S_2395", "m72"),
        ("130_S_6043", "bl"),
        ("031_S_2018", "bl"),
        ("027_S_5170", "m72"),
        ("135_S_6284", "m12"),
        ("068_S_0127", "m180"),
        ("068_S_2187", "m120"),
        # Several output images
        ("114_S_6039", "bl"),
    ]

    # Removing known exceptions from images to convert
    if not flair_df.empty:
        error_ind = flair_df.index[
            flair_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        flair_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(flair_df, source_dir, "FLAIR", "S", "Series_ID")
    images.to_csv(path.join(conversion_dir, "flair_paths.tsv"), sep="\t", index=False)

    return images


def flair_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
    """
    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary

    Args:
        subject_id: Subject identifier
        timepoint: Visit code
        visit_str: Visit name
        visit_mri_list: List of images metadata
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: dictionary - contains image metadata
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
        select_image_qc,
    )

    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if sel_image is None:
        return None

    sel_scan = visit_mri_list[visit_mri_list.IMAGEUID == sel_image].iloc[0]

    image_dict = {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_sequence_chars(sel_scan.SEQUENCE),
        "Scan_Date": sel_scan["SCANDATE"],
        "Study_ID": str(int(sel_scan.STUDYID)),
        "Series_ID": str(int(sel_scan.SERIESID)),
        "Image_ID": str(int(sel_scan.IMAGEUID)),
        "Field_Strength": sel_scan.MAGSTRENGTH,
    }

    return image_dict
