"""Module for converting fMRI of ADNI."""
from os import PathLike
from typing import List, Optional


def convert_adni_fmri(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: Optional[List[str]] = None,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert fMR images of ADNI into BIDS format.

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
        f"Calculating paths of fMRI images. Output will be stored in {conversion_dir}."
    )
    images = compute_fmri_path(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of fMRI images found. Exporting images into BIDS ...")
    # fmri_paths_to_bids(dest_dir, images)
    paths_to_bids(
        images, destination_dir, "fmri", mod_to_update=mod_to_update, n_procs=n_procs
    )
    cprint(msg="fMRI conversion done.", lvl="debug")


def compute_fmri_path(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute the paths to fMR images.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns:
        pandas Dataframe containing the path for each fmri
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        load_clinical_csv,
        visits_to_timepoints,
    )

    fmri_col = [
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Field_Strength",
        "Series_ID",
        "Image_ID",
    ]
    fmri_df = pd.DataFrame(columns=fmri_col)
    fmri_dfs_list = []

    # Loading needed .csv files
    adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")

    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_12_08_15")

    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "fMRI"]
    mayo_mri_qc.columns = [x.upper() for x in mayo_mri_qc.columns]

    mayo_mri_qc3 = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_QUALITY_ADNI3")
    mayo_mri_qc3 = mayo_mri_qc3[mayo_mri_qc3.SERIES_TYPE == "EPB"]

    # Concatenating visits in both QC files
    mayo_mri_qc = pd.concat(
        [mayo_mri_qc, mayo_mri_qc3], axis=0, ignore_index=True, sort=False
    )
    mri_list = load_clinical_csv(csv_dir, "MRILIST")

    # Selecting only fMRI images that are not Multiband
    mri_list = mri_list[
        mri_list.SEQUENCE.str.contains("MRI")
    ]  # 'MRI' includes all fMRI and fMRI scans, but not others
    unwanted_sequences = ["MB"]
    mri_list = mri_list[
        mri_list.SEQUENCE.map(
            lambda x: not any(subs in x for subs in unwanted_sequences)
        )
    ]

    # We will convert the images for each subject in the subject list
    for subj in subjs_list:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subj, mri_list_subj, adnimerge_subj, "fMRI")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            image = fmri_image(
                subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )

            if image is not None:
                row_to_append = pd.DataFrame(
                    image,
                    index=[
                        "i",
                    ],
                )
                fmri_dfs_list.append(row_to_append)

    if fmri_dfs_list:
        fmri_df = pd.concat(fmri_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [
        ("006_S_4485", "m84"),
        ("123_S_4127", "m96"),
        # Eq_1
        ("094_S_4503", "m24"),
        ("009_S_4388", "m72"),
        ("036_S_6088", "bl"),
        ("036_S_6134", "bl"),
        ("016_S_6802", "bl"),
        ("016_S_6816", "bl"),
        ("126_S_4891", "m84"),
        # Multiple images
        ("029_S_2395", "m72"),
    ]

    # Removing known exceptions from images to convert
    if not fmri_df.empty:
        error_ind = fmri_df.index[
            fmri_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        fmri_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(fmri_df, source_dir, "fMRI", "S", "Series_ID")
    images.to_csv(path.join(conversion_dir, "fmri_paths.tsv"), sep="\t", index=False)

    return images


def fmri_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
    """
    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary.

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

    mri_qc_subj.columns = [x.lower() for x in mri_qc_subj.columns]
    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if not sel_image:
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
