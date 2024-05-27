"""Module for converting FLAIR of ADNI."""

from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

__all__ = ["convert_flair"]


def convert_flair(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    mod_to_update: bool = False,
    n_procs: int = 1,
):
    """Convert FLAIR images of ADNI into BIDS format.

    Parameters
    ----------
    source_dir : Path
        The path to the ADNI directory.

    csv_dir : Path
        The path to the clinical data directory.

    destination_dir : Path
        The path to the destination BIDS directory.

    conversion_dir : Path
        The path to the TSV files including the paths to original images.

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
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        ADNIModalityConverter,
        load_clinical_csv,
        paths_to_bids,
    )
    from clinica.utils.stream import cprint

    cprint(
        f"Calculating paths of {ADNIModalityConverter.FLAIR} images. Output will be stored in {conversion_dir}.",
        lvl="info",
    )
    images = _compute_flair_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint(
        f"Paths of {ADNIModalityConverter.FLAIR} images found. Exporting images into BIDS ...",
        lvl="info",
    )
    # flair_paths_to_bids(images, dest_dir)
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.FLAIR,
        mod_to_update=mod_to_update,
        n_procs=n_procs,
    )
    cprint(msg="FLAIR conversion done.", lvl="debug")


def _compute_flair_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute the paths to the FLAIR images and store them in a TSV file.

    Parameters
    ----------
    source_dir : Path
        The path to the ADNI directory.

    csv_dir : Path
        The path to the clinical data directory.

    subjects : list of str
        The subjects list.

    conversion_dir : Path
        The path to the TSV files including the paths to original images.

    Returns
    -------
    images : pd.DataFrame
        A dataframe with all the paths to the FLAIR images that will be converted into BIDS.
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    from ._image_path_utils import find_image_path
    from ._visits_utils import visits_to_timepoints

    flair_df = _initialize_flair_df()
    flair_dfs_list = []
    # Loading needed .csv files
    adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")
    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_12_08_15")
    mri_list = load_clinical_csv(csv_dir, "MRILIST")

    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "AFL"]

    # Selecting FLAIR DTI images that are not MPR
    mri_list = mri_list[mri_list.SEQUENCE.str.contains("flair", case=False, na=False)]
    mri_list = mri_list[
        mri_list.SEQUENCE.map(lambda x: not any(subs in x for subs in ("_MPR_",)))
    ]
    for subject in subjects:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subject]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")

        mri_list_subj = mri_list[mri_list.SUBJECT == subject]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subject[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subject, mri_list_subj, adnimerge_subj, "FLAIR")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]
            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            flair = _select_flair_image(
                subject, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )
            if flair is not None:
                row_to_append = pd.DataFrame(flair, index=["i"])
                flair_dfs_list.append(row_to_append)
    if flair_dfs_list:
        flair_df = pd.concat(flair_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not flair_df.empty:
        error_ind = flair_df.index[
            flair_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        flair_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(flair_df, source_dir, "FLAIR")
    images.to_csv(conversion_dir / "flair_paths.tsv", sep="\t", index=False)

    return images


def _initialize_flair_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
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
    )


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [
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


def _select_flair_image(
    subject_id: str,
    timepoint,
    visit_str: str,
    visit_mri_list: list,
    mri_qc_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select a FLAIR image based on QC information.

    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary

    Parameters
    ----------
    subject_id : str
        The subject identifier.

    timepoint : ??
        The visit code.

    visit_str : str
        The visit name.

    visit_mri_list : list
        A list of images metadata.

    mri_qc_subj : pd.DataFrame
        A Dataframe containing list of QC of scans for the subject.

    Returns
    -------
    dict or None :
        Contains image metadata.
    """
    from clinica.iotools.converter_utils import replace_sequence_chars

    from ._qc_utils import select_image_qc

    selected_image_id = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if selected_image_id is None:
        return None
    selected_scan = visit_mri_list[visit_mri_list.IMAGEUID == selected_image_id].iloc[0]

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_sequence_chars(selected_scan.SEQUENCE),
        "Scan_Date": selected_scan["SCANDATE"],
        "Study_ID": str(int(selected_scan.STUDYID)),
        "Series_ID": str(int(selected_scan.SERIESID)),
        "Image_ID": str(int(selected_scan.IMAGEUID)),
        "Field_Strength": selected_scan.MAGSTRENGTH,
    }
