"""Module for converting fMRI of ADNI."""

from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

__all__ = ["convert_fmri"]


def convert_fmri(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    force_new_extraction: bool = False,
    n_procs: int = 1,
    convert_multiband: bool = True,
):
    """Convert fMR images of ADNI into BIDS format.

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

    force_new_extraction : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.

    convert_multiband : bool, optional
        If True, the converter will also convert multiband fmri scans.
        Otherwise, they will be ignored.
        Default=True.
    """
    from clinica.utils.stream import cprint

    from .._utils import ADNIModalityConverter, paths_to_bids

    cprint(
        f"Calculating paths of {ADNIModalityConverter.FMRI} images. Output will be stored in {conversion_dir}.",
        lvl="info",
    )
    images = _compute_fmri_path(
        source_dir, csv_dir, subjects, conversion_dir, convert_multiband
    )
    cprint(
        f"Paths of {ADNIModalityConverter.FMRI} images found. Exporting images into BIDS ...",
        lvl="info",
    )
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.FMRI,
        force_new_extraction=force_new_extraction,
        n_procs=n_procs,
    )
    cprint(msg=f"{ADNIModalityConverter.FMRI} conversion done.", lvl="debug")


def _compute_fmri_path(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
    convert_multiband: bool,
) -> pd.DataFrame:
    """Compute the paths to fMR images.

    Parameters
    ----------
    source_dir : Path
        The path to the ADNI directory.

    csv_dir : Path
        The path to the clinical data directory.

    subjects : list of str
        The list of subjects.

    conversion_dir : Path
        The path to the TSV files including the paths to original images.

    convert_multiband : bool, optional
        If True, the converter will also convert multiband fmri scans.
        Otherwise, they will be ignored.
        Default=True.

    Returns
    -------
    images : pd.DataFrame
        Pandas Dataframe containing the path for each fmri.
    """
    from .._utils import load_clinical_csv
    from ._image_path_utils import find_image_path
    from ._visits_utils import visits_to_timepoints

    fmri_dfs_list = []
    fmri_df = _initialize_fmri_df()
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

    # Selecting fMRI images
    mri_list = mri_list[mri_list.SEQUENCE.str.contains("MRI")]
    unwanted_sequences = [] if convert_multiband else ["MB"]
    mri_list = mri_list[
        mri_list.SEQUENCE.map(
            lambda x: not any(subs in x for subs in unwanted_sequences)
        )
    ]
    # We will convert the images for each subject in the subject list
    for subject in subjects:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subject]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")
        mri_list_subj = mri_list[mri_list.SUBJECT == subject]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")
        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subject[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(
            subject, mri_list_subj, adnimerge_subj, modality="fMRI"
        )

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]
            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            selected_image = _select_fmri_image(
                subject, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )
            if selected_image is not None:
                fmri_dfs_list.append(pd.DataFrame(selected_image, index=["i"]))
    if fmri_dfs_list:
        fmri_df = pd.concat(fmri_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not fmri_df.empty:
        error_ind = fmri_df.index[
            fmri_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        fmri_df.drop(error_ind, inplace=True)
    # Checking for images paths in filesystem
    images = find_image_path(fmri_df, source_dir, modality="fMRI")
    images.to_csv(conversion_dir / "fmri_paths.tsv", sep="\t", index=False)

    return images


def _initialize_fmri_df() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
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
    )


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [
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


def _select_fmri_image(
    subject_id: str,
    timepoint: str,
    visit_str: str,
    visit_mri_list: pd.DataFrame,
    mri_qc_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select the fMRI image for the given subject and visit based on QC.

    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary.

    Parameters
    ----------
    subject_id : str
        The subject identifier.

    timepoint : str
        The visit code.

    visit_str : str
        The visit name.

    visit_mri_list : pd.DataFrame
        List of images metadata.

    mri_qc_subj : pd.DataFrame
        A Dataframe containing list of QC of scans for the subject.

    Returns
    -------
    dict or None :
        A dictionary which contains image metadata.
        If None, no image was found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol

    from ._qc_utils import select_image_qc

    mri_qc_subj.columns = [x.lower() for x in mri_qc_subj.columns]
    selected_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if not selected_image:
        return None
    selected_scan = visit_mri_list[visit_mri_list.IMAGEUID == selected_image].iloc[0]

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_special_characters_with_symbol(
            selected_scan.SEQUENCE, symbol="_"
        ),
        "Scan_Date": selected_scan["SCANDATE"],
        "Study_ID": str(int(selected_scan.STUDYID)),
        "Series_ID": str(int(selected_scan.SERIESID)),
        "Image_ID": str(int(selected_scan.IMAGEUID)),
        "Field_Strength": selected_scan.MAGSTRENGTH,
    }
