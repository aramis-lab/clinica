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
    force_new_extraction: bool = False,
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

    force_new_extraction : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.utils.stream import cprint

    from .._utils import ADNIModalityConverter, paths_to_bids

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
        force_new_extraction=force_new_extraction,
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
    from .._utils import load_clinical_csv
    from ._image_path_utils import find_image_path
    from ._mri_utils import load_all_images_metadata

    flair_df = _initialize_flair_df()
    flair_dfs_list = []
    # Loading needed .csv files
    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_05_07_15")
    all_images = load_all_images_metadata(csv_dir)

    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "AFL"]

    # Selecting FLAIR DTI images that are not MPR
    all_images = all_images[
        all_images.image_description.str.contains("FLAIR", case=False, na=False)
    ]
    all_images = all_images[
        all_images.image_description.map(
            lambda x: not any(subs in x for subs in ("MPR",))
        )
    ]
    for subject in subjects:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        all_images_subj = all_images[
            (all_images.subject_id == subject) & (all_images.modality == "MRI")
        ]
        all_images_subj = all_images_subj.sort_values("image_date")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subject[-4:])]

        for _, row in all_images_subj.iterrows():
            timepoint = row["image_visit"]
            visit_all_images = all_images_subj[all_images_subj.image_visit == timepoint]
            flair = _select_flair_image(
                subject, timepoint, visit_all_images, mayo_mri_qc_subj
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
    timepoint: str,
    visit_all_images: list,
    mri_qc_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select a FLAIR image based on QC information.

    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary

    Parameters
    ----------
    subject_id : str
        The subject identifier.

    timepoint : str
        The visit code.

    visit_all_images : list
        A list of images metadata.

    mri_qc_subj : pd.DataFrame
        A Dataframe containing list of QC of scans for the subject.

    Returns
    -------
    dict or None :
        Contains image metadata.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol

    from ._qc_utils import select_image_qc

    selected_image_id = select_image_qc(list(visit_all_images.image_id), mri_qc_subj)
    if selected_image_id is None:
        return None
    selected_scan = visit_all_images[
        visit_all_images.image_id == selected_image_id
    ].iloc[0]

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": replace_special_characters_with_symbol(
            selected_scan.image_description, symbol="_"
        ),
        "Scan_Date": selected_scan["image_date"],
        "Study_ID": str(int(selected_scan.study_id)),
        "Series_ID": str(int(selected_scan.series_id)),
        "Image_ID": str(int(selected_scan.image_id)),
        "Field_Strength": selected_scan.MagneticFieldStrength,
    }
