"""Module for converting AV45 and Florbetaben PET of ADNI."""

from pathlib import Path
from typing import Iterable

import pandas as pd

__all__ = ["convert_av45_fbb_pet"]


def convert_av45_fbb_pet(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    mod_to_update: bool = False,
    n_procs: int = 1,
):
    """Convert AV-45 and Florbetaben PET images of ADNI into BIDS format.

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

    n_procs : int, default=1
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
        f"Calculating paths of AV45 and Florbetaben PET images. Output will be stored in {conversion_dir}."
    )
    images = _compute_av45_fbb_pet_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint(
        "Paths of AV45 and Florbetaben PET images found. Exporting images into BIDS ..."
    )
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.PET_AV45,
        mod_to_update=mod_to_update,
        n_procs=n_procs,
    )
    cprint(msg="AV45 and Florbetaben PET conversion done.", lvl="debug")


def _compute_av45_fbb_pet_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute the paths to the AV45 and Florbetaben PET images and store them in a TSV file.

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
        Dataframe with all the paths to the PET images that will be converted into BIDS.
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    from ._image_path_utils import find_image_path
    from ._pet_utils import get_images_pet

    pet_amyloid_col = [
        "Phase",
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Series_ID",
        "Image_ID",
        "Original",
        "Tracer",
    ]
    pet_amyloid_df = pd.DataFrame(columns=pet_amyloid_col)
    pet_amyloid_dfs_list = []

    # Loading needed .csv files
    av45qc = load_clinical_csv(csv_dir, "AV45QC")
    amyqc = load_clinical_csv(csv_dir, "AMYQC")
    pet_meta_list = load_clinical_csv(csv_dir, "PET_META_LIST")

    for subject in subjects:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subject]

        if subject_pet_meta.empty:
            continue
        # QC for AV45 PET images for ADNI 1, GO and 2
        av45_qc_subj = av45qc[(av45qc.PASS == 1) & (av45qc.RID == int(subject[-4:]))]

        # QC for Amyloid PET images for ADNI 3
        amy_qc_subj = amyqc[(amyqc.SCANQLTY == 1) & (amyqc.RID == int(subject[-4:]))]
        amy_qc_subj.insert(0, "EXAMDATE", amy_qc_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        amyloid_qc_subj = pd.concat(
            [av45_qc_subj, amy_qc_subj], axis=0, ignore_index=True, sort=False
        )
        sequences_preprocessing_step = [
            "AV45 Co-registered, Averaged",
            "FBB Co-registered, Averaged",
        ]
        subj_dfs_list = get_images_pet(
            subject,
            amyloid_qc_subj,
            subject_pet_meta,
            pet_amyloid_col,
            "Amyloid-PET",
            sequences_preprocessing_step,
        )
        if subj_dfs_list:
            pet_amyloid_dfs_list += subj_dfs_list

    if pet_amyloid_dfs_list:
        pet_amyloid_df = pd.concat(pet_amyloid_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not pet_amyloid_df.empty:
        error_ind = pet_amyloid_df.index[
            pet_amyloid_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        pet_amyloid_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_amyloid_df, source_dir, "Amyloid")
    images.to_csv(conversion_dir / "amyloid_pet_paths.tsv", sep="\t", index=False)

    return images


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [
        ("128_S_2220", "m48"),
        # Several output images
        ("098_S_4275", "m84"),
    ]
