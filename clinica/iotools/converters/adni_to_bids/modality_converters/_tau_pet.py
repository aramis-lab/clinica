"""Module for converting Tau PET of ADNI."""

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

__all__ = ["convert_tau_pet"]


def convert_tau_pet(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    mod_to_update: bool = False,
    n_procs: int = 1,
):
    """Convert Tau PET images of ADNI into BIDS format.

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
        Default=1
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        ADNIModalityConverter,
        load_clinical_csv,
        paths_to_bids,
    )
    from clinica.utils.stream import cprint

    cprint(
        (
            f"Calculating paths of {ADNIModalityConverter.PET_TAU.value} images. "
            f"Output will be stored in {conversion_dir}."
        ),
        lvl="info",
    )
    images = _compute_tau_pet_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint(
        f"Paths of {ADNIModalityConverter.PET_TAU.value} images found. Exporting images into BIDS ...",
        lvl="info",
    )
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.PET_TAU,
        mod_to_update=mod_to_update,
        n_procs=n_procs,
    )
    cprint(msg=f"{ADNIModalityConverter.PET_TAU.value} conversion done.", lvl="debug")


def _compute_tau_pet_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute the paths to Tau PET images.

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
    pd.DataFrame :
        A pandas Dataframe containing the path for each Tau PET image.
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv
    from clinica.utils.pet import Tracer

    from ._image_path_utils import find_image_path
    from ._pet_utils import get_images_pet

    pet_tau_df = pd.DataFrame(columns=_get_tau_pet_df_columns())
    pet_tau_dfs_list = []
    tauqc = load_clinical_csv(csv_dir, "TAUQC")
    tauqc3 = load_clinical_csv(csv_dir, "TAUQC3")
    pet_meta_list = load_clinical_csv(csv_dir, "PET_META_LIST")

    for subject in subjects:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subject]
        if subject_pet_meta.empty:
            continue
        # QC for TAU PET images for ADNI 2
        tau_qc2_subj = tauqc[(tauqc.SCANQLTY == 1) & (tauqc.RID == int(subject[-4:]))]
        # QC for TAU PET images for ADNI 3
        tau_qc3_subj = tauqc3[
            (tauqc3.SCANQLTY == 1) & (tauqc3.RID == int(subject[-4:]))
        ]
        # Concatenating visits in both QC files
        tau_qc_subj = pd.concat(
            [tau_qc2_subj, tau_qc3_subj], axis=0, ignore_index=True, sort=False
        )
        tau_qc_subj.rename(columns={"SCANDATE": "EXAMDATE"}, inplace=True)
        subj_dfs_list = get_images_pet(
            subject=subject,
            pet_qc_subj=tau_qc_subj,
            subject_pet_meta=subject_pet_meta,
            df_cols=_get_tau_pet_df_columns(),
            modality="TAU-PET",
            sequences_preprocessing_step=["AV1451 Co-registered, Averaged"],
        )
        if subj_dfs_list:
            pet_tau_dfs_list += subj_dfs_list
    if pet_tau_dfs_list:
        pet_tau_df = pd.concat(pet_tau_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not pet_tau_df.empty:
        error_ind = pet_tau_df.index[
            pet_tau_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        pet_tau_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(pet_tau_df, source_dir, "TAU")
    images.to_csv(
        conversion_dir / f"{Tracer.AV1451.value}_pet_paths.tsv",
        sep="\t",
        index=False,
    )

    return images


def _get_tau_pet_df_columns() -> list[str]:
    return [
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
    ]


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [("098_S_4275", "m84")]
