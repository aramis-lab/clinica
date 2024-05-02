"""Module for converting PIB PET of ADNI."""

from pathlib import Path
from typing import List, Optional

import pandas as pd

__all__ = ["convert_adni_pib_pet"]


def convert_adni_pib_pet(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: List[str],
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert PIB PET images of ADNI into BIDS format.

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
        load_clinical_csv,
        paths_to_bids,
    )
    from clinica.utils.stream import cprint

    cprint(
        f"Calculating paths of PIB PET images. Output will be stored in {conversion_dir}."
    )
    images = _compute_pib_pet_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of PIB PET images found. Exporting images into BIDS ...")
    paths_to_bids(
        images, destination_dir, "pib", mod_to_update=mod_to_update, n_procs=n_procs
    )
    cprint(msg="PIB PET conversion done.", lvl="debug")


def _compute_pib_pet_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: list[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute the paths to the PIB PET images and store them in a TSV file.

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
        A dataframe with all the paths to the PET images that will be converted into BIDS.
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        get_images_pet,
        load_clinical_csv,
    )
    from clinica.utils.pet import Tracer

    pet_pib_dfs_list = []
    pet_pib_df = pd.DataFrame(columns=_get_pib_pet_columns())
    pibqc = load_clinical_csv(csv_dir, "PIBQC")
    pet_meta_list = load_clinical_csv(csv_dir, "PET_META_LIST")

    for subject in subjects:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subject]
        if subject_pet_meta.empty:
            continue
        # QC for PIB PET images
        pet_qc_subj = pibqc[(pibqc.PASS == 1) & (pibqc.RID == int(subject[-4:]))]
        subj_dfs_list = get_images_pet(
            subject=subject,
            pet_qc_subj=pet_qc_subj,
            subject_pet_meta=subject_pet_meta,
            df_cols=_get_pib_pet_columns(),
            modality="PIB-PET",
            sequences_preprocessing_step=["PIB Co-registered, Averaged"],
            viscode_field="VISCODE",
        )
        if subj_dfs_list:
            pet_pib_dfs_list += subj_dfs_list
    if pet_pib_dfs_list:
        pet_pib_df = pd.concat(pet_pib_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not pet_pib_df.empty:
        error_ind = pet_pib_df.index[
            pet_pib_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        pet_pib_df.drop(error_ind, inplace=True)
    images = find_image_path(pet_pib_df, source_dir, modality="PIB")
    images.to_csv(
        conversion_dir / f"{Tracer.PIB.value}_pet_paths.tsv",
        sep="\t",
        index=False,
    )

    return images


def _get_pib_pet_columns() -> list[str]:
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


def _get_known_conversion_errors() -> List[tuple[str, str]]:
    return []
