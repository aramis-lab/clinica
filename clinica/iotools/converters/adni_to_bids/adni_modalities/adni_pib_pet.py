"""Module for converting PIB PET of ADNI."""
from os import PathLike
from typing import List, Optional


def convert_adni_pib_pet(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: Optional[List[str]] = None,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert PIB PET images of ADNI into BIDS format.

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
        f"Calculating paths of PIB PET images. Output will be stored in {conversion_dir}."
    )
    images = compute_pib_pet_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of PIB PET images found. Exporting images into BIDS ...")
    paths_to_bids(
        images, destination_dir, "pib", mod_to_update=mod_to_update, n_procs=n_procs
    )
    cprint(msg="PIB PET conversion done.", lvl="debug")


def compute_pib_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute the paths to the PIB PET images and store them in a TSV file.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns:
        images: a dataframe with all the paths to the PET images that will be converted into BIDS
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        get_images_pet,
        load_clinical_csv,
    )
    from clinica.utils.pet import Tracer

    pet_pib_col = [
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
    pet_pib_df = pd.DataFrame(columns=pet_pib_col)
    pet_pib_dfs_list = []

    # Loading needed .csv files
    pibqc = load_clinical_csv(csv_dir, "PIBQC")
    pet_meta_list = load_clinical_csv(csv_dir, "PET_META_LIST")

    for subj in subjs_list:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for PIB PET images
        pet_qc_subj = pibqc[(pibqc.PASS == 1) & (pibqc.RID == int(subj[-4:]))]

        sequences_preprocessing_step = ["PIB Co-registered, Averaged"]
        subj_dfs_list = get_images_pet(
            subj,
            pet_qc_subj,
            subject_pet_meta,
            pet_pib_col,
            "PIB-PET",
            sequences_preprocessing_step,
            viscode_field="VISCODE",
        )
        if subj_dfs_list:
            pet_pib_dfs_list += subj_dfs_list

    if pet_pib_dfs_list:
        pet_pib_df = pd.concat(pet_pib_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = []

    # Removing known exceptions from images to convert
    if not pet_pib_df.empty:
        error_ind = pet_pib_df.index[
            pet_pib_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        pet_pib_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_pib_df, source_dir, "PIB", "I", "Image_ID")
    images.to_csv(
        path.join(conversion_dir, f"{Tracer.PIB}_pet_paths.tsv"), sep="\t", index=False
    )

    return images
