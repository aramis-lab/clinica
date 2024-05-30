"""Module for converting Tau PET of ADNI."""
from os import PathLike
from typing import List, Optional


def convert_adni_tau_pet(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: List[str],
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert Tau PET images of ADNI into BIDS format.

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
        Default=1
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        load_clinical_csv,
        paths_to_bids,
    )
    from clinica.utils.stream import cprint

    cprint(
        f"Calculating paths of TAU PET images. Output will be stored in {conversion_dir}."
    )
    images = compute_tau_pet_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of TAU PET images found. Exporting images into BIDS ...")
    paths_to_bids(
        images, destination_dir, "tau", mod_to_update=mod_to_update, n_procs=n_procs
    )
    cprint(msg="TAU PET conversion done.", lvl="debug")


def compute_tau_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute the paths to Tau PET images.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns: pandas Dataframe containing the path for each Tau PET image
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        get_images_pet,
        load_clinical_csv,
    )
    from clinica.utils.pet import Tracer

    pet_tau_col = [
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
    pet_tau_df = pd.DataFrame(columns=pet_tau_col)
    pet_tau_dfs_list = []

    # Loading needed .csv files
    tauqc = load_clinical_csv(csv_dir, "TAUQC")
    tauqc3 = load_clinical_csv(csv_dir, "TAUQC3")
    pet_meta_list = load_clinical_csv(csv_dir, "PET_META_LIST")

    for subj in subjs_list:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for TAU PET images for ADNI 2
        tau_qc2_subj = tauqc[(tauqc.SCANQLTY == 1) & (tauqc.RID == int(subj[-4:]))]

        # QC for TAU PET images for ADNI 3
        tau_qc3_subj = tauqc3[(tauqc3.SCANQLTY == 1) & (tauqc3.RID == int(subj[-4:]))]

        # Concatenating visits in both QC files
        tau_qc_subj = pd.concat(
            [tau_qc2_subj, tau_qc3_subj], axis=0, ignore_index=True, sort=False
        )
        tau_qc_subj.rename(columns={"SCANDATE": "EXAMDATE"}, inplace=True)

        sequences_preprocessing_step = ["AV1451 Co-registered, Averaged"]
        subj_dfs_list = get_images_pet(
            subj,
            tau_qc_subj,
            subject_pet_meta,
            pet_tau_col,
            "TAU-PET",
            sequences_preprocessing_step,
        )
        if subj_dfs_list:
            pet_tau_dfs_list += subj_dfs_list

    if pet_tau_dfs_list:
        pet_tau_df = pd.concat(pet_tau_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [("098_S_4275", "m84")]  # Multiple output images

    # Removing known exceptions from images to convert
    if not pet_tau_df.empty:
        error_ind = pet_tau_df.index[
            pet_tau_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        pet_tau_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(pet_tau_df, source_dir, "TAU", "I", "Image_ID")
    images.to_csv(
        path.join(conversion_dir, f"{Tracer.AV1451.value}_pet_paths.tsv"),
        sep="\t",
        index=False,
    )

    return images
