# coding: utf-8

"""Module for converting Tau PET of ADNI."""


def convert_adni_tau_pet(
    source_dir, csv_dir, dest_dir, conversion_dir, subjs_list=None, mod_to_update=False
):
    """Convert Tau PET images of ADNI into BIDS format.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        conversion_dir: path to the TSV files including the paths to original images
        subjs_list: subjects list
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted again.
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids
    from clinica.utils.stream import cprint

    if not subjs_list:
        adni_merge_path = path.join(csv_dir, "ADNIMERGE.csv")
        adni_merge = pd.read_csv(adni_merge_path, sep=",", low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint(
        f"Calculating paths of TAU PET images. Output will be stored in {conversion_dir}."
    )
    images = compute_tau_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir)
    cprint("Paths of TAU PET images found. Exporting images into BIDS ...")
    paths_to_bids(images, dest_dir, "tau", mod_to_update=mod_to_update)
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
    tauqc = pd.read_csv(path.join(csv_dir, "TAUQC.csv"), sep=",", low_memory=False)
    tauqc3 = pd.read_csv(path.join(csv_dir, "TAUQC3.csv"), sep=",", low_memory=False)
    pet_meta_list = pd.read_csv(
        path.join(csv_dir, "PET_META_LIST.csv"), sep=",", low_memory=False
    )

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
        path.join(conversion_dir, f"{Tracer.AV1451}_pet_paths.tsv"),
        sep="\t",
        index=False,
    )

    return images
