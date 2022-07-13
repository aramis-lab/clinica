# coding: utf-8

"""Module for converting FDG PET of ADNI."""


def convert_adni_fdg_pet(
    source_dir, csv_dir, dest_dir, conversion_dir, subjs_list=None, mod_to_update=False
):
    """Convert FDG PET images of ADNI into BIDS format.

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
        f"Calculating paths of FDG PET images. Output will be stored in {conversion_dir}."
    )
    images = compute_fdg_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir)
    cprint("Paths of FDG PET images found. Exporting images into BIDS ...")
    paths_to_bids(images, dest_dir, "fdg", mod_to_update=mod_to_update)
    cprint(msg="FDG PET conversion done.", lvl="debug")


def compute_fdg_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute the paths to the FDG PET images and store them in a TSV file.

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
    )
    from clinica.utils.pet import Tracer

    pet_fdg_col = [
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
    pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
    pet_fdg_dfs_list = []

    # Loading needed .csv files
    petqc = pd.read_csv(path.join(csv_dir, "PETQC.csv"), sep=",", low_memory=False)
    petqc3 = pd.read_csv(path.join(csv_dir, "PETC3.csv"), sep=",", low_memory=False)
    pet_meta_list = pd.read_csv(
        path.join(csv_dir, "PET_META_LIST.csv"), sep=",", low_memory=False
    )

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for FDG PET images for ADNI 1, GO and 2
        pet_qc_1go2_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subj[-4:]))]

        # QC for FDG PET images for ADNI 3
        pet_qc3_subj = petqc3[(petqc3.SCANQLTY == 1) & (petqc3.RID == int(subj[-4:]))]
        pet_qc3_subj.insert(0, "EXAMDATE", pet_qc3_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        pet_qc_subj = pd.concat(
            [pet_qc_1go2_subj, pet_qc3_subj], axis=0, ignore_index=True, sort=False
        )

        sequences_preprocessing_step = ["Co-registered, Averaged"]
        subj_dfs_list = get_images_pet(
            subj,
            pet_qc_subj,
            subject_pet_meta,
            pet_fdg_col,
            "FDG-PET",
            sequences_preprocessing_step,
        )
        if subj_dfs_list:
            pet_fdg_dfs_list += subj_dfs_list

    if pet_fdg_dfs_list:
        # Concatenating dataframes into one
        pet_fdg_df = pd.concat(pet_fdg_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # NONAME.nii
        ("031_S_0294", "bl"),
        ("037_S_1421", "m36"),
        ("037_S_1078", "m36"),
        # Empty folders
        ("941_S_1195", "m48"),
        ("005_S_0223", "m12"),
    ]

    # Removing known exceptions from images to convert
    if not pet_fdg_df.empty:
        error_ind = pet_fdg_df.index[
            pet_fdg_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        pet_fdg_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_fdg_df, source_dir, "FDG", "I", "Image_ID")
    images.to_csv(
        path.join(conversion_dir, f"{Tracer.FDG}_pet_paths.tsv"), sep="\t", index=False
    )

    return images
