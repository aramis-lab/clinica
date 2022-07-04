# coding: utf-8

"""Module for converting PIB PET of ADNI."""


def convert_adni_pib_pet(
    source_dir, csv_dir, dest_dir, conversion_dir, subjs_list=None, mod_to_update=False
):
    """Convert PIB PET images of ADNI into BIDS format.

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
        f"Calculating paths of PIB PET images. Output will be stored in {conversion_dir}."
    )
    images = compute_pib_pet_paths(source_dir, csv_dir, subjs_list, conversion_dir)
    cprint("Paths of PIB PET images found. Exporting images into BIDS ...")
    paths_to_bids(images, dest_dir, "pib", mod_to_update=mod_to_update)
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
    pibqc = pd.read_csv(path.join(csv_dir, "PIBQC.csv"), sep=",", low_memory=False)
    pet_meta_list = pd.read_csv(
        path.join(csv_dir, "PET_META_LIST.csv"), sep=",", low_memory=False
    )

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
