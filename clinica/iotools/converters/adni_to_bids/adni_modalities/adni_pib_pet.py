# coding: utf-8

"""
 Module for converting PIB PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_pib_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """Convert PIB PET images of ADNI into BIDS format

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    """

    import pandas as pd
    from os import path
    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids
    from colorama import Fore

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of PIB PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_pib_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of PIB PET images found. Exporting images into BIDS ...')
    paths_to_bids(images, dest_dir, 'pib')
    cprint(Fore.GREEN + 'PIB PET conversion done.' + Fore.RESET)


def compute_pib_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """Compute the paths to the PIB PET images and store them in a tsv file

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns:
        images: a dataframe with all the paths to the PET images that will be converted into BIDS

    """

    import pandas as pd
    import os
    from os import path
    from clinica.iotools.converters.adni_to_bids.adni_utils import get_images_pet, find_image_path

    pet_pib_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']
    pet_pib_df = pd.DataFrame(columns=pet_pib_col)
    pet_pib_dfs_list = []

    # Loading needed .csv files
    pibqc = pd.read_csv(path.join(csv_dir, 'PIBQC.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for PIB PET images
        pet_qc_subj = pibqc[(pibqc.PASS == 1) & (pibqc.RID == int(subj[-4:]))]

        sequences_preprocessing_step = ['PIB Co-registered, Averaged']
        subj_dfs_list = get_images_pet(subj, pet_qc_subj, subject_pet_meta, pet_pib_col, 'PIB-PET',
                                       sequences_preprocessing_step, viscode_field="VISCODE")
        if subj_dfs_list:
            pet_pib_dfs_list += subj_dfs_list

    if pet_pib_dfs_list:
        pet_pib_df = pd.concat(pet_pib_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = []

    # Removing known exceptions from images to convert
    if not pet_pib_df.empty:
        error_ind = pet_pib_df.index[pet_pib_df.apply(lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors),
                                                      axis=1)]
        pet_pib_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_pib_df, source_dir, 'PIB', 'I', 'Image_ID')

    pib_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(pib_csv_path):
        os.mkdir(pib_csv_path)
    images.to_csv(path.join(pib_csv_path, 'pib_pet_paths.tsv'), sep='\t', index=False)

    return images
