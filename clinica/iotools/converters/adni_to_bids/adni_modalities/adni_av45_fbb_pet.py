# coding: utf-8

"""
 Module for converting AV45 and Florbetaben PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_av45_fbb_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert AV-45 and Florbetaben PET images of ADNI into BIDS format

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

    cprint('Calculating paths of AV45 and Florbetaben PET images. Output will be stored in %s.' %
           path.join(dest_dir, 'conversion_info'))
    images = compute_av45_fbb_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of AV45 and Florbetaben PET images found. Exporting images into BIDS ...')
    paths_to_bids(images, dest_dir, 'av45_fbb')
    cprint(Fore.GREEN + 'AV45 and Florbetaben PET conversion done.' + Fore.RESET)


def compute_av45_fbb_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to the AV45 and Florbetaben PET images and store them in a tsv file

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

    pet_amyloid_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID', 'Series_ID',
                       'Image_ID', 'Original', 'Tracer']
    pet_amyloid_df = pd.DataFrame(columns=pet_amyloid_col)
    pet_amyloid_dfs_list = []

    # Loading needed .csv files
    av45qc = pd.read_csv(path.join(csv_dir, 'AV45QC.csv'), sep=',', low_memory=False)
    amyqc = pd.read_csv(path.join(csv_dir, 'AMYQC.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for AV45 PET images for ADNI 1, GO and 2
        av45_qc_subj = av45qc[(av45qc.PASS == 1) & (av45qc.RID == int(subj[-4:]))]

        # QC for Amyloid PET images for ADNI 3
        amy_qc_subj = amyqc[(amyqc.SCANQLTY == 1) & (amyqc.RID == int(subj[-4:]))]
        amy_qc_subj.insert(0, 'EXAMDATE', amy_qc_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        amyloid_qc_subj = pd.concat([av45_qc_subj, amy_qc_subj], axis=0, ignore_index=True, sort=False)

        sequences_preprocessing_step = ['AV45 Co-registered, Averaged',
                                        'FBB Co-registered, Averaged']
        subj_dfs_list = get_images_pet(subj, amyloid_qc_subj, subject_pet_meta, pet_amyloid_col, 'Amyloid-PET',
                                       sequences_preprocessing_step)
        if subj_dfs_list:
            pet_amyloid_dfs_list += subj_dfs_list

    if pet_amyloid_dfs_list:
        pet_amyloid_df = pd.concat(pet_amyloid_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # Eq_1
                         ('128_S_2220', 'm48'),
                         # Several output images
                         ('098_S_4275', 'm84')]

    # Removing known exceptions from images to convert
    if not pet_amyloid_df.empty:
        error_ind = pet_amyloid_df.index[pet_amyloid_df.apply(
            lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1)]
        pet_amyloid_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_amyloid_df, source_dir, 'Amyloid', 'I', 'Image_ID')

    amyloid_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(amyloid_csv_path):
        os.mkdir(amyloid_csv_path)
    images.to_csv(path.join(amyloid_csv_path, 'amyloid_pet_paths.tsv'), sep='\t', index=False)

    return images
