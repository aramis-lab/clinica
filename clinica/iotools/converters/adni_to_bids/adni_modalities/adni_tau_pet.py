# coding: utf-8

"""
 Module for converting Tau PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_tau_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert Tau PET images of ADNI into BIDS format

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

    cprint('Calculating paths of TAU PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_tau_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of TAU PET images found. Exporting images into BIDS ...')
    paths_to_bids(images, dest_dir, 'tau')
    cprint(Fore.GREEN + 'TAU PET conversion done.' + Fore.RESET)


def compute_tau_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to Tau PET images

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns: pandas Dataframe containing the path for each Tau PET image

    """

    import pandas as pd
    import os
    from os import path
    from clinica.iotools.converters.adni_to_bids.adni_utils import get_images_pet, find_image_path

    pet_tau_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']
    pet_tau_df = pd.DataFrame(columns=pet_tau_col)
    pet_tau_dfs_list = []

    # Loading needed .csv files
    tauqc = pd.read_csv(path.join(csv_dir, 'TAUQC.csv'), sep=',', low_memory=False)
    tauqc3 = pd.read_csv(path.join(csv_dir, 'TAUQC3.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.empty:
            continue

        # QC for TAU PET images for ADNI 2
        tau_qc2_subj = tauqc[(tauqc.SCANQLTY == 1) & (tauqc.RID == int(subj[-4:]))]

        # QC for TAU PET images for ADNI 3
        tau_qc3_subj = tauqc3[(tauqc3.SCANQLTY == 1) & (tauqc3.RID == int(subj[-4:]))]

        # Concatenating visits in both QC files
        tau_qc_subj = pd.concat([tau_qc2_subj, tau_qc3_subj], axis=0, ignore_index=True, sort=False)
        tau_qc_subj.rename(columns={"SCANDATE": "EXAMDATE"}, inplace=True)

        sequences_preprocessing_step = ['AV1451 Co-registered, Averaged']
        subj_dfs_list = get_images_pet(subj, tau_qc_subj, subject_pet_meta, pet_tau_col, 'TAU-PET',
                                       sequences_preprocessing_step)
        if subj_dfs_list:
            pet_tau_dfs_list += subj_dfs_list

    if pet_tau_dfs_list:
        pet_tau_df = pd.concat(pet_tau_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # Multiple output images
                         ('098_S_4275', 'm84')]

    # Removing known exceptions from images to convert
    if not pet_tau_df.empty:
        error_ind = pet_tau_df.index[pet_tau_df.apply(lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors),
                                                      axis=1)]
        pet_tau_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(pet_tau_df, source_dir, 'TAU', 'I', 'Image_ID')

    tau_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(tau_csv_path):
        os.mkdir(tau_csv_path)
    images.to_csv(path.join(tau_csv_path, 'tau_pet_paths.tsv'), sep='\t', index=False)

    return images
