# coding: utf-8

"""
 Module for converting FLAIR of ADNI
"""

__author__ = "Simona Bottani and Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_flair(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert FLAIR images of ADNI into BIDS format

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    """

    import pandas as pd
    from os import path
    from clinica.utils.stream import cprint
    from colorama import Fore
    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of FLAIR images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_flair_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of FLAIR images found. Exporting images into BIDS ...')
    # flair_paths_to_bids(images, dest_dir)
    paths_to_bids(images, dest_dir, 'flair')
    cprint(Fore.GREEN + 'FLAIR conversion done.' + Fore.RESET)


def compute_flair_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to the FLAIR images and store them in a tsv file

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns:
        images: a dataframe with all the paths to the FLAIR images that will be converted into BIDS

    """

    from os import path, mkdir
    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import find_image_path, visits_to_timepoints

    flair_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                    'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner']
    flair_df = pd.DataFrame(columns=flair_col_df)
    flair_dfs_list = []

    # Loading needed .csv files
    adni_merge = pd.read_csv(path.join(csv_dir, 'ADNIMERGE.csv'), sep=',', low_memory=False)

    mayo_mri_qc = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv'), sep=',', low_memory=False)
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'AFL']

    mri_list = pd.read_csv(path.join(csv_dir, 'MRILIST.csv'), sep=',', low_memory=False)

    # Selecting FLAIR DTI images that are not MPR
    mri_list = mri_list[mri_list.SEQUENCE.str.contains('flair', case=False, na=False)]
    unwanted_sequences = ['_MPR_']
    mri_list = mri_list[mri_list.SEQUENCE.map(lambda x: not any(subs in x for subs in unwanted_sequences))]

    for subj in subjs_list:

        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values('SCANDATE')

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subj, mri_list_subj, adnimerge_subj, 'FLAIR')

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            flair = flair_image(subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj)

            if flair is not None:
                row_to_append = pd.DataFrame(flair, index=['i', ])
                flair_dfs_list.append(row_to_append)

    if flair_dfs_list:
        flair_df = pd.concat(flair_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # Eq_1 images
                         ('141_S_0767', 'm84'),
                         ('067_S_5205', 'bl'),
                         ('127_S_4928', 'm24'),
                         ('024_S_4674', 'm06'),
                         ('123_S_2363', 'm24'),
                         ('053_S_4578', 'm48'),
                         ('128_S_4586', 'm48'),
                         ('053_S_4813', 'm48'),
                         ('053_S_5272', 'm24')]

    # Removing known exceptions from images to convert
    if not flair_df.empty:
        error_ind = flair_df.index[flair_df.apply(lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1)]
        flair_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(flair_df, source_dir, 'FLAIR', 'S', 'Series_ID')

    flair_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(flair_tsv_path):
        mkdir(flair_tsv_path)
    images.to_csv(path.join(flair_tsv_path, 'flair_paths.tsv'), sep='\t', index=False)

    return images


def flair_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
    """
    One image among those in the input list is chosen according to QC
    and then correspoding metadata is extracted to a dictionary

    Args:
        subject_id: Subject identifier
        timepoint: Visit code
        visit_str: Visit name
        visit_mri_list: List of images metadata
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: dictionary - contains image metadata

    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars, select_image_qc

    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if sel_image is None:
        return None

    sel_scan = visit_mri_list[visit_mri_list.IMAGEUID == sel_image].iloc[0]

    image_dict = {'Subject_ID': subject_id,
                  'VISCODE': timepoint,
                  'Visit': visit_str,
                  'Sequence': replace_sequence_chars(sel_scan.SEQUENCE),
                  'Scan_Date': sel_scan['SCANDATE'],
                  'Study_ID': str(int(sel_scan.STUDYID)),
                  'Series_ID': str(int(sel_scan.SERIESID)),
                  'Image_ID': str(int(sel_scan.IMAGEUID)),
                  'Field_Strength': sel_scan.MAGSTRENGTH}

    return image_dict
