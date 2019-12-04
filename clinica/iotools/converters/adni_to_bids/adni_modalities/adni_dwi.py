# coding: utf-8

"""
 Module for converting DWI of ADNI
"""

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_dwi(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert DW images of ADNI into BIDS format

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

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of DWI images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_dwi_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of DWI images found. Exporting images into BIDS ...')
    dwi_paths_to_bids(images, dest_dir)
    cprint(Fore.GREEN + 'DWI conversion done.' + Fore.RESET)


def compute_dwi_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute paths to DW images to convert to BIDS

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns:
        images: pandas dataframe that contains the path to all the DW images to convert

    """

    from os import path, mkdir
    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import find_image_path, visits_to_timepoints_mrilist

    dwi_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                  'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength']
    dwi_df = pd.DataFrame(columns=dwi_col_df)
    dwi_dfs_list = []

    # Loading needed .csv files
    adni_merge = pd.read_csv(path.join(csv_dir, 'ADNIMERGE.csv'), sep=',', low_memory=False)

    mayo_mri_qc = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv'), sep=',', low_memory=False)
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'DTI']

    mri_list = pd.read_csv(path.join(csv_dir, 'MRILIST.csv'), sep=',', low_memory=False)

    # Selecting only DTI images that are not Multiband, processed or enchanced images
    mri_list = mri_list[mri_list.SEQUENCE.str.contains('dti', case=False, na=False)]
    unwanted_sequences = ['MB', 'ADC', 'FA', 'TRACEW', 'Enhanced', 'Reg']
    mri_list = mri_list[mri_list.SEQUENCE.map(lambda x: not any(subs in x for subs in unwanted_sequences))]

    for subj in subjs_list:

        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values('SCANDATE')

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints_mrilist(subj, mri_list_subj, adnimerge_subj, "DWI")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            axial = dwi_image(subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj)

            if axial is not None:
                row_to_append = pd.DataFrame(axial, index=['i', ])
                dwi_dfs_list.append(row_to_append)

    if dwi_dfs_list:
        dwi_df = pd.concat(dwi_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [('029_S_2395', 'm60'),
                         ('029_S_0824', 'm108'),
                         ('029_S_0914', 'm108'),
                         ('027_S_2219', 'm36'),
                         ('129_S_2332', 'm12'),
                         ('029_S_4384', 'm48'),
                         ('029_S_4385', 'm48'),
                         ('029_S_4585', 'm48'),
                         ('016_S_4591', 'm24'),
                         ('094_S_4630', 'm06'),
                         ('094_S_4649', 'm06'),
                         ('029_S_5219', 'm24'),
                         ('094_S_2238', 'm48'),
                         ('129_S_4287', 'bl'),
                         ('007_S_4611', 'm03'),
                         ('016_S_4638', 'bl'),
                         ('027_S_5118', 'bl'),
                         ('098_S_4018', 'bl'),
                         ('098_S_4003', 'm12'),
                         ('016_S_4584', 'm24'),
                         ('016_S_5007', 'm12'),
                         ('129_S_2347', 'm06'),
                         ('129_S_4220', 'bl'),
                         ('007_S_2058', 'm12'),
                         ('016_S_2007', 'm06'),
                         ('020_S_6358', 'bl'),
                         ('114_S_6039', 'm12'),
                         ('114_S_6057', 'bl'),
                         ('153_S_6274', 'bl'),
                         ('006_S_4485', 'm84'),
                         ('153_S_6237', 'bl'),
                         ('153_S_6336', 'bl'),
                         ('153_S_6450', 'bl'),
                         ('003_S_4441', 'm12')]

    # Removing known exceptions from images to convert
    if dwi_df.shape[0] > 0:
        error_ind = dwi_df.index[dwi_df.apply(lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1)]
        dwi_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(dwi_df, source_dir, 'DWI', 'S', 'Series_ID')

    dwi_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(dwi_tsv_path):
        mkdir(dwi_tsv_path)
    images.to_csv(path.join(dwi_tsv_path, 'dwi_paths.tsv'), sep='\t', index=False)

    return images


def dwi_paths_to_bids(images, dest_dir, mod_to_update=False):
    """
    Convert DWI images

    Args:
        images: dataframe returned by the method compute_dwi_paths
        dest_dir: path to the destination directory
        mod_to_update: if is true and an image is already existing it will overwrite the old version

    """
    from multiprocessing import cpu_count, Pool, Value
    from functools import partial

    counter = None

    def init(args):
        ''' store the counter for later use '''
        global counter
        counter = args

    subjs_list = images['Subject_ID'].unique()
    counter = Value('i', 0)
    total = len(images)
    partial_generate_subject_files = partial(generate_subject_files,
                                             images=images,
                                             dest_dir=dest_dir,
                                             mod_to_update=mod_to_update)
    poolrunner = Pool(cpu_count(), initializer=init, initargs=(counter,))
    poolrunner.map(partial_generate_subject_files, subjs_list)
    del counter


def dwi_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
    """

    Args:
        subject_id:
        timepoint:
        visit_str:
        visit_mri_list:
        mri_qc_subj:

    Returns:

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


def generate_subject_files(subj, images, dest_dir, mod_to_update):
    """

    Args:
        subj:
        images:
        dest_dir:
        mod_to_update:

    Returns:

    """
    import clinica.iotools.bids_utils as bids
    import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
    from clinica.utils.stream import cprint
    import subprocess
    import os
    import shutil
    from os import path
    from glob import glob
    import pandas as pd
    from colorama import Fore

    global counter
    alpha_id = bids.remove_space_and_symbols(subj)
    bids_id = 'sub-ADNI' + alpha_id
    # Extract the list of sessions available from the dwi paths files, removing the duplicates
    sess_list = images[(images['Subject_ID'] == subj)]['VISCODE'].unique()

    if not os.path.exists(path.join(dest_dir, bids_id)):
        os.mkdir(path.join(dest_dir, bids_id))

    # For each session available, create the folder if doesn't exist and convert the files
    for ses in sess_list:
        with counter.get_lock():
            counter.value += 1
        cprint('[DWI] Processing subject ' + str(subj)
               + ' - session ' + ses + ', ' + str(counter.value)
               + ' / ' + str(len(images)))
        ses_bids = adni_utils.viscode_to_session(ses)
        bids_ses_id = 'ses-' + ses_bids
        bids_file_name = bids_id + '_ses-' + ses_bids
        ses_path = path.join(dest_dir, bids_id, bids_ses_id)

        if mod_to_update:
            if os.path.exists(path.join(ses_path, 'dwi')):
                shutil.rmtree(path.join(ses_path, 'dwi'))

        dwi_info = images[(images['Subject_ID'] == subj) & (images['VISCODE'] == ses)]

        # For the same subject, same session there could be multiple dwi with different acq label
        for j in range(len(dwi_info)):

            dwi_subj = dwi_info.iloc[j]

            if not pd.isnull(dwi_subj['Path']) and dwi_subj['Path']:

                if not os.path.exists(path.join(ses_path, 'dwi')):
                    os.makedirs(path.join(ses_path, 'dwi'))
                dwi_path = dwi_subj['Path']

                bids_name = bids_file_name + '_acq-axial_dwi'

                bids_dest_dir = path.join(ses_path, 'dwi')

                if not os.path.exists(bids_dest_dir):
                    os.makedirs(bids_dest_dir)
                command = 'dcm2niix -b y -z y -o ' + bids_dest_dir + ' -f ' + bids_name + ' ' + dwi_path
                subprocess.run(command,
                               shell=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)

                # Removing ADC images
                adc_image = path.join(bids_dest_dir, bids_name + '_ADC.nii.gz')
                if os.path.exists(adc_image):
                    os.remove(adc_image)

                # If dcm2niix didn't work use dcm2nii
                # We do not chek if .json exists given that dcm2nii does not generate .json file anyways
                if (not os.path.exists(path.join(bids_dest_dir, bids_name + '.nii.gz')))\
                        or (not os.path.exists(path.join(bids_dest_dir, bids_name + '.bvec')))\
                        or (not os.path.exists(path.join(bids_dest_dir, bids_name + '.bval'))):

                    cprint('WARNING: Conversion with dcm2niix failed, trying with dcm2nii '
                           'for subject ' + subj + ' and session ' + ses)

                    # Find all the files eventually created by dcm2niix and remove them
                    dwi_dcm2niix = glob(path.join(bids_dest_dir, bids_name + '*'))

                    for d in dwi_dcm2niix:
                        os.remove(d)

                    command = 'dcm2nii -a n -d n -e n -i y -g y -p n -m n -r n -x n -o %s %s' \
                              % (bids_dest_dir, dwi_path)
                    subprocess.run(command,
                                   shell=True,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)

                    nii_file = path.join(bids_dest_dir, subj.replace('_', '') + '.nii.gz')
                    bvec_file = path.join(bids_dest_dir, subj.replace('_', '') + '.bvec')
                    bval_file = path.join(bids_dest_dir, subj.replace('_', '') + '.bval')

                    if os.path.exists(bvec_file) and os.path.exists(bval_file):
                        os.rename(bvec_file, path.join(bids_dest_dir,
                                                       bids_name + '.bvec'))
                        os.rename(bval_file, path.join(bids_dest_dir,
                                                       bids_name + '.bval'))
                    else:
                        cprint('WARNING: bvec and bval not generated by dcm2nii '
                               'for subject ' + subj + ' and session ' + ses)

                    if os.path.exists(nii_file):
                        os.rename(nii_file, path.join(bids_dest_dir,
                                                      bids_name + '.nii.gz'))
                    else:
                        cprint(Fore.RED + 'WARNING: Conversion of the dicom failed '
                               'for subject ' + subj + ' and session ' + ses + '. '
                               'Image path: ' + dwi_path + Fore.RESET)

                        # Removing folder and generated files
                        if os.path.exists(bids_dest_dir):
                            shutil.rmtree(bids_dest_dir, ignore_errors=True)

                # Case when JSON file was not generated by dcm2niix
                elif not os.path.exists(path.join(bids_dest_dir, bids_name + '.json')):
                    cprint('WARNING: JSON file not generated by dcm2niix '
                           'for subject ' + subj + ' and session ' + ses)
