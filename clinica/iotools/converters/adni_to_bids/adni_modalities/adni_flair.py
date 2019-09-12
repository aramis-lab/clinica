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

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of FLAIR images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_flair_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of FLAIR images found. Exporting images into BIDS ...')
    flair_paths_to_bids(images, dest_dir)
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

    import operator
    from os import path, mkdir
    from functools import reduce
    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import find_image_path, visits_to_timepoints_mrilist

    flair_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                    'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner']
    flair_df = pd.DataFrame(columns=flair_col_df)

    # Loading needed .csv files
    adni_merge = pd.read_csv(path.join(csv_dir, 'ADNIMERGE.csv'), sep=',', low_memory=False)

    mayo_mri_qc = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv'), sep=',', low_memory=False)
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'AFL']

    mri_list = pd.read_csv(path.join(csv_dir, 'MRILIST.csv'), sep=',', low_memory=False)

    # Selecting FLAIR DTI images that are not TODO images
    mri_list = mri_list[mri_list.SEQUENCE.map(lambda x: x.lower().find('flair') > -1)]
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
        visits = visits_to_timepoints_mrilist(subj, mri_list_subj, adnimerge_subj, 'FLAIR')

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            flair = flair_image(subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj)

            if flair is not None:
                row_to_append = pd.DataFrame(flair, index=['i', ])
                # TODO Replace iteratively appending by pandas.concat
                flair_df = flair_df.append(row_to_append, ignore_index=True, sort=False)

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
                         ('053_S_5272', 'm24'),

    ]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((flair_df.Subject_ID == conv_error[0])
                             & (flair_df.VISCODE == conv_error[1]))
    if error_indices:
        indices_to_remove = flair_df.index[reduce(operator.or_, error_indices, False)]
        flair_df.drop(indices_to_remove, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(flair_df, source_dir, 'FLAIR', 'S', 'Series_ID')

    flair_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(flair_tsv_path):
        mkdir(flair_tsv_path)
    images.to_csv(path.join(flair_tsv_path, 'flair_paths.tsv'), sep='\t', index=False)

    return images


def flair_paths_to_bids(images, dest_dir, mod_to_update=False):
    """
    Convert FLAIR images

    Args:
        images: dataframe returned by the method compute_flair_paths
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

    subjs_list = [sub for sub in images['Subject_ID'].unique()
                  if sub == sub and 'S' in sub]

    counter = Value('i', 0)
    partial_generate_subject_files = partial(generate_subject_files,
                                             images=images,
                                             dest_dir=dest_dir,
                                             mod_to_update=mod_to_update)
    poolrunner = Pool(cpu_count(), initializer=init, initargs=(counter,))
    poolrunner.map(partial_generate_subject_files, subjs_list)
    del counter


def flair_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
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
    import clinica.iotools.bids_utils as bids
    import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
    from clinica.utils.stream import cprint
    import subprocess
    import os
    import shutil
    from os import path
    from glob import glob

    alpha_id = bids.remove_space_and_symbols(subj)
    bids_id = 'sub-ADNI' + alpha_id
    # Extract the list of sessions available from the flair paths files, removing the duplicates
    sess_list = images[(images['Subject_ID'] == subj)]['VISCODE'].unique()

    if not os.path.exists(path.join(dest_dir, bids_id)):
        os.mkdir(path.join(dest_dir, bids_id))

    # For each session available, create the folder if doesn't exist and convert the files
    for ses in sess_list:
        with counter.get_lock():
            counter.value += 1
        cprint('[FLAIR] Processing subject ' + str(subj)
               + ' - session ' + ses + ', ' + str(counter.value)
               + ' / ' + str(len(images)))
        ses_bids = adni_utils.viscode_to_session(ses)
        bids_ses_id = 'ses-' + ses_bids
        bids_file_name = bids_id + '_ses-' + ses_bids
        ses_path = path.join(dest_dir, bids_id, bids_ses_id)

        if mod_to_update:
            if os.path.exists(path.join(ses_path, 'FLAIR')):
                shutil.rmtree(path.join(ses_path, 'FLAIR'))

        if not os.path.exists(ses_path):
            os.mkdir(ses_path)

        flair_info = images[
            (images['Subject_ID'] == subj) & (images['VISCODE'] == ses)]

        # For the same subject, same session there could be multiple flair with different acq label
        for j in range(len(flair_info)):
            flair_subj = flair_info.iloc[j]
            if type(flair_subj['Path']) != float and flair_subj['Path'] != '':
                if not os.path.exists(path.join(ses_path, 'FLAIR')):
                    os.mkdir(path.join(ses_path, 'FLAIR'))
                flair_path = flair_subj['Path']
                bids_name = bids_file_name + '_FLAIR'
                bids_dest_dir = path.join(ses_path, 'FLAIR')
                image_id = flair_subj.Image_ID

                # If the original image is a DICOM, check if contains two DICOM
                # inside the same folder
                if flair_subj.Is_Dicom:
                    flair_path = adni_utils.check_two_dcm_folder(flair_path,
                                                                 bids_dest_dir,
                                                                 image_id)

                if not os.path.exists(bids_dest_dir):
                    os.mkdir(dest_dir)
                command = 'dcm2niix -b y -z y -o ' + bids_dest_dir + ' -f ' + bids_name + ' ' + flair_path
                subprocess.run(command,
                               shell=True,
                               stderr=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL)

                # If dcm2niix didn't work use dcm2nii
                if not os.path.exists(path.join(bids_dest_dir, bids_name + '.nii.gz')):
                    cprint('\tConversion with dcm2niix failed, trying with dcm2nii')

                    # Find all the files eventually created by dcm2niix and remove them
                    flair_dcm2niix = glob(
                        path.join(bids_dest_dir, bids_name + '*'))
                    for d in flair_dcm2niix:
                        os.remove(d)

                    command = 'dcm2nii -a n -d n -e n -i y -g y -p n -m n -r n -x n -o ' + bids_dest_dir + ' ' + flair_path
                    subprocess.run(command,
                                   shell=True,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                    nii_file = path.join(bids_dest_dir,
                                         subj.replace('_', '') + '.nii.gz')
                    if os.path.exists(nii_file):
                        os.rename(nii_file, path.join(bids_dest_dir,
                                                      bids_name + '.nii.gz'))
                    else:
                        cprint('WARNING: CONVERSION FAILED...')


def check_exceptions(bids_dir):
    from os import path
    import pandas as pd
    from glob import glob

    flair_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'flair_paths.tsv'), sep='\t')

    flair_paths = flair_paths[flair_paths.Path.notnull()]

    flair_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in flair_paths.Subject_ID.to_list()]
    flair_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in flair_paths.VISCODE.to_list()]

    count = 0
    count_wrong = 0
    count_eq = 0

    for r in flair_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'FLAIR')
        image_pattern = path.join(image_dir, '%s_%s_*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)

        if sum(['Eq' in f for f in files_list]):
            # print("Eq images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            # print("('%s', '%s')," % (image.BIDS_SubjID[-8:].replace('S', '_S_'),
            #                          image.BIDS_Session[4:].replace('M00', 'bl').lower()))
            count_eq += 1
            continue

        if not files_list:
            # print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 2:
            print("Wrong files count for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)
            count_wrong += 1

    print(count)
    print(count_wrong)
    print(count_eq)
