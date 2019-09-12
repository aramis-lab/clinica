# coding: utf-8

"""
 Module for converting fMRI of ADNI
"""

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_fmri(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert fMR images of ADNI into BIDS format

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

    cprint('Calculating paths of fMRI images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_fmri_path(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of fMRI images found. Exporting images into BIDS ...')
    fmri_paths_to_bids(dest_dir, images)
    cprint(Fore.GREEN + 'fMRI conversion done.' + Fore.RESET)


def compute_fmri_path(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to fMR images

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns: pandas Dataframe containing the path for each fmri

    """

    import operator
    from functools import reduce
    from os import path, mkdir
    import pandas as pd
    from clinica.iotools.converters.adni_to_bids.adni_utils import find_image_path, visits_to_timepoints_mrilist

    fmri_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                'Study_ID', 'Field_Strength', 'Series_ID', 'Image_ID']
    fmri_df = pd.DataFrame(columns=fmri_col)

    # Loading needed .csv files
    adni_merge = pd.read_csv(path.join(csv_dir, 'ADNIMERGE.csv'), sep=',', low_memory=False)

    mayo_mri_qc = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv'), sep=',', low_memory=False)
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'fMRI']
    mayo_mri_qc.columns = [x.upper() for x in mayo_mri_qc.columns]

    mayo_mri_qc3 = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_QUALITY_ADNI3.csv'), sep=',', low_memory=False)
    mayo_mri_qc3 = mayo_mri_qc3[mayo_mri_qc3.SERIES_TYPE == 'EPB']

    # Concatenating visits in both QC files
    mayo_mri_qc = pd.concat([mayo_mri_qc, mayo_mri_qc3], axis=0, ignore_index=True, sort=False)

    mri_list = pd.read_csv(path.join(csv_dir, 'MRILIST.csv'), sep=',', low_memory=False)

    # Selecting only fMRI images that are not Multiband
    mri_list = mri_list[mri_list.SEQUENCE.str.contains('MRI')]  # 'MRI' includes all fMRI and fMRI scans, but not others
    unwanted_sequences = ['MB']
    mri_list = mri_list[mri_list.SEQUENCE.map(lambda x: not any(subs in x for subs in unwanted_sequences))]

    # We will convert the images for each subject in the subject list
    for subj in subjs_list:

        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values('SCANDATE')

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints_mrilist(subj, mri_list_subj, adnimerge_subj, "fMRI")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            image = fmri_image(subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj)

            if image is not None:
                row_to_append = pd.DataFrame(image, index=['i', ])
                # TODO Replace iteratively appending by pandas.concat
                fmri_df = fmri_df.append(row_to_append, ignore_index=True, sort=False)

    # Exceptions
    # ==========
    conversion_errors = [('006_S_4485', 'm84')]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((fmri_df.Subject_ID == conv_error[0])
                             & (fmri_df.VISCODE == conv_error[1]))

    if error_indices:
        indices_to_remove = fmri_df.index[reduce(operator.or_, error_indices, False)]
        fmri_df.drop(indices_to_remove, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(fmri_df, source_dir, 'fMRI', 'S', 'Series_ID')

    fmri_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(fmri_tsv_path):
        mkdir(fmri_tsv_path)
    images.to_csv(path.join(fmri_tsv_path, 'fmri_paths.tsv'), sep='\t', index=False)

    return images


def fmri_paths_to_bids(dest_dir, fmri_paths, mod_to_update=False):
    """
    Convert the fmri extracted from the fmri_paths to BIDS

    Args:
        dest_dir: path to the input directory
        fmri_paths: path to the BIDS directory
        mod_to_update:  if True overwrite (or create if is missing) all the existing fmri

    """
    from multiprocessing import cpu_count, Pool, Value
    from functools import partial

    counter = None

    def init(args):
        # store the counter for later use
        global counter
        counter = args

    subjs_list = fmri_paths['Subject_ID'].drop_duplicates().values

    counter = Value('i', 0)
    partial_generate_subject_files = partial(generate_subject_files,
                                             fmri_paths=fmri_paths,
                                             dest_dir=dest_dir,
                                             mod_to_update=mod_to_update)
    poolrunner = Pool(cpu_count(), initializer=init, initargs=(counter,))
    poolrunner.map(partial_generate_subject_files, subjs_list)
    del counter


def generate_subject_files(subj, fmri_paths, dest_dir, mod_to_update):
    import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
    import clinica.iotools.bids_utils as bids_utils
    from clinica.utils.stream import cprint
    import os
    from os import path
    from glob import glob

    sess_list = fmri_paths[(fmri_paths['Subject_ID'] == subj)]['VISCODE'].values
    alpha_id = adni_utils.remove_space_and_symbols(subj)
    bids_id = 'sub-ADNI' + alpha_id

    # For each session available, create the folder if doesn't exist and convert the files
    for ses in sess_list:
        with counter.get_lock():
            counter.value += 1
        cprint('[fMRI] Processing subject ' + str(subj)
               + ' - session ' + ses + ', ' + str(counter.value)
               + ' / ' + str(len(fmri_paths)))
        ses_bids = adni_utils.viscode_to_session(ses)
        bids_ses_id = 'ses-' + ses_bids
        bids_file_name = bids_id + '_ses-' + ses_bids
        ses_path = path.join(dest_dir, bids_id, bids_ses_id)

        # If the fmri already exist
        existing_fmri = glob(path.join(ses_path, 'func', '*_bold*'))
        # if mod_to_add:
        #     if len(existing_fmri) > 0:
        #         print 'Fmri already existing. Skipped.'
        #         continue

        if mod_to_update and len(existing_fmri) > 0:
            # print 'Removing the old fmri folder...'
            os.remove(existing_fmri[0])

        if not os.path.exists(ses_path):
            if not os.path.exists(path.join(dest_dir, bids_id)):
                os.mkdir(path.join(dest_dir, bids_id))
            os.mkdir(path.join(dest_dir, bids_id, bids_ses_id))

        fmri_info = fmri_paths[(fmri_paths['Subject_ID'] == subj) & (fmri_paths['VISCODE'] == ses)]
        if not fmri_info['Path'].values[0] == '':
            if type(fmri_info['Path'].values[0]) != float:
                if not os.path.exists(path.join(ses_path, 'func')):
                    os.mkdir(path.join(ses_path, 'func'))
                fmri_path = fmri_info['Path'].values[0]
                dcm_to_convert = adni_utils.check_two_dcm_folder(fmri_path,
                                                                 dest_dir,
                                                                 fmri_info['Image_ID'].values[0])
                if not os.path.isfile(os.path.join(ses_path, 'func', bids_file_name + '_task-rest_bold.nii.gz')):
                    bids_utils.convert_fmri(dcm_to_convert, path.join(ses_path, 'func'), bids_file_name)
                else:
                    cprint("Images already converted")

                # Delete the temporary folder used for copying fmri with 2 subjects inside the DICOM folder
                adni_utils.remove_tmp_dmc_folder(dest_dir, fmri_info['Image_ID'].values[0])


def fmri_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
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

    mri_qc_subj.columns = [x.lower() for x in mri_qc_subj.columns]
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


def check_exceptions(bids_dir):

    from os import path
    import pandas as pd
    from glob import glob

    fmri_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t')

    fmri_paths = fmri_paths[fmri_paths.Path.notnull()]

    fmri_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in fmri_paths.Subject_ID.to_list()]
    fmri_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in fmri_paths.VISCODE.to_list()]

    count = 0
    count_wrong = 0
    name_wrong = 0
    no_dir = 0

    for r in fmri_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'func')

        if not path.isdir(image_dir):
            no_dir += 1
            # continue

        image_pattern = path.join(image_dir, '%s_%s_*bold*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)

        if not files_list:
            # print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) != 2:
            print("Wrong files count for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)
            count_wrong += 1
        elif sum([not f.endswith(('_task-rest_bold.json', '_task-rest_bold.nii.gz')) for f in files_list]) > 0:
            name_wrong += 1
            print(files_list)

    print(count)
    print(count_wrong)
    print(name_wrong)
    print(no_dir)
