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

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination directory
        subjs_list: subjects list

    Returns:

    """
    import pandas as pd
    from os import path
    from clinica.utils.stream import cprint
    from colorama import Fore

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of fMRI images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    if path.isfile(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv')):
        images = pd.io.parsers.read_csv(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t')
    else:

        images = compute_fmri_path(source_dir, csv_dir, dest_dir, subjs_list)

    cprint('Paths of fMRI images found. Exporting images into BIDS ...')

    fmri_paths_to_bids(dest_dir, images)
    cprint(Fore.GREEN + 'fMRI conversion done.' + Fore.RESET)


def compute_fmri_path(source_dir, clinical_dir, dest_dir, subjs_list):
    """
    Compute the paths to fmri images.

    The fmri images to convert into BIDS are chosen in the following way:
        - Extract the list of subjects from MAYOADIRL_MRI_FMRI_09_15_16.csv
        - Select the only the scans that came from PHILIPS machine (field Scanner from IDA_MR_Metadata_Listing.csv)
        - Discard all the subjects with column  series_quality = 4  (4 means that the scan is not usable) in MAYOADIRL_MRI_IMAGEQC_12_08_15.csv

    In case of multiple scans for the same session, same date the one to convert is chosen with the following criteria:
        - Check if in the file MAYOADIRL_MRI_IMAGEQC_12_08_15.csv there is a single scan with the field series_selected == 1
        - If yes choose the one with series_selected == 1
        - If no choose the scan with the best quality

    Args:
        source_dir: path to the ADNI image folder
        clinical_dir: path to the directory with all the clinical data od ADNI
        dest_dir:  path to the output_folder
        subjs_list: subjects list

    Returns: pandas Dataframe containing the path for each fmri

    """
    import os
    from os import path
    from os import walk
    import pandas as pd
    import logging
    from clinica.iotools.converters.adni_to_bids import adni_utils
    from clinica.utils.stream import cprint

    fmri_col = ['Subject_ID', 'VISCODE', 'Visit', 'IMAGEUID', 'Sequence', 'Scan Date', 'LONIUID', 'Scanner',
                'MagStregth', 'Path']

    fmri_df = pd.DataFrame(columns=fmri_col)

    # Load the requested clinical data
    mayo_mri_fmri_path = path.join(clinical_dir, 'MAYOADIRL_MRI_FMRI_09_15_16.csv')
    mayo_mri_imageqc_path = path.join(clinical_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')
    ida_mr_metadata_path = path.join(clinical_dir, 'IDA_MR_Metadata_Listing.csv')

    mayo_mri_fmri = pd.io.parsers.read_csv(mayo_mri_fmri_path, sep=',')
    ida_mr_metadata = pd.io.parsers.read_csv(ida_mr_metadata_path, sep=',')
    mayo_mri_imageqc = pd.io.parsers.read_csv(mayo_mri_imageqc_path, sep=',')

    for subj in subjs_list:
        # print subj
        fmri_subjs_info = mayo_mri_fmri[(mayo_mri_fmri.RID == int(subj[-4:]))]
        # Extract visits available
        visits_list = fmri_subjs_info['VISCODE2'].tolist()
        # Removing duplicates
        visits_list = list(set(visits_list))

        if len(visits_list) != 0:
            for viscode in visits_list:
                visit = ''
                image_path = ''

                fmri_subj = fmri_subjs_info[fmri_subjs_info['VISCODE2'] == viscode]

                if not fmri_subj.empty:

                    # If there are multiple scans for the same session same subject, check what is the one selected for the usage (field 'series_selected') or
                    # choose the one with the best quality
                    if len(fmri_subj) > 1:
                        fmri_imageuid = fmri_subj['IMAGEUID'].tolist()
                        loni_uid_list = ['I' + str(imageuid) for imageuid in fmri_imageuid]
                        images_qc = mayo_mri_imageqc[mayo_mri_imageqc.loni_image.isin(loni_uid_list)]
                        series_selected_values = images_qc['series_selected'].tolist()
                        sum_series_selected = sum(series_selected_values)
                        if sum_series_selected == 1:
                            imageuid_to_select = images_qc[images_qc['series_selected'] > 0]['loni_image'].iloc[
                                0].replace('I', '')
                        else:
                            imageuid_to_select = select_image_qc(fmri_imageuid, images_qc)

                        fmri_subj = fmri_subj[fmri_subj['IMAGEUID'] == int(imageuid_to_select)].iloc[0]
                    else:
                        fmri_subj = fmri_subj.iloc[0]

                    fmri_imageuid = fmri_subj['IMAGEUID']

                    # Discard scans made with non Philips scanner and with a bad quality
                    fmri_metadata = ida_mr_metadata[ida_mr_metadata['IMAGEUID'] == fmri_imageuid]

                    if not fmri_metadata.empty:
                        fmri_metadata = fmri_metadata.iloc[0]

                        if 'Philips' not in fmri_metadata['Scanner']:
                            cprint('No Philips scanner for ' + subj + ' visit ' + viscode + '. Skipped.')
                            continue

                        elif 4 in mayo_mri_imageqc[mayo_mri_imageqc['loni_image'] == 'I' + str(fmri_imageuid)]['series_quality'].values:
                            cprint('Bad scan quality for ' + subj + ' visit ' + viscode + '. Skipped.')
                            continue

                        scan_date = fmri_subj.SCANDATE
                        sequence = adni_utils.replace_sequence_chars(fmri_subj.SERDESC)
                        scanner = fmri_metadata['Scanner']
                        loni_uid = fmri_metadata['LONIUID']
                        visit = fmri_metadata['Visit']
                        mag_strenght = fmri_metadata['MagStrength']

                        # Calculate the path
                        seq_path = path.join(source_dir, str(subj), sequence)
                        for (dirpath, dirnames, filenames) in walk(seq_path):
                            found = False
                            for d in dirnames:
                                if d == 'S' + str(loni_uid):
                                    image_path = path.join(dirpath, d)
                                    # Check if the path exists
                                    if not os.path.isdir(image_path):
                                        cprint('Path not existing for subject ' + subj + ' visit ' + visit)
                                    found = True
                                    break
                            if found:
                                break

                        # The session scmri correspond to the baseline
                        if viscode == 'scmri':
                            viscode = 'bl'
                    else:
                        cprint('Missing visit, sequence, scan date and loniuid for subject ' + subj + ' visit ' + visit)
                        continue

                    row_to_append = pd.DataFrame(
                        [[subj, str(viscode), visit, str(fmri_imageuid), sequence, scan_date, str(loni_uid),
                          scanner, mag_strenght, image_path]], columns=fmri_col)

                    fmri_df = fmri_df.append(row_to_append, ignore_index=True)
                else:
                    logging.info('Missing fMRI for ', subj, 'visit', visit)

    fmri_df.to_csv(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t', index=False)
    return fmri_df


def fmri_paths_to_bids(dest_dir, fmri_paths, mod_to_update=False):
    """
    Convert the fmri extracted from the fmri_paths to BIDS

    Args:
        dest_dir: path to the input directory
        fmri_paths: path to the BIDS directory
        mod_to_add: if True add the fmri only where is missing
        mod_to_update:  if True overwrite (or create if is missing) all the existing fmri

    """
    from multiprocessing import cpu_count, Pool, Value
    from functools import partial

    counter = None

    def init(args):
        ''' store the counter for later use '''
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
    import subprocess
    import os
    import shutil
    from os import path
    from glob import glob

    sess_list = fmri_paths[(fmri_paths['Subject_ID'] == subj)]['VISCODE'].values
    alpha_id = adni_utils.remove_space_and_symbols(subj)
    bids_id = 'sub-ADNI' + alpha_id

    # For each session available, create the folder if doesn't exist and convert the files
    for ses in sess_list:
        with counter.get_lock():
            counter.value += 1
        cprint('[FLAIR] Processing subject ' + str(subj)
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
                                                                 fmri_info['IMAGEUID'].values[0])
                if not os.path.isfile(os.path.join(ses_path, 'func', bids_file_name + '_task-rest_bold.nii.gz')):
                    bids_utils.convert_fmri(dcm_to_convert, path.join(ses_path, 'func'), bids_file_name)
                else:
                    cprint("Images already converted")

                # Delete the temporary folder used for copying fmri with 2 subjects inside the DICOM folder
                adni_utils.remove_tmp_dmc_folder(dest_dir, fmri_info['IMAGEUID'].values[0])


def select_image_qc(id_list, mri_qc_subj):
    """

    Args:
        id_list:
        mri_qc_subj:

    Returns:

    """
    import numpy as np

    if len(id_list) == 0:
        return None

    selected_image = None
    image_ids = ['I' + str(imageuid) for imageuid in id_list]
    int_ids = [int(imageuid) for imageuid in id_list]
    images_qc = mri_qc_subj[mri_qc_subj.loni_image.isin(image_ids)]

    if images_qc.shape[0] < 1:
        return min(int_ids)

    if np.sum(images_qc.series_selected) == 1:
        selected_image = images_qc[images_qc.series_selected == 1].iloc[0].loni_image[1:]
    else:
        images_not_rejected = images_qc[images_qc.series_quality < 4]

        if images_not_rejected.shape[0] < 1:

            # There are no images that passed the qc
            # so we'll try to see if there are other images without qc,
            # otherwise return None
            qc_ids = set([int(qc_id[1:]) for qc_id in images_qc.loni_image.unique()])
            no_qc_ids = list(set(int_ids) - qc_ids)

            if len(no_qc_ids) == 0:
                return None
            else:
                return min(no_qc_ids)

        series_quality = [q if q > 0 else 4 for q in list(images_not_rejected.series_quality)]
        best_q = np.amin(series_quality)

        images_best_qc = images_not_rejected[images_not_rejected.series_quality == best_q]
        if images_best_qc.shape[0] == 1:
            selected_image = images_best_qc.iloc[0].loni_image[1:]
        else:
            selected_image = min(int_ids)

    return int(selected_image)


def compute_fmri_path_refactoring(source_dir, clinical_dir, dest_dir, subjs_list):
    """
    Compute the paths to fmri images.

    The fmri images to convert into BIDS are chosen in the following way:
        - Extract the list of subjects from MAYOADIRL_MRI_FMRI_09_15_16.csv
        - Select the only the scans that came from PHILIPS machine (field Scanner from IDA_MR_Metadata_Listing.csv)
        - Discard all the subjects with column  series_quality = 4  (4 means that the scan is not usable) in MAYOADIRL_MRI_IMAGEQC_12_08_15.csv

    In case of multiple scans for the same session, same date the one to convert is chosen with the following criteria:
        - Check if in the file MAYOADIRL_MRI_IMAGEQC_12_08_15.csv there is a single scan with the field series_selected == 1
        - If yes choose the one with series_selected == 1
        - If no choose the scan with the best quality

    Args:
        source_dir: path to the ADNI image folder
        clinical_dir: path to the directory with all the clinical data od ADNI
        dest_dir:  path to the output_folder
        subjs_list: subjects list

    Returns: pandas Dataframe containing the path for each fmri

    """
    import os
    from os import path
    from os import walk
    import pandas as pd
    import logging
    from clinica.iotools.converters.adni_to_bids import adni_utils
    from clinica.utils.stream import cprint

    fmri_col = ['Subject_ID', 'VISCODE', 'Visit', 'IMAGEUID', 'Sequence', 'Scan Date', 'LONIUID', 'Scanner',
                'MagStregth', 'Path']

    fmri_df = pd.DataFrame(columns=fmri_col)

    # Load the requested clinical data
    mayo_mri_fmri_path = path.join(clinical_dir, 'MAYOADIRL_MRI_FMRI_09_15_16.csv')
    mayo_mri_imageqc_path = path.join(clinical_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    mayo_mri_fmri = pd.io.parsers.read_csv(mayo_mri_fmri_path, sep=',')
    mayo_mri_imageqc = pd.io.parsers.read_csv(mayo_mri_imageqc_path, sep=',')

    for subj in subjs_list:
        # print subj
        fmri_subjs_info = mayo_mri_fmri[(mayo_mri_fmri.RID == int(subj[-4:]))]
        # Extract visits available
        visits_list = fmri_subjs_info['VISCODE2'].tolist()
        # Removing duplicates
        visits_list = list(set(visits_list))

        if len(visits_list) != 0:
            for viscode in visits_list:
                visit = ''
                image_path = ''

                fmri_subj = fmri_subjs_info[fmri_subjs_info['VISCODE2'] == viscode]

                if not fmri_subj.empty:

                    # If there are multiple scans for the same session same subject, check what is the one selected for the usage (field 'series_selected') or
                    # choose the one with the best quality
                    if len(fmri_subj) > 1:
                        fmri_imageuid = fmri_subj['IMAGEUID'].tolist()
                        loni_uid_list = ['I' + str(imageuid) for imageuid in fmri_imageuid]
                        images_qc = mayo_mri_imageqc[mayo_mri_imageqc.loni_image.isin(loni_uid_list)]
                        series_selected_values = images_qc['series_selected'].tolist()
                        sum_series_selected = sum(series_selected_values)
                        if sum_series_selected == 1:
                            imageuid_to_select = images_qc[images_qc['series_selected'] > 0]['loni_image'].iloc[
                                0].replace('I', '')
                        else:
                            imageuid_to_select = select_image_qc(fmri_imageuid, images_qc)

                        fmri_subj = fmri_subj[fmri_subj['IMAGEUID'] == int(imageuid_to_select)].iloc[0]
                    else:
                        fmri_subj = fmri_subj.iloc[0]

                    fmri_imageuid = fmri_subj['IMAGEUID']

                    # Discard scans made with non Philips scanner and with a bad quality
                    fmri_metadata = ida_mr_metadata[ida_mr_metadata['IMAGEUID'] == fmri_imageuid]

                    if not fmri_metadata.empty:
                        fmri_metadata = fmri_metadata.iloc[0]

                        if 'Philips' not in fmri_metadata['Scanner']:
                            cprint('No Philips scanner for ' + subj + ' visit ' + viscode + '. Skipped.')
                            continue

                        elif 4 in mayo_mri_imageqc[mayo_mri_imageqc['loni_image'] == 'I' + str(fmri_imageuid)]['series_quality'].values:
                            cprint('Bad scan quality for ' + subj + ' visit ' + viscode + '. Skipped.')
                            continue

                        scan_date = fmri_subj.SCANDATE
                        sequence = adni_utils.replace_sequence_chars(fmri_subj.SERDESC)
                        scanner = fmri_metadata['Scanner']
                        loni_uid = fmri_metadata['LONIUID']
                        visit = fmri_metadata['Visit']
                        mag_strenght = fmri_metadata['MagStrength']

                        # Calculate the path
                        seq_path = path.join(source_dir, str(subj), sequence)
                        for (dirpath, dirnames, filenames) in walk(seq_path):
                            found = False
                            for d in dirnames:
                                if d == 'S' + str(loni_uid):
                                    image_path = path.join(dirpath, d)
                                    # Check if the path exists
                                    if not os.path.isdir(image_path):
                                        cprint('Path not existing for subject ' + subj + ' visit ' + visit)
                                    found = True
                                    break
                            if found:
                                break

                        # The session scmri correspond to the baseline
                        if viscode == 'scmri':
                            viscode = 'bl'
                    else:
                        cprint('Missing visit, sequence, scan date and loniuid for subject ' + subj + ' visit ' + visit)
                        continue

                    row_to_append = pd.DataFrame(
                        [[subj, str(viscode), visit, str(fmri_imageuid), sequence, scan_date, str(loni_uid),
                          scanner, mag_strenght, image_path]], columns=fmri_col)

                    fmri_df = fmri_df.append(row_to_append, ignore_index=True)
                else:
                    logging.info('Missing fMRI for ', subj, 'visit', visit)

    fmri_df.to_csv(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t', index=False)
    return fmri_df
