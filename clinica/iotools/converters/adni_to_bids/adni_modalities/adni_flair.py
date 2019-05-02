# coding: utf-8

"""
 Module for converting FLAIR of ADNI
"""

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_flair(source_dir, csv_dir, dest_dir, subjs_list=None):
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

    cprint('Calculating paths of FLAIR images. Output will be stored in '
           + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_flair_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of FLAIR images found. Exporting images into BIDS ...')
    flair_paths_to_bids(images, dest_dir)
    cprint(Fore.GREEN + 'FLAIR conversion done.' + Fore.RESET)
    pass


def compute_flair_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute paths to FLAIR images to convert to BIDS

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination directory
        subjs_list: subjects list

    Returns:
        images: pandas dataframe that contains the path to all the FLAIR images to convert

    """
    import pandas as pd
    from os import path, walk, mkdir
    from clinica.utils.stream import cprint

    flair_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                    'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner', 'Enhanced']

    flair_df = pd.DataFrame(columns=flair_col_df)
    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')

    ida_meta_path = path.join(csv_dir, 'MRILIST.csv')
    mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
    ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')

    ida_meta = ida_meta[ida_meta.SEQUENCE.map(lambda x: x.lower().find('flair') > -1)]
    mri_qc = pd.io.parsers.read_csv(mri_qc_path, sep=',')
    mri_qc = mri_qc[mri_qc.series_type == 'AFL']

    for subj in subjs_list:
        # print 'Computing path for subj', subj

        adnimerge_subj = adni_merge[adni_merge.PTID == subj]

        # Sort the values by examination date
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')
        ida_meta_subj = ida_meta[ida_meta.SUBJECT == subj]
        ida_meta_subj = ida_meta_subj.sort_values('SCANDATE')

        mri_qc_subj = mri_qc[mri_qc.RID == int(subj[-4:])]
        visits = visits_to_timepoints_flair(subj, ida_meta_subj, adnimerge_subj)
        keys = visits.keys()
        # What is supposed to do the line below ????
        # keys.sort()

        for visit_info in visits.keys():
            visit_str = visits[visit_info]
            visit_ida_meta = ida_meta_subj[ida_meta_subj.VISIT == visit_str]
            axial_ida_meta = visit_ida_meta[visit_ida_meta.SEQUENCE.map(lambda x: x.lower().find('enhanced') < 0)]
            axial = flair_image(subj, visit_info[0], visits[visit_info], axial_ida_meta, mri_qc_subj, False)
            row_to_append = pd.DataFrame(axial, index=['i', ])
            flair_df = flair_df.append(row_to_append, ignore_index=True)

    images = flair_df
    is_dicom = []
    nifti_paths = []
    count = 0

    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), str(image.Sequence))
        count += 1
        series_path = ''
        s = 'S' + str(image.Series_ID)
        for (dirpath, dirnames, filenames) in walk(seq_path):
            found = False
            for d in dirnames:
                if d == s:
                    series_path = path.join(dirpath, d)
                    found = True
                    break
            if found:
                break

        nifti_path = series_path
        dicom = True
        is_dicom.append(dicom)
        nifti_paths.append(nifti_path)

    images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
    images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

    # Drop all the lines that have the Path section empty
    # images = images.drop(images[images.Path == ''].index)
    # Store the paths inside a file called t1_paths inside the input directory

    flair_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(flair_tsv_path):
        mkdir(flair_tsv_path)

    cprint('\tDone! Saving the results into' + path.join(flair_tsv_path, 'flair_paths.tsv'))
    images.to_csv(path.join(flair_tsv_path, 'flair_paths.tsv'),
                  sep='\t',
                  index=False)
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
        # For the same subject, same session there could be multiple flar with different acq label
        for j in range(len(flair_info)):
            flair_subj = flair_info.iloc[j]
            # TODO For now in CLINICA we ignore Enhanced FLAIR.
            if flair_subj['Enhanced']:
                continue
            if type(flair_subj['Path']) != float and flair_subj['Path'] != '':
                if not os.path.exists(path.join(ses_path, 'FLAIR')):
                    os.mkdir(path.join(ses_path, 'FLAIR'))
                flair_path = flair_subj['Path']

                bids_name = bids_file_name + '_FLAIR'
                # bids.dcm_to_nii(dwi_path, path.join(ses_path, 'dwi'), bids_name)

                bids_dest_dir = path.join(ses_path, 'FLAIR')

                if not os.path.exists(bids_dest_dir):
                    os.mkdir(dest_dir)
                command = 'dcm2niix -b y -z y -o ' + bids_dest_dir + ' -f ' + bids_name + ' ' + flair_path
                subprocess.run(command,
                               shell=True,
                               stderr=subprocess.DEVNULL,
                               stdout=subprocess.DEVNULL)

                # If dcm2niix didn't work use dcm2nii
                # print path.join(dest_dir, bids_name + '.nii.gz')
                if not os.path.exists(
                        path.join(bids_dest_dir, bids_name + '.nii.gz')):
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
        if best_q == 4:
            best_q = -1
        images_best_qc = images_not_rejected[images_not_rejected.series_quality == best_q]
        if images_best_qc.shape[0] == 1:
            selected_image = images_best_qc.iloc[0].loni_image[1:]
        else:
            best_ids = [int(x[1:]) for x in images_best_qc.loni_image.unique()]
            selected_image = min(best_ids)

    return int(selected_image)


def visits_to_timepoints_flair(subject, ida_meta_subj, adnimerge_subj):
    """

    Args:
        subject:
        ida_meta_subj:
        adnimerge_subj:

    Returns:

    """
    from datetime import datetime
    from clinica.iotools.converters.adni_to_bids.adni_utils import days_between
    from clinica.utils.stream import cprint

    visits = dict()
    unique_visits = list(ida_meta_subj.VISIT.unique())
    pending_timepoints = []

    # We try to obtain the corresponding image Visit for a given VISCODE
    for adni_row in adnimerge_subj.iterrows():
        visit = adni_row[1]
        if visit.ORIGPROT == 'ADNI2':
            if visit.VISCODE == 'bl':
                preferred_visit_name = 'ADNI2 Screening MRI-New Pt'
            elif visit.VISCODE == 'm03':
                preferred_visit_name = 'ADNI2 Month 3 MRI-New Pt'
            elif visit.VISCODE == 'm06':
                preferred_visit_name = 'ADNI2 Month 6-New Pt'
            else:
                year = str(int(visit.VISCODE[1:]) / 12)
                preferred_visit_name = 'ADNI2 Year ' + year + ' Visit'
        else:
            if visit.VISCODE == 'bl':
                if visit.ORIGPROT == 'ADNI1':
                    preferred_visit_name = 'ADNI Screening'
                else:  # ADNIGO
                    preferred_visit_name = 'ADNIGO Screening MRI'
            elif visit.VISCODE == 'm03':  # Only for ADNIGO Month 3
                preferred_visit_name = 'ADNIGO Month 3 MRI'
            else:
                month = int(visit.VISCODE[1:])
                if month < 54:
                    preferred_visit_name = 'ADNI1/GO Month ' + str(month)
                else:
                    preferred_visit_name = 'ADNIGO Month ' + str(month)

        if preferred_visit_name in unique_visits:
            key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
            if key_preferred_visit not in visits.keys():
                visits[key_preferred_visit] = preferred_visit_name
            elif visits[key_preferred_visit] != preferred_visit_name:
                cprint('Multiple visits for one timepoint!')
                cprint(subject)
                cprint(key_preferred_visit)
                cprint(visits[key_preferred_visit])
                cprint(visit)
            unique_visits.remove(preferred_visit_name)
            continue

        pending_timepoints.append(visit)

    # Then for images.Visit non matching the expected labels we find the closest date in visits list
    for visit in unique_visits:
        image = (ida_meta_subj[ida_meta_subj.VISIT == visit]).iloc[0]
        min_db = 100000
        min_db2 = 0
        min_visit = None
        min_visit2 = None

        for timepoint in pending_timepoints:
            db = days_between(image['SCANDATE'], timepoint.EXAMDATE)
            if db < min_db:
                min_db2 = min_db
                min_visit2 = min_visit

                min_db = db
                min_visit = timepoint

        if min_visit is None:
            cprint('No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.VISIT)
            cprint(image)
            continue

        if min_visit2 is not None and min_db > 90:
            cprint('More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.VISIT + ' on ' + image.SCANDATE)
            cprint('Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                min_db) + ' days)')
            cprint('Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                min_db2) + ' days)')

            # If image is too close to the date between two visits we prefer the earlier visit
            if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                    > datetime.strptime(image.SCANDATE, "%Y-%m-%d")
                    > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                dif = days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                if abs((dif / 2.0) - min_db) < 30:
                    min_visit = min_visit2
            cprint('We prefer ' + min_visit.VISCODE)

        key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
        if key_min_visit not in visits.keys():
            visits[key_min_visit] = image.VISIT
        elif visits[key_min_visit] != image.VISIT:
            cprint('Multiple visits for one timepoint!')

    return visits


def flair_image(subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):
    """

    Args:
        subject_id:
        timepoint:
        visit_str:
        ida_meta_scans:
        mri_qc_subj:
        enhanced:

    Returns:

    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars

    sel_image = select_image_qc(list(ida_meta_scans.IMAGEUID), mri_qc_subj)
    if sel_image is None:
        return None

    sel_scan = ida_meta_scans[ida_meta_scans.IMAGEUID == sel_image].iloc[0]
    sequence = sel_scan.SEQUENCE
    sequence = replace_sequence_chars(sequence)
    image_dict = {'Subject_ID': subject_id,
                  'VISCODE': timepoint,
                  'Visit': visit_str,
                  'Sequence': sequence,
                  'Scan_Date': sel_scan['SCANDATE'],
                  'Study_ID': str(int(sel_scan.STUDYID)),
                  'Series_ID': str(int(sel_scan.SERIESID)),
                  'Image_ID': str(int(sel_scan.IMAGEUID)),
                  'Field_Strength': sel_scan.MAGSTRENGTH,
                  'Enhanced': enhanced}

    return image_dict
