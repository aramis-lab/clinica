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

    cprint('Calculating paths of DWI images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_dwi_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of DWI images found. Exporting images into BIDS ...')
    dwi_paths_to_bids(images, dest_dir)
    cprint(Fore.GREEN + 'DWI conversion done.' + Fore.RESET)


def compute_dwi_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute paths to DWI images to convert to BIDS

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination directory
        subjs_list: subjects list

    Returns:
        images: pandas dataframe that contains the path to all the dwi images to convert

    """
    import pandas as pd
    import operator
    from os import path, walk, mkdir
    from functools import reduce

    dwi_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                  'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner', 'Enhanced']

    dwi_df = pd.DataFrame(columns=dwi_col_df)

    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')

    if path.exists(path.join(csv_dir, 'IDA_MR_Metadata_Listing.csv')):
        new_download = False
    else:
        new_download = True

    # If IDA_MR_Metadata_Listing.csv has been added manually, that will cause a
    # bug
    new_download = True

    if new_download:
        ida_meta_path = path.join(csv_dir, 'MRILIST.csv')
    else:
        ida_meta_path = path.join(csv_dir, 'IDA_MR_Metadata_Listing.csv')
    mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
    ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')

    if new_download:
        ida_meta = ida_meta[ida_meta.SEQUENCE.map(lambda x: x.lower().find('dti') > -1)]
    else:
        ida_meta = ida_meta[ida_meta.Sequence.map(lambda x: x.lower().find('dti') > -1)]

    mri_qc = pd.io.parsers.read_csv(mri_qc_path, sep=',')
    mri_qc = mri_qc[mri_qc.series_type == 'DTI']

    for subj in subjs_list:
        # print 'Computing path for subj', subj

        adnimerge_subj = adni_merge[adni_merge.PTID == subj]

        # Sort the values by examination date
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')
        if new_download:
            ida_meta_subj = ida_meta[ida_meta.SUBJECT == subj]
            ida_meta_subj = ida_meta_subj.sort_values('SCANDATE')
        else:
            ida_meta_subj = ida_meta[ida_meta.Subject == subj]
            ida_meta_subj = ida_meta_subj.sort_values('Scan Date')
        mri_qc_subj = mri_qc[mri_qc.RID == int(subj[-4:])]

        if new_download:
            visits = visits_to_timepoints_dwi_refactoring(subj, ida_meta_subj, adnimerge_subj)
        else:
            visits = visits_to_timepoints_dwi(subj, ida_meta_subj, adnimerge_subj)

        keys = visits.keys()
        # What is the purpose of the following line ?
        # keys.sort()

        for visit_info in visits.keys():
            visit_str = visits[visit_info]

            if new_download:
                visit_ida_meta = ida_meta_subj[ida_meta_subj.VISIT == visit_str]
                axial_ida_meta = visit_ida_meta[visit_ida_meta.SEQUENCE.map(lambda x: x.lower().find('enhanced') < 0)]
                enhanced_ida_meta = visit_ida_meta[visit_ida_meta.SEQUENCE.map(lambda x: x.lower().find('enhanced') > -1)]
                axial = dwi_image_refactoring(subj, visit_info[0], visits[visit_info], axial_ida_meta, mri_qc_subj, False)
                enhanced = dwi_image_refactoring(subj, visit_info[0], visits[visit_info], enhanced_ida_meta, mri_qc_subj, True)
            else:
                visit_ida_meta = ida_meta_subj[ida_meta_subj.Visit == visit_str]
                axial_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') < 0)]
                enhanced_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') > -1)]
                axial = dwi_image(subj, visit_info[0], visits[visit_info], axial_ida_meta, mri_qc_subj, False)
                enhanced = dwi_image(subj, visit_info[0], visits[visit_info], enhanced_ida_meta, mri_qc_subj, True)

            if axial is not None:
                row_to_append = pd.DataFrame(axial, index=['i', ])
                dwi_df = dwi_df.append(row_to_append, ignore_index=True)

            if enhanced is not None:
                row_to_append = pd.DataFrame(enhanced, index=['i', ])
                dwi_df = dwi_df.append(row_to_append, ignore_index=True)

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
                         ('016_S_2007', 'm06')]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((dwi_df.Enhanced is False)
                             & (dwi_df.Subject_ID == conv_error[0])
                             & (dwi_df.VISCODE == conv_error[1]))

    indices_to_remove = dwi_df.index[reduce(operator.or_, error_indices, False)]
    dwi_df.drop(indices_to_remove, inplace=True)

    images = dwi_df
    is_dicom = []
    nifti_paths = []
    count = 0

    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)
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


def generate_subject_files(subj, images, dest_dir, mod_to_update):
    import clinica.iotools.bids_utils as bids
    import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
    from clinica.utils.stream import cprint
    import subprocess
    import os
    import shutil
    from os import path
    from glob import glob

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

        if not os.path.exists(ses_path):
            os.mkdir(ses_path)

        dwi_info = images[
            (images['Subject_ID'] == subj) & (images['VISCODE'] == ses)]

        # For the same subject, same session there could be multiple dwi with different acq label
        for j in range(len(dwi_info)):
            dwi_subj = dwi_info.iloc[j]
            # TODO For now in CLINICA we ignore Enhanced DWI.
            if dwi_subj['Enhanced']:
                continue
            if type(dwi_subj['Path']) != float and dwi_subj['Path'] != '':
                if not os.path.exists(path.join(ses_path, 'dwi')):
                    os.mkdir(path.join(ses_path, 'dwi'))
                dwi_path = dwi_subj['Path']

                bids_name = bids_file_name + '_acq-' + ('axialEnhanced'
                                                        if dwi_subj['Enhanced'] else 'axial') + '_dwi'

                bids_dest_dir = path.join(ses_path, 'dwi')

                if not os.path.exists(bids_dest_dir):
                    os.mkdir(dest_dir)
                command = 'dcm2niix -b n -z y -o ' + bids_dest_dir + ' -f ' + bids_name + ' ' + dwi_path
                subprocess.run(command,
                               shell=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)

                # If dcm2niix didn't work use dcm2nii
                # print path.join(dest_dir, bids_name + '.nii.gz')
                if not os.path.exists(path.join(bids_dest_dir,
                                                bids_name + '.nii.gz')) or not os.path.exists(
                        path.join(bids_dest_dir,
                                  bids_name + '.bvec') or not os.path.exists(
                            path.join(bids_dest_dir,
                                      bids_name + '.bval'))):
                    cprint('\tConversion with dcm2niix failed, trying with dcm2nii')

                    # Find all the files eventually created by dcm2niix and remove them
                    dwi_dcm2niix = glob(path.join(bids_dest_dir, bids_name + '*'))

                    for d in dwi_dcm2niix:
                        # print 'Removing the old', d
                        os.remove(d)

                    command = 'dcm2nii -a n -d n -e n -i y -g y -p n -m n -r n -x n -o ' + bids_dest_dir + ' ' + dwi_path
                    subprocess.run(command,
                                   shell=True,
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                    nii_file = path.join(bids_dest_dir,
                                         subj.replace('_', '') + '.nii.gz')
                    bvec_file = path.join(bids_dest_dir,
                                          subj.replace('_', '') + '.bvec')
                    bval_file = path.join(bids_dest_dir,
                                          subj.replace('_', '') + '.bval')

                    if os.path.exists(bvec_file) and os.path.exists(
                            bval_file):
                        os.rename(bvec_file, path.join(bids_dest_dir,
                                                       bids_name + '.bvec'))
                        os.rename(bval_file, path.join(bids_dest_dir,
                                                       bids_name + '.bval'))
                    else:
                        cprint('WARNING: bvec and bval not generated by dcm2nii'
                               + ' for subject ' + subj + ' and session ' + ses)

                    if os.path.exists(nii_file):
                        os.rename(nii_file, path.join(bids_dest_dir,
                                                      bids_name + '.nii.gz'))
                    else:
                        cprint('WARNING: CONVERSION FAILED...'
                               + ' for subject ' + subj + ' and session ' + ses)


def dwi_image(subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):
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

    sequence = sel_scan.Sequence
    sequence = replace_sequence_chars(sequence)

    image_dict = {'Subject_ID': subject_id,
                  'VISCODE': timepoint,
                  'Visit': visit_str,
                  'Sequence': sequence,
                  'Scan_Date': sel_scan['Scan Date'],
                  'Study_ID': str(int(sel_scan.LONISID)),
                  'Series_ID': str(int(sel_scan.LONIUID)),
                  'Image_ID': str(int(sel_scan.IMAGEUID)),
                  'Field_Strength': sel_scan.MagStrength,
                  'Scanner': sel_scan.Scanner,
                  'Enhanced': enhanced}

    return image_dict


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


def visits_to_timepoints_dwi(subject, ida_meta_subj, adnimerge_subj):
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
    unique_visits = list(ida_meta_subj.Visit.unique())
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
        image = (ida_meta_subj[ida_meta_subj.Visit == visit]).iloc[0]
        min_db = 100000
        min_db2 = 0
        min_visit = None
        min_visit2 = None

        for timepoint in pending_timepoints:
            db = days_between(image['Scan Date'], timepoint.EXAMDATE)
            if db < min_db:
                min_db2 = min_db
                min_visit2 = min_visit

                min_db = db
                min_visit = timepoint

        if min_visit is None:
            cprint('No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit)
            cprint(image)
            continue

        if min_visit2 is not None and min_db > 90:
            cprint('More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit + ' on ' + image.ScanDate)
            cprint('Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                min_db) + ' days)')
            cprint('Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                min_db2) + ' days)')

            # If image is too close to the date between two visits we prefer the earlier visit
            if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                    > datetime.strptime(image.ScanDate, "%Y-%m-%d")
                    > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                dif = days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                if abs((dif / 2.0) - min_db) < 30:
                    min_visit = min_visit2

            cprint('We prefer ' + min_visit.VISCODE)

        key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
        if key_min_visit not in visits.keys():
            visits[key_min_visit] = image.Visit
        elif visits[key_min_visit] != image.Visit:
            cprint('Multiple visits for one timepoint!')
            cprint(subject)
            cprint(key_min_visit)
            cprint(visits[key_min_visit])
            cprint(image.Visit)

    return visits


def visits_to_timepoints_dwi_refactoring(subject, ida_meta_subj, adnimerge_subj):
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
            # cprint(subject)
            # cprint(key_min_visit)
            # cprint(visits[key_min_visit])
            # cprint(image.Visit)

    return visits


def dwi_image_refactoring(subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):
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
