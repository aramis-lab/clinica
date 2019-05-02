# coding: utf-8

"""
 Module for converting T1 of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_t1(source_dir, csv_dir, dest_dir, subjs_list=None):
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
    from clinica.iotools.converters.adni_to_bids.adni_utils import t1_pet_paths_to_bids
    from colorama import Fore

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
        subjs_list = list(adni_merge.PTID.unique())
    # new_download: indicates the version of csv Downloaded
    if path.isfile(path.join(csv_dir, 'IDA_MR_Metadata_Listing.csv')):
        new_download = False
    else:
        new_download = True
    if new_download is True:
        if '127_S_5200' in subjs_list:
            subjs_list.remove('127_S_5200')
        if '027_S_5083' in subjs_list:
            subjs_list.remove('027_S_5083')
        if '941_S_4365' in subjs_list:
            subjs_list.remove('941_S_4365')
        if '037_S_4028' in subjs_list:
            subjs_list.remove('037_S_4028')

    cprint('Calculating paths of T1 images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list, new_download)
    cprint('Paths of T1 images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 't1')
    cprint(Fore.GREEN + 'T1 conversion done.' + Fore.RESET)


def compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list, new_download):
    """

    Compute paths to t1 images of ADNI.

    :param source_dir:
    :param csv_dir:
    :param dest_dir:
    :param subjs_list:
    :return: a pandas dataframe
    """

    import pandas as pd
    import operator
    from os import path, walk, mkdir
    from functools import reduce

    t1_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                 'Study_ID', 'Field_Strength', 'Series_ID', 'Image_ID', 'Original']

    t1_df = pd.DataFrame(columns=t1_col_df)
    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
    # adni_screening_path = path.join(clinical_dir, 'ADNI_ScreeningList_8_22_12.csv')
    if new_download is False:
        ida_meta_path = path.join(csv_dir, 'IDA_MR_Metadata_Listing.csv')
        ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',', low_memory=False)

    mprage_meta_path = path.join(csv_dir, 'MPRAGEMETA.csv')
    mri_quality_path = path.join(csv_dir, 'MRIQUALITY.csv')
    mayo_mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',', low_memory=False)
    mprage_meta = pd.io.parsers.read_csv(mprage_meta_path, sep=',', low_memory=False)
    mri_quality = pd.io.parsers.read_csv(mri_quality_path, sep=',', low_memory=False)
    mayo_mri_qc = pd.io.parsers.read_csv(mayo_mri_qc_path, sep=',', low_memory=False)
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'T1']

    for subj in subjs_list:
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        # Sort the values by examination date
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mprage_meta_subj = mprage_meta[mprage_meta.SubjectID == subj]
        mprage_meta_subj = mprage_meta_subj.sort_values('ScanDate')

        mri_quality_subj = mri_quality[mri_quality.RID == int(subj[-4:])]
        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
        visits = visits_to_timepoints_t1(subj, mprage_meta_subj, adnimerge_subj)

        keys = list(visits.keys())
        keys.sort()
        for visit_info in visits.keys():
            if visit_info[1] == 'ADNI1':
                if not new_download:
                    ida_meta_subj = ida_meta[ida_meta.Subject == subj]
                    image_dict = adni1_image(
                            subj, visit_info[0], visits[visit_info],
                            mprage_meta_subj, ida_meta_subj,
                            mri_quality_subj, mayo_mri_qc_subj)
                else:
                    image_dict = adni1_image_refactoring(csv_dir, adni_merge, subj, visit_info[0], visits[visit_info],
                                                         mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj)
            elif visit_info[1] == 'ADNIGO':
                if not new_download:
                    image_dict = adnigo_image(
                            subj, visit_info[0], visits[visit_info],
                            mprage_meta_subj, ida_meta_subj, mri_quality_subj,
                            mayo_mri_qc_subj, visit_info[2])
                else:
                    image_dict = adnigo_image_refactoring(
                            csv_dir, adni_merge, subj, visit_info[0],
                            visits[visit_info], mprage_meta_subj,
                            mri_quality_subj, mayo_mri_qc_subj, visit_info[2])

            else:  # ADNI2
                image_dict = adni2_image(
                        subj, visit_info[0], visits[visit_info],
                        mprage_meta_subj_orig, mayo_mri_qc_subj)

            if image_dict is None:
                image_dict = {'Subject_ID': subj,
                              'VISCODE': visit_info[0],
                              'Visit': visits[visit_info],
                              'Sequence': '',
                              'Scan_Date': '',
                              'Study_ID': '',
                              'Series_ID': '',
                              'Image_ID': '',
                              'Field_Strength': '',
                              'Original': True}

            row_to_append = pd.DataFrame(image_dict, index=['i', ])
            t1_df = t1_df.append(row_to_append, ignore_index=True)

    # Exceptions
    # ==========
    # Subject 141_S_0726
    if not new_download:
        # conversion errors known
        subj_ind = (t1_df.Subject_ID == '141_S_0726') & (t1_df.VISCODE == 'bl')
        t1_df.loc[subj_ind, 'Sequence'] = 'MPR-R__GradWarp__B1_Correction'
        t1_df.loc[subj_ind, 'Series_ID'] = '18355'
        t1_df.loc[subj_ind, 'Image_ID'] = '94827'

        conversion_errors = [  # Eq_1
                         ('031_S_0830', 'm48'),
                         ('100_S_0995', 'm18'),
                         ('031_S_0867', 'm48'),
                         ('100_S_0892', 'm18'),
                         # Empty folders
                         ('029_S_0845', 'm24'),
                         ('094_S_1267', 'm24'),
                         ('029_S_0843', 'm24'),
                         ('027_S_0307', 'm48'),
                         ('057_S_1269', 'm24'),
                         ('036_S_4899', 'm03')]

        error_indices = []
        for conv_error in conversion_errors:
            error_indices.append((t1_df.Subject_ID == conv_error[0])
                                 & (t1_df.VISCODE == conv_error[1]))

        indices_to_remove = t1_df.index[reduce(operator.or_, error_indices, False)]
        t1_df.drop(indices_to_remove, inplace=True)

    images = t1_df
    is_dicom = []
    nifti_paths = []
    count = 0

    for row in images.iterrows():

        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)

        count += 1
        series_path = ''
        for (dirpath, dirnames, filenames) in walk(seq_path):
            found = False
            for d in dirnames:
                if d == 'S' + str(image.Series_ID):
                    series_path = path.join(dirpath, d)
                    found = True
                    break
            if found:
                break

        nifti_path = series_path
        dicom = True

        for (dirpath, dirnames, filenames) in walk(series_path):
            for f in filenames:
                if f.endswith(".nii"):
                    dicom = False
                    nifti_path = path.join(dirpath, f)
                    break

        is_dicom.append(dicom)
        nifti_paths.append(nifti_path)

    images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
    images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

    if new_download:
        images = images.drop_duplicates()
        for j in images.Path:
            if not path.exists(j):
                images.Path[images.Path == j] = ''

        images = images[images.Path != '']

    # Store the paths inside a file called conversion_info inside the input directory
    t1_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(t1_tsv_path):
        mkdir(t1_tsv_path)
    images.to_csv(path.join(t1_tsv_path, 't1_paths.tsv'), sep='\t', index=False)

    return images


def adni1_image(subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj, mri_quality_subj, mayo_mri_qc_subj):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint
    # Get the preferred scan (image series that has been Scaled)
    filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                       & (mprage_meta_subj.Visit == visit_str)
                                       & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]

    # If there is not a preferred image we use ADNI2 processing
    # (get the best qc if available, otherwise the original) preferring 1.5T images
    if filtered_mprage.shape[0] < 1:
        mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
        return adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, mayo_mri_qc_subj, preferred_field_strength=1.5)

    filtered_mprage_mag = filtered_mprage
    if len(filtered_mprage.MagStrength.unique()) > 1:
        filtered_mprage_mag = filtered_mprage[filtered_mprage.MagStrength == 1.5]  # Select 1.5T images

    scan = filtered_mprage_mag.iloc[0]
    series_id = scan.SeriesID

    qc_passed = True
    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # print 'QC found but NOT passed'
        # print 'Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID)
        mprage_meta_subj_alt = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                                & (mprage_meta_subj.Visit == visit_str)
                                                & (mprage_meta_subj.SeriesID != series_id)]

        qc_prev_sequence = scan.Sequence
        scan = mprage_meta_subj_alt.iloc[0]
        series_id = scan.SeriesID
        qc_passed = False

    # TODO replace by the function to find if the scanner is_philips

    filtered_scan = ida_meta_subj[ida_meta_subj.LONIUID == series_id]

    if filtered_scan.shape[0] < 1:
        # If no IDA_META for 1.5T try for 3T
        filtered_mprage_mag = filtered_mprage[filtered_mprage.MagStrength == 3.0]
        scan = filtered_mprage_mag.iloc[0]
        series_id = scan.SeriesID

        filtered_scan = ida_meta_subj[ida_meta_subj.LONIUID == series_id]

        if filtered_scan.shape[0] < 1:
            # TODO - LOG THIS
            cprint('NO IDA Meta: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
            return None
    ida_scan = filtered_scan.iloc[0]

    original = True

    # TODO replace the condition with if is_philips
    if ida_scan.Scanner.find('Philips') > -1:

        scan = (mprage_meta_subj[
                    (mprage_meta_subj['Orig/Proc'] == 'Original') & (mprage_meta_subj.SeriesID == series_id)]).iloc[
            0]
        sequence = scan.Sequence

    else:  # scan already selected above
        sequence = scan.Sequence[:scan.Sequence.find('N3') - 2]
        original = False

    if not qc_passed:
        if scan.Sequence == 'MP-RAGE':
            original_img_seq = 'MPR'
        else:  # 'MP-RAGE REPEAT'
            original_img_seq = 'MPR-R'

        processing_seq = qc_prev_sequence[qc_prev_sequence.find(';'):qc_prev_sequence.find('N3') - 2]
        sequence = original_img_seq + processing_seq
        # print sequence

    sequence = replace_sequence_chars(sequence)

    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # TODO - LOG THIS
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))

    return {'Subject_ID': subject_id,
            'VISCODE': timepoint,
            'Visit': visit_str,
            'Sequence': sequence,
            'Scan_Date': scan.ScanDate,
            'Study_ID': str(scan.StudyID),
            'Series_ID': str(scan.SeriesID),
            'Image_ID': str(scan.ImageUID),
            'Field_Strength': scan.MagStrength,
            'Original': original}


def adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, mayo_mri_qc_subj, preferred_field_strength=3.0):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint

    cond_mprage = ((mprage_meta_subj_orig.Visit == visit_str)
                   & mprage_meta_subj_orig.Sequence.map(
                       lambda x: ((x.lower().find('mprage') > -1)
                                  | (x.lower().find('mp-rage') > -1)
                                  | (x.lower().find('mp rage') > -1))
                       & (x.find('2') < 0)
                       )
                   )

    cond_spgr = ((mprage_meta_subj_orig.Visit == visit_str) & mprage_meta_subj_orig.Sequence.map(
        lambda x: (x.lower().find('spgr') > -1) & (x.lower().find('acc') < 0)))

    filtered_scan = mprage_meta_subj_orig[cond_mprage | cond_spgr]

    if filtered_scan.shape[0] < 1:
        # TODO - LOG THIS
        cprint('NO MPRAGE Meta2: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
        return None

    scan = select_scan_qc_adni2(filtered_scan, mayo_mri_qc_subj, preferred_field_strength)

    sequence = replace_sequence_chars(scan.Sequence)

    return {'Subject_ID': subject_id,
            'VISCODE': timepoint,
            'Visit': visit_str,
            'Sequence': sequence,
            'Scan_Date': scan.ScanDate,
            'Study_ID': str(scan.StudyID),
            'Series_ID': str(scan.SeriesID),
            'Image_ID': str(scan.ImageUID),
            'Field_Strength': scan.MagStrength,
            'Original': True}


def adnigo_image(subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj, mri_quality_subj, mayo_mri_qc_subj, original_phase):

    if original_phase == 'ADNI1':
        filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                           & (mprage_meta_subj.MagStrength == 1.5)
                                           & (mprage_meta_subj.Visit == visit_str)
                                           & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]
        if filtered_mprage.shape[0] > 0:
            return adni1_image(subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj, mri_quality_subj, mayo_mri_qc_subj)

    mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
    return adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, mayo_mri_qc_subj)


def visits_to_timepoints_t1(subject, mprage_meta_subj_orig, adnimerge_subj):
    from datetime import datetime
    from clinica.iotools.converters.adni_to_bids.adni_utils import days_between
    from clinica.utils.stream import cprint

    mprage_meta_subj_orig = mprage_meta_subj_orig[mprage_meta_subj_orig['Visit'] != 'ADNI Baseline']

    visits = dict()

    unique_visits = list(mprage_meta_subj_orig.Visit.unique())

    pending_timepoints = []

    # We try to obtain the corresponding image Visit for a given VISCODE
    for adni_row in adnimerge_subj.iterrows():  # (adnimerge_subj[adnimerge_subj.FLDSTRENG.map(lambda x: x is not '')]).iterrows():
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
                cprint('[T1] Subject ' + subject + ' has multiple visits for one timepoint. ')
                # cprint(key_preferred_visit)
                # cprint(visits[key_preferred_visit])
                # cprint(visit)
            unique_visits.remove(preferred_visit_name)
            continue

        pending_timepoints.append(visit)

    # Then for images.Visit non matching the expected labels we find the closest date in visits list
    for visit in unique_visits:
        image = (mprage_meta_subj_orig[mprage_meta_subj_orig.Visit == visit]).iloc[0]
        min_db = 100000
        min_db2 = 0
        min_visit = None
        min_visit2 = None

        for timepoint in pending_timepoints:
            db = days_between(image.ScanDate, timepoint.EXAMDATE)
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
            cprint('[T1] Subject ' + subject + ' has multiple visits for one timepoint.')
            # cprint(key_min_visit)
            # cprint(visits[key_min_visit])
            # cprint(image.Visit)
    return visits


def select_scan_no_qc(scans_meta):

    selected_scan = scans_meta[scans_meta.Sequence.map(
        lambda x: (x.lower().find('repeat') < 0))]
    if selected_scan.shape[0] < 1:
        selected_scan = scans_meta

    scan = selected_scan.iloc[0]
    return scan


def select_scan_qc_adni2(scans_meta, mayo_mri_qc_subj, preferred_field_strength):
    import numpy as np

    multiple_mag_strength = False
    # Select preferred_field_strength images
    if len(scans_meta.MagStrength.unique()) > 1:
        multiple_mag_strength = True
        not_preferred_scan = scans_meta[scans_meta.MagStrength != preferred_field_strength]
        scans_meta = scans_meta[scans_meta.MagStrength == preferred_field_strength]

    if scans_meta.MagStrength.unique()[0] == 3.0:

        id_list = scans_meta.ImageUID.unique()
        image_ids = ['I' + str(imageuid) for imageuid in id_list]
        int_ids = [int(imageuid) for imageuid in id_list]
        images_qc = mayo_mri_qc_subj[mayo_mri_qc_subj.loni_image.isin(image_ids)]

        if images_qc.shape[0] > 0:
            selected_image = None
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
                    # If none of images passed qc the scan is None, otherwise:
                    if len(no_qc_ids) > 0:
                        no_qc_scans_meta = scans_meta[scans_meta.ImageUID.isin(no_qc_ids)]
                        return select_scan_no_qc(no_qc_scans_meta)
                else:
                    series_quality = [q if q > 0 else 4 for q in list(images_not_rejected.series_quality)]
                    best_q = np.amin(series_quality)
                    if best_q == 4:
                        best_q = -1
                    images_best_qc = images_not_rejected[images_not_rejected.series_quality == best_q]
                    if images_best_qc.shape[0] == 1:
                        selected_image = images_best_qc.iloc[0].loni_image[1:]
                    else:
                        best_ids = [int(x[1:]) for x in images_best_qc.loni_image.unique()]
                        best_qc_meta = scans_meta[scans_meta.ImageUID.isin(best_ids)]
                        return select_scan_no_qc(best_qc_meta)

            if selected_image is None and multiple_mag_strength:
                scans_meta = not_preferred_scan
            else:
                scan = scans_meta[scans_meta.ImageUID == int(selected_image)].iloc[0]
                return scan

    # 1.5T or 3.0T without QC
    return select_scan_no_qc(scans_meta)


def adni1_select_scanner(subj, csv_dir, adnimerge, timepoint):
    import pandas as pd
    import os

    com_philips = ['phil', 'n/a', 'not', 'phl', 'ND', 'na', 'n"a']
    sub_discarded = ['refused', 'declined', 'not completed', 'error', 'not on 3T scanner', '']

    mrimeta_tsv = pd.io.parsers.read_csv(os.path.join(csv_dir, 'MRIMETA.csv'), sep=',')
    mri3meta_tsv = pd.io.parsers.read_csv(os.path.join(csv_dir, 'MRI3META.csv'), sep=',')
    rid = list(adnimerge[adnimerge.PTID == subj].RID.unique())[0]

    # 1.5 T
    mrimeta_rid = mrimeta_tsv[mrimeta_tsv.RID == rid]
    mrimeta = mrimeta_rid[mrimeta_rid.VISCODE2 == timepoint]
    list_scanner_not_found = []

    if len(mrimeta) == 2:
        ix = mrimeta.index
        mrimeta = mrimeta.loc[ix[0]]

    if mrimeta.empty and timepoint == 'bl':
        timepoint = 'sc'
        mrimeta = mrimeta_rid[mrimeta_rid.VISCODE2 == timepoint]
        if len(mrimeta) > 1:

            mrimeta.drop(mrimeta.index[0], inplace=True)
        #    print mrimeta

    if not mrimeta.empty:
        scan = mrimeta.MMB1HEAD.item()
        com = mrimeta.MMB1BDCOM

        if isinstance(com, str):
            pass
        elif isinstance(com, float):
            pass
        else:
            com = com.item()

        if scan == '-4' or not scan:
            scanner = 'Philips'
        else:
            if not isinstance(com, float) and com.lower() in com_philips:
                scanner = 'Philips'
            else:
                scanner = 'GE_Siemens'
    if mrimeta.empty:
        # 3T
        mri3meta_rid = mri3meta_tsv[mri3meta_tsv.RID == rid]
        mri3meta = mri3meta_rid[mri3meta_rid.VISCODE2 == timepoint]

        if mri3meta.empty and timepoint == 'bl':
            timepoint = 'sc'
            mri3meta = mri3meta_rid[mri3meta_rid.VISCODE2 == timepoint]
        if not mri3meta.empty:
            scan = mri3meta.MMB1HEAD.item()
            com = mri3meta.MMB1BDCOM.item()
            if scan == '-4' or not scan:
                scanner = 'Philips'
            else:
                if com.lower() in com_philips:
                    scanner = 'Philips'
                else:
                    scanner = 'GE_Siemens'
        if mri3meta.empty:
            scanner = 'PROBLEM_SCANNER_NOT_FOUND'

    return scanner


def adni1_image_refactoring(csv_dir, adnimerge, subject_id, timepoint, visit_str, mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint
    # Get the preferred scan (image series that has been Scaled)
    filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                       & (mprage_meta_subj.Visit == visit_str)
                                       & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]

    # If there is not a preferred image we use ADNI2 processing
    # (get the best qc if available, otherwise the original) preferring 1.5T images
    if filtered_mprage.shape[0] < 1:
        mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
        return adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, mayo_mri_qc_subj, preferred_field_strength=1.5)

    filtered_mprage_mag = filtered_mprage
    if len(filtered_mprage.MagStrength.unique()) > 1:
        filtered_mprage_mag = filtered_mprage[filtered_mprage.MagStrength == 1.5]  # Select 1.5T images

    scan = filtered_mprage_mag.iloc[0]
    series_id = scan.SeriesID

    qc_passed = True
    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # print 'QC found but NOT passed'
        # print 'Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID)
        mprage_meta_subj_alt = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                                & (mprage_meta_subj.Visit == visit_str)
                                                & (mprage_meta_subj.SeriesID != series_id)]

        qc_prev_sequence = scan.Sequence
        scan = mprage_meta_subj_alt.iloc[0]
        series_id = scan.SeriesID
        qc_passed = False

    # TODO replace by the function to find if the scanner is_philips

    scanner = adni1_select_scanner(subject_id, csv_dir, adnimerge, timepoint)

    original = True

    # TODO replace the condition with if is_philips
    if scanner == 'Philips':

        scan = (mprage_meta_subj[
                    (mprage_meta_subj['Orig/Proc'] == 'Original') & (mprage_meta_subj.SeriesID == series_id)]).iloc[
            0]
        sequence = scan.Sequence

    else:  # scan already selected above
        sequence = scan.Sequence[:scan.Sequence.find('N3') - 2]
        original = False

    if not qc_passed:
        if scan.Sequence == 'MP-RAGE':
            original_img_seq = 'MPR'
        else:  # 'MP-RAGE REPEAT'
            original_img_seq = 'MPR-R'

        processing_seq = qc_prev_sequence[qc_prev_sequence.find(';'):qc_prev_sequence.find('N3') - 2]
        sequence = original_img_seq + processing_seq
        # print sequence

    sequence = replace_sequence_chars(sequence)

    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # TODO - LOG THIS
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))

    return {'Subject_ID': subject_id,
            'VISCODE': timepoint,
            'Visit': visit_str,
            'Sequence': sequence,
            'Scan_Date': scan.ScanDate,
            'Study_ID': str(scan.StudyID),
            'Series_ID': str(scan.SeriesID),
            'Image_ID': str(scan.ImageUID),
            'Field_Strength': scan.MagStrength,
            'Original': original}


def adnigo_image_refactoring(csv_dir, adnimerge, subject_id, timepoint, visit_str, mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj, original_phase):

    if original_phase == 'ADNI1':
        filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                           & (mprage_meta_subj.MagStrength == 1.5)
                                           & (mprage_meta_subj.Visit == visit_str)
                                           & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]
        if filtered_mprage.shape[0] > 0:
            return adni1_image_refactoring(csv_dir, adnimerge, subject_id, timepoint, visit_str, mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj)

    mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
    return adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, mayo_mri_qc_subj)
