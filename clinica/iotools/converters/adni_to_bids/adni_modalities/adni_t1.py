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
    from os import path

    from pandas.io import parsers
    from colorama import Fore

    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids.adni_utils import t1_pet_paths_to_bids

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = parsers.read_csv(adni_merge_path, sep=',')
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of T1 images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of T1 images found. Exporting images into BIDS ...')
    # t1_pet_paths_to_bids(images, dest_dir, 't1')
    # cprint(Fore.GREEN + 'T1 conversion done.' + Fore.RESET)


def compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """

    Compute paths to t1 images of ADNI.

    :param source_dir:
    :param csv_dir:
    :param dest_dir:
    :param subjs_list:
    :return: a pandas dataframe
    """

    from os import path, walk, makedirs

    import pandas as pd

    from clinica.utils.stream import cprint

    t1_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                 'Study_ID', 'Field_Strength', 'Series_ID', 'Image_ID', 'Original']

    t1_df = pd.DataFrame(columns=t1_col_df)

    # Getting paths for needed .csv files
    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
    mprage_meta_path = path.join(csv_dir, 'MPRAGEMETA.csv')
    mri_quality_path = path.join(csv_dir, 'MRIQUALITY.csv')
    mayo_mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    # Loading needed .csv files
    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',', low_memory=False)
    mprage_meta = pd.io.parsers.read_csv(mprage_meta_path, sep=',', low_memory=False)
    mri_quality = pd.io.parsers.read_csv(mri_quality_path, sep=',', low_memory=False)
    mayo_mri_qc = pd.io.parsers.read_csv(mayo_mri_qc_path, sep=',', low_memory=False)
    # Keep only T1 scans
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'T1']

    # We will convert the images for each subject in the subject list
    for subj in subjs_list:
        cprint(subj)

        # Filter ADNIMERGE, MPRAGE METADATA and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mprage_meta_subj = mprage_meta[mprage_meta.SubjectID == subj]
        mprage_meta_subj = mprage_meta_subj.sort_values('ScanDate')

        mri_quality_subj = mri_quality[mri_quality.RID == int(subj[-4:])]
        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints_t1(subj, mprage_meta_subj, adnimerge_subj)

        keys = list(visits.keys())
        keys.sort()

        for visit_info in keys:
            cohort = visit_info[1]
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            image_dict = None

            if cohort in ('ADNI1', 'ADNIGO', 'ADNI2'):
                image_dict = adni1GO2_image(subj, timepoint, visit_str, mprage_meta_subj, mri_quality_subj,
                                            mayo_mri_qc_subj,
                                            preferred_field_strength=1.5 if cohort == 'ADNI1' else 3.0)
            elif cohort == 'ADNI3':
                image_dict = adni3_image(subj, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj)
            else:
                cprint("Subject %s visit %s belongs to an unknown cohort: %s" % (subj, visit_str, cohort))

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

    images = t1_df
    is_dicom = []
    nifti_paths = []
    count = 0

    # For each image a path will be reconstructed and checked if it exists
    for row in images.iterrows():

        image = row[1]
        # Path to folder with the corresponding image sequence inside subject directory
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)

        # Find path to image folder with the corresponding image series inside image sequence directory
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

        # Check if there is a NIFTI image inside the folder. Otherwise assume it is a DICOM image
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

    # Columns with if image is DICOM and path information are added to the dataframe
    images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
    images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

    # Store the paths inside a file called conversion_info inside the input directory
    t1_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(t1_tsv_path):
        makedirs(t1_tsv_path)
    images.to_csv(path.join(t1_tsv_path, 't1_paths.tsv'), sep='\t', index=False)

    return images


def adni1GO2_image(subject_id, timepoint, visit_str, mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj,
                   preferred_field_strength=3.0):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars

    # Get the preferred scan (image series that has been Scaled)
    filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                       & (mprage_meta_subj.Visit == visit_str)
                                       & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]

    # If no preferred image found, get N3 processed image (N3m)
    if filtered_mprage.shape[0] < 1:
        filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                           & (mprage_meta_subj.Visit == visit_str)
                                           & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('N3m')))]

        # If no N3 processed image found (it means there are no processed images at all), get best original image
        if filtered_mprage.shape[0] < 1:
            return original_image(subject_id, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj,
                                  preferred_field_strength)

    # If there are images with different magnetic field strength, prefer 1.5T images for ADNI1, 3.0T otherwise
    if len(filtered_mprage.MagStrength.unique()) > 1:
        filtered_mprage = filtered_mprage[filtered_mprage.MagStrength == preferred_field_strength]

    # Sort by Series ID in case there are several images, so we keep the one acquired first
    filtered_mprage = filtered_mprage.sort_values('SeriesID')

    scan = filtered_mprage.iloc[0]

    # Check if selected scan passes QC (if QC exists)
    if not check_qc(scan, subject_id, visit_str, mprage_meta_subj, mri_quality_subj):
        return None

    n3 = scan.Sequence.find('N3')
    # Sequence ends in 'N3' or in 'N3m'
    sequence = scan.Sequence[:n3 + 2 + int(scan.Sequence[n3 + 2] == 'm')]
    sequence = replace_sequence_chars(sequence)

    return {'Subject_ID': subject_id,
            'VISCODE': timepoint,
            'Visit': visit_str,
            'Sequence': sequence,
            'Scan_Date': scan.ScanDate,
            'Study_ID': str(scan.StudyID),
            'Series_ID': str(scan.SeriesID),
            'Image_ID': str(scan.ImageUID),
            'Field_Strength': scan.MagStrength,
            'Original': False}


def original_image(subject_id, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj, preferred_field_strength=3.0):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint

    mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']

    cond_mprage = ((mprage_meta_subj_orig.Visit == visit_str)
                   & mprage_meta_subj_orig.Sequence.map(
                       lambda x: ((x.lower().find('mprage') > -1)
                                  | (x.lower().find('mp-rage') > -1)
                                  | (x.lower().find('mp rage') > -1))
                       & (x.find('2') < 0)
                       )
                   )

    cond_spgr = ((mprage_meta_subj_orig.Visit == visit_str)
                 & mprage_meta_subj_orig.Sequence.map(
                    lambda x: (x.lower().find('spgr') > -1)
                    & (x.lower().find('acc') < 0)
                    )
                 )

    filtered_scan = mprage_meta_subj_orig[cond_mprage | cond_spgr]

    if filtered_scan.shape[0] < 1:
        # TODO - LOG THIS
        cprint('NO MPRAGE Meta: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
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


def adni3_image(subject_id, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj):

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint

    filtered_scan = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                     & (mprage_meta_subj.Visit == visit_str)
                                     & mprage_meta_subj.Sequence.map(
        lambda x: (x.lower().find('accel') > -1)
                  & ~(x.lower().endswith('_ND')))]

    if filtered_scan.shape[0] < 1:
        # TODO - LOG THIS
        cprint('NO MPRAGE Meta for ADNI3: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
        return None

    scan = select_scan_qc_adni2(filtered_scan, mayo_mri_qc_subj, preferred_field_strength=3.0)
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


def check_qc(scan, subject_id, visit_str, mprage_meta_subj, mri_quality_subj):
    from clinica.utils.stream import cprint

    series_id = scan.SeriesID
    qc_passed = True

    # Check if QC exists for image series
    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]

    # If QC exists and failed we keep the other scan (in case 2 scans were performed)
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))

        mprage_meta_subj_alt = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                                & (mprage_meta_subj.Visit == visit_str)
                                                & (mprage_meta_subj.SeriesID != series_id)]

        qc_prev_sequence = scan.Sequence
        scan = mprage_meta_subj_alt.iloc[0]
        series_id = scan.SeriesID
        qc_passed = False

    if not qc_passed:
        if scan.Sequence == 'MP-RAGE':
            original_img_seq = 'MPR'
        else:  # 'MP-RAGE REPEAT'
            original_img_seq = 'MPR-R'

        processing_seq = qc_prev_sequence[qc_prev_sequence.find(';'):qc_prev_sequence.find('Scaled') - 2]
        sequence = original_img_seq + processing_seq
        # print sequence

    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # TODO - LOG THIS
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))
        return False

    return True


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

        if visit.ORIGPROT == 'ADNI3':
            if visit.VISCODE == 'bl':
                preferred_visit_name = 'ADNI Screening'
            else:
                year = str(int(visit.VISCODE[1:]) / 12)
                preferred_visit_name = 'ADNI3 Year ' + year + ' Visit'
        elif visit.ORIGPROT == 'ADNI2':
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
