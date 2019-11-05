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
    """Convert T1 MR images of ADNI into BIDS format

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    """

    from os import path

    from pandas.io import parsers
    from colorama import Fore

    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids.adni_utils import t1_pet_paths_to_bids

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = parsers.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of T1 images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of T1 images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 't1')
    cprint(Fore.GREEN + 'T1 conversion done.' + Fore.RESET)


def compute_t1_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """Compute the paths to T1 MR images and store them in a tsv file

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns:
        images: a dataframe with all the paths to the T1 MR images that will be converted into BIDS

    """

    from os import path, makedirs

    import pandas as pd

    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids.adni_utils import find_image_path

    t1_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                 'Study_ID', 'Field_Strength', 'Series_ID', 'Image_ID', 'Original']
    t1_df = pd.DataFrame(columns=t1_col_df)
    t1_dfs_list = []

    # Loading needed .csv files
    adni_merge = pd.read_csv(path.join(csv_dir, 'ADNIMERGE.csv'), sep=',', low_memory=False)
    mprage_meta = pd.read_csv(path.join(csv_dir, 'MPRAGEMETA.csv'), sep=',', low_memory=False)
    mri_quality = pd.read_csv(path.join(csv_dir, 'MRIQUALITY.csv'), sep=',', low_memory=False)
    mayo_mri_qc = pd.read_csv(path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv'), sep=',', low_memory=False)
    # Keep only T1 scans
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == 'T1']

    # We will convert the images for each subject in the subject list
    for subj in subjs_list:

        # Filter ADNIMERGE, MPRAGE METADATA and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mprage_meta_subj = mprage_meta[mprage_meta.SubjectID == subj]
        mprage_meta_subj = mprage_meta_subj.sort_values('ScanDate')

        mri_quality_subj = mri_quality[mri_quality.RID == int(subj[-4:])]
        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints_t1(subj, mprage_meta_subj, adnimerge_subj)

        for visit_info in visits.keys():
            cohort = visit_info[1]
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            image_dict = None

            if cohort in ('ADNI1', 'ADNIGO', 'ADNI2'):
                image_dict = adni1go2_image(subj, timepoint, visit_str, mprage_meta_subj, mri_quality_subj,
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
            t1_dfs_list.append(row_to_append)

    t1_df = pd.concat(t1_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
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
                         ('036_S_4899', 'm03'),
                         ('033_S_1016', 'm120'),
                         ('130_S_4984', 'm12'),
                         ('027_S_4802', 'm06'),
                         ('131_S_0409', 'bl'),
                         ('082_S_4224', 'm24'),
                         ('006_S_4960', 'bl'),
                         ('006_S_4960', 'm03'),
                         ('006_S_4960', 'm06'),
                         ('006_S_4960', 'm12'),
                         ('006_S_4960', 'm24'),
                         ('006_S_4960', 'm36'),
                         ('006_S_4960', 'm72'),
                         ('022_S_5004', 'bl'),
                         ('022_S_5004', 'm03'),
                         # T1wa
                         ('006_S_4485', 'm84')]

    # Removing known exceptions from images to convert
    error_indices = t1_df.index[t1_df.apply(lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1)]
    t1_df.drop(error_indices, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(t1_df, source_dir, 'T1', 'S', 'Series_ID')

    # Store the paths inside a file called conversion_info inside the input directory
    t1_tsv_path = path.join(dest_dir, 'conversion_info')
    if not path.exists(t1_tsv_path):
        makedirs(t1_tsv_path)
    images.to_csv(path.join(t1_tsv_path, 't1_paths.tsv'), sep='\t', index=False)

    return images


def visits_to_timepoints_t1(subject, mprage_meta_subj_orig, adnimerge_subj):
    """Finds the corresponding visit names and visit codes in months for a subject visits

    Args:
        subject: string containing subject ID
        mprage_meta_subj_orig: DataFrame of MPRAGE metadata of original images corresponding to the subject
        adnimerge_subj: DataFrame of ADNIMERGE informations corresponding to the subject

    Returns:
        Dictionary containing, for each visit, an entry with corresponding information

    """
    from datetime import datetime
    from clinica.iotools.converters.adni_to_bids.adni_utils import days_between
    from clinica.utils.stream import cprint

    mprage_meta_subj_orig = mprage_meta_subj_orig[mprage_meta_subj_orig['Visit'] != 'ADNI Baseline']

    visits = dict()
    unique_visits = list(mprage_meta_subj_orig.Visit.unique())
    pending_timepoints = []

    # We try to obtain the corresponding image Visit for a given VISCODE
    for adni_row in adnimerge_subj.iterrows():
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
            cprint('More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' +
                   image.Visit + ' on ' + image.ScanDate)
            cprint('Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE +
                   ' (Distance: ' + str(min_db) + ' days)')
            cprint('Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE +
                   ' (Distance: ' + str(min_db2) + ' days)')

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

    return visits


def adni1go2_image(subject_id, timepoint, visit_str, mprage_meta_subj, mri_quality_subj, mayo_mri_qc_subj,
                   preferred_field_strength=3.0):
    """Selects the preferred scan for a subject in a visit, given the subject belongs to ADNI 1, Go or 2 cohorts

    Args:

        subject_id: string containing subject ID
        timepoint: string of visit code in months
        visit_str: string of visit name
        mprage_meta_subj: DataFrame of MPRAGE metadata of images corresponding to the subject
        mri_quality_subj: DatFrame of MR image quality of images corresponding to the subject
        mayo_mri_qc_subj: DatFrame of MAYO Clinic MR image quality of images corresponding to the subject
        preferred_field_strength: Field strength that is preferred in case there are several image acquisitions

    Returns:
        Dictionary containing selected scan information

    """

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars

    # Get the preferred scan (image series that has been Scaled)
    filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                       & (mprage_meta_subj.Visit == visit_str)
                                       & mprage_meta_subj.Sequence.str.endswith('Scaled', na=False)]

    # If no preferred image found, get N3 processed image (N3m)
    if filtered_mprage.shape[0] < 1:
        filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                           & (mprage_meta_subj.Visit == visit_str)
                                           & mprage_meta_subj.Sequence.str.endswith('N3m', na=False)]

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


def adni3_image(subject_id, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj):
    """Selects the preferred scan for a subject in a visit, given the subject belongs to ADNI 3 cohort

    Args:
        subject_id: string containing subject ID
        timepoint: string of visit code in months
        visit_str: string of visit name
        mprage_meta_subj: DataFrame of MPRAGE metadata of images corresponding to the subject
        mayo_mri_qc_subj: DatFrame of MAYO Clinic MR image quality of images corresponding to the subject

    Returns:
        Dictionary containing selected scan information

    """

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint

    filtered_scan = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                     & (mprage_meta_subj.Visit == visit_str)
                                     & mprage_meta_subj.Sequence.str.contains('accel', case=False, na=False)
                                     & ~mprage_meta_subj.Sequence.str.lower().str.endswith('_nd', na=False)]

    if filtered_scan.shape[0] < 1:
        # TODO - LOG THIS
        cprint('NO MPRAGE Meta for ADNI3: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
        return None

    scan = select_scan_from_qc(filtered_scan, mayo_mri_qc_subj, preferred_field_strength=3.0)
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


def original_image(subject_id, timepoint, visit_str, mprage_meta_subj, mayo_mri_qc_subj, preferred_field_strength=3.0):
    """Selects the preferred scan for a subject in a visit, given the scan is not preprocessed

    Args:
        subject_id:string containing subject ID
        timepoint: string of visit code in months
        visit_str: string of visit name
        mprage_meta_subj: DataFrame of MPRAGE metadata of images corresponding to the subject
        mayo_mri_qc_subj: DatFrame of MAYO Clinic MR image quality of images corresponding to the subject
        preferred_field_strength: Field strength that is preferred in case there are several image acquisitions


    Returns:
        Dictionary containing selected scan information

    """

    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint

    mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']

    cond_mprage = ((mprage_meta_subj_orig.Visit == visit_str)
                   & mprage_meta_subj_orig.Sequence.str.contains("mprage|mp-rage|mp rage", case=False, na=False)
                   & ~mprage_meta_subj_orig.Sequence.str.contains("2", na=False))

    cond_spgr = ((mprage_meta_subj_orig.Visit == visit_str)
                 & mprage_meta_subj_orig.Sequence.str.contains("spgr", case=False, na=False)
                 & ~mprage_meta_subj_orig.Sequence.str.contains("acc", case=False, na=False))

    filtered_scan = mprage_meta_subj_orig[cond_mprage | cond_spgr]

    if filtered_scan.shape[0] < 1:
        # TODO - LOG THIS
        cprint('NO MPRAGE Meta: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str)
        return None

    scan = select_scan_from_qc(filtered_scan, mayo_mri_qc_subj, preferred_field_strength)
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


def select_scan_from_qc(scans_meta, mayo_mri_qc_subj, preferred_field_strength):
    """Selects a scan from a list of scans taking into account available QC

    Args:
        scans_meta: DataFrame containing the metadata for images to choose among
        mayo_mri_qc_subj: DatFrame of MAYO Clinic MR image quality of images corresponding to the subject
        preferred_field_strength: Field strength that is preferred in case there are several image acquisitions

    Returns:
        DataFrame row containing selected scan

    """
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


def select_scan_no_qc(scans_meta):
    """Selects a scan from available scans in case there is not an available QC

    Args:
        scans_meta: DataFrame containing the metadata for images to choose among

    Returns:
        DataFrame row containing selected scan

    """

    # We choose the first scan (not containing repeat in name)
    selected_scan = scans_meta[scans_meta.Sequence.str.contains("repeat", case=False, na=False)]

    if selected_scan.shape[0] < 1:
        selected_scan = scans_meta

    scan = selected_scan.iloc[0]
    return scan


def check_qc(scan, subject_id, visit_str, mprage_meta_subj, mri_quality_subj):
    """Checks if a scan has passed check quality, if it was performed

    Args:
        scan: DataFrame row containing scan information
        subject_id: string containing subject ID
        visit_str: string of visit name
        mprage_meta_subj: DataFrame of MPRAGE metadata of images corresponding to the subject
        mri_quality_subj: DatFrame of MR image quality of images corresponding to the subject

    Returns: boolean, True if image passed QC or if there is not available QC, False elsewhere

    """
    from clinica.utils.stream import cprint

    series_id = scan.SeriesID

    # Check if QC exists for image series
    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]

    # If QC exists and failed we keep the other scan (in case 2 scans were performed)
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))

        mprage_meta_subj_alt = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Original')
                                                & (mprage_meta_subj.Visit == visit_str)
                                                & (mprage_meta_subj.SeriesID != series_id)]

        scan = mprage_meta_subj_alt.iloc[0]

    qc = mri_quality_subj[mri_quality_subj.LONIUID == 'S' + str(scan.SeriesID)]
    if qc.shape[0] > 0 and qc.iloc[0].PASS != 1:
        # TODO - LOG THIS
        cprint('QC found but NOT passed')
        cprint('Subject ' + subject_id + ' - Series: ' + str(scan.SeriesID) + ' - Study: ' + str(scan.StudyID))
        return False

    return True
