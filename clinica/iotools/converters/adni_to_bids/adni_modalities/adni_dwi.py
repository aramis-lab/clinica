__author__ = "Jorge Samper and Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = ""
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper"
__email__ = ""
__status__ = "Development"

def compute_dti_paths(adni_dir, csv_dir, dest_dir, subjs_list):
    '''
    Computes paths to DWI images of ADNI

    The pandas datataframe returned as output cointains the following information:
        - Subject_ID = id of the subjects
        - VISCODE = id of the session
        - Visit =
        - Sequence
        - Scan_Date = acquisition date
        - Study_ID
        - Series_ID
        - Image_ID = id of the image
        - Scanner = name of the scanner used for the acquisition
        - Enhanced

    :param adni_dir: path to the input adni directory
    :param csv_dir: path to the directory containing the clinical data
    :param dest_dir: path to the destination list
    :param subjs_list: a list containing the subjects to process
    :return: a pandas dataframe
    '''

    import pandas as pd
    from os import path, walk

    dti_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                  'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner', 'Enhanced']

    dti_df = pd.DataFrame(columns=dti_col_df)

    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
    ida_meta_path = path.join(csv_dir, 'IDA_MR_Metadata_Listing.csv')
    mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
    ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
    ida_meta = ida_meta[ida_meta.Sequence.map(lambda x: x.lower().find('dti') > -1)]

    mri_qc = pd.io.parsers.read_csv(mri_qc_path, sep=',')
    mri_qc = mri_qc[mri_qc.series_type == 'DTI']
    print '=================================================='

    for subj in subjs_list:
        print 'Computing path for subj', subj

        adnimerge_subj = adni_merge[adni_merge.PTID == subj]

        # Sort the values by examination date
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')
        ida_meta_subj = ida_meta[ida_meta.Subject == subj]
        ida_meta_subj = ida_meta_subj.sort_values('Scan Date')
        mri_qc_subj = mri_qc[mri_qc.RID == int(subj[-4:])]

        visits = visits_to_timepoints_dti(subj, ida_meta_subj, adnimerge_subj)

        keys = visits.keys()
        keys.sort()

        for visit_info in visits.keys():

            # visit_info = (VISCODE, COLPROT, ORIGPROT)

            visit_str = visits[visit_info]
            visit_ida_meta = ida_meta_subj[ida_meta_subj.Visit == visit_str]

            axial_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') < 0)]
            enhanced_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') > -1)]

            axial = dti_image(subj, visit_info[0], visits[visit_info], axial_ida_meta, mri_qc_subj, False)
            enhanced = dti_image(subj, visit_info[0], visits[visit_info], enhanced_ida_meta, mri_qc_subj, True)

            if axial is not None:
                row_to_append = pd.DataFrame(axial, index=['i', ])
                dti_df = dti_df.append(row_to_append, ignore_index=True)

            if enhanced is not None:
                row_to_append = pd.DataFrame(enhanced, index=['i', ])
                dti_df = dti_df.append(row_to_append, ignore_index=True)

    images = dti_df
    is_dicom = []
    nifti_paths = []
    count = 0

    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(adni_dir, str(image.Subject_ID), image.Sequence)

        count += 1
        # print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
        #     count) + ' / ' + str(total)

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

        # for (dirpath, dirnames, filenames) in walk(series_path):
        #     for f in filenames:
        #         if f.endswith(".nii"):
        #             dicom = False
        #             nifti_path = path.join(dirpath, f)
        #             break

        is_dicom.append(dicom)
        nifti_paths.append(nifti_path)

    images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
    images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

    # Drop all the lines that have the Path section empty
    # images = images.drop(images[images.Path == ''].index)
    # Store the paths inside a file called t1_paths inside the input directory
    print '\nDone! Saving the results into',path.join(dest_dir, 'conversion_info', 'dwi_paths.tsv')
    images.to_csv(path.join(dest_dir, 'conversion_info', 'dwi_paths.tsv'), sep='\t', index=False)
    print '=================================================='

    return images


def convert_dwi(dest_dir, dwi_paths, mod_to_add=False, mod_to_update=False):
    '''


    :param dest_dir: path to the BIDS ADNI directory
    :param dwi_paths: pandas dataframe containing all the paths to the dwi images
    :param mod_to_add:
    :param mod_to_update:
    :return:
    '''

    import clinica.iotools.bids_utils as bids
    import clinica.iotools.converters.adni_utils as adni_utils
    from os import path
    import os
    from glob import glob
    import shutil

    subjs_list = dwi_paths['Subject_ID'].drop_duplicates().values

    for i in range(0, len(subjs_list)):
        print '--Converting dwi for subject', subjs_list[i], '--'
        alpha_id = bids.remove_space_and_symbols(subjs_list[i])
        bids_id = 'sub-ADNI' + alpha_id
        # Extract the list of sessions available from the dwi paths files, removing the duplicates
        sess_list = dwi_paths[(dwi_paths['Subject_ID'] == subjs_list[i])]['VISCODE'].drop_duplicates().values

        if not os.path.exists(path.join(dest_dir, bids_id)):
            os.mkdir(path.join(dest_dir, bids_id))

        # For each session available, create the folder if doesn't exist and convert the files
        for ses in sess_list:

            ses_bids = adni_utils.viscode_to_session(ses)
            bids_ses_id = 'ses-' + ses_bids
            bids_file_name = bids_id + '_ses-' + ses_bids
            ses_path = path.join(dest_dir, bids_id, bids_ses_id)

            existing_nii = glob(path.join(ses_path, 'func', '*nii*'))

            if mod_to_add:
                if len(existing_nii) > 0:
                    print 'DWI folder already existing. Skipped.'
                    continue

            if mod_to_update:
                if os.path.exists(path.join(ses_path, 'dwi')):
                    print 'Removing the old dwi file for session', ses
                    shutil.rmtree(path.join(ses_path, 'dwi'))
                else:
                    print 'Adding a new dwi file for session', ses

            if not os.path.exists(ses_path):
                os.mkdir(ses_path)

            dwi_info = dwi_paths[(dwi_paths['Subject_ID'] == subjs_list[i]) & (dwi_paths['VISCODE'] == ses)]

            # For the same subject, same session there could be multiple dwi with different acq label
            for j in range(0, len(dwi_info)):
                dwi_subj = dwi_info.iloc[j]
                if type(dwi_subj['Path']) != float and dwi_subj['Path'] != '':
                    if not os.path.exists(path.join(ses_path, 'dwi')):
                        os.mkdir(path.join(ses_path, 'dwi'))
                    dwi_path = dwi_subj['Path']
                    bids_name = bids_file_name + '_acq-' + (
                        'axialEnhanced' if dwi_subj['Enhanced'] else 'axial') + '_dwi'
                    # bids.dcm_to_nii(dwi_path, path.join(ses_path, 'dwi'), bids_name)

                    bids_dest_dir = path.join(ses_path, 'dwi')

                    if not os.path.exists(bids_dest_dir):
                        os.mkdir(dest_dir)
                    os.system('dcm2niix -b n -z y -o ' + bids_dest_dir + ' -f ' + bids_name + ' ' + dwi_path)

                    # If dcm2niix didn't work use dcm2nii
                    print path.join(dest_dir, bids_name + '.nii.gz')
                    if not os.path.exists(path.join(bids_dest_dir, bids_name + '.nii.gz')) or not os.path.exists(
                                    path.join(bids_dest_dir, bids_name + '.bvec') or not os.path.exists(
                                    path.join(bids_dest_dir, bids_name + '.bval'))):
                        print '\nConversion with dcm2niix failed, trying with dcm2nii'

                        # Find all the files eventually created by dcm2niix and remove them
                        dwi_dcm2niix = glob(path.join(bids_dest_dir, bids_name + '*'))

                        for d in dwi_dcm2niix:
                            print 'Removing the old', d
                            os.remove(d)

                        os.system(
                            'dcm2nii -a n -d n -e n -i y -g y -p n -m n -r n -x n -o ' + bids_dest_dir + ' ' + dwi_path)
                        nii_file = path.join(bids_dest_dir, subjs_list[i].replace('_', '') + '.nii.gz')
                        bvec_file = path.join(bids_dest_dir, subjs_list[i].replace('_', '') + '.bvec')
                        bval_file = path.join(bids_dest_dir, subjs_list[i].replace('_', '') + '.bval')

                        if os.path.exists(bvec_file) and os.path.exists(bval_file):
                            os.rename(bvec_file, path.join(bids_dest_dir, bids_name + '.bvec'))
                            os.rename(bval_file, path.join(bids_dest_dir, bids_name + '.bval'))
                        else:
                            print 'WARNING: bvec and bval not generated by dcm2nii'

                        if os.path.exists(nii_file):
                            os.rename(nii_file, path.join(bids_dest_dir, bids_name + '.nii.gz'))
                        else:
                            print 'WARNING: CONVERSION FAILED...'

        print '--Conversion finished--\n'


def dti_image(subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):

    from clinica.iotools.converters.adni_utils import replace_sequence_chars

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
            # TODO verify if new code above works. Previously we had this:
            # selected_image = min(int_ids)

    return int(selected_image)


def visits_to_timepoints_dti(subject, ida_meta_subj, adnimerge_subj):
    from datetime import datetime
    from clinica.iotools.converters.adni_utils import days_between

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
                print 'Multiple visits for one timepoint!'
                print subject
                print key_preferred_visit
                print visits[key_preferred_visit]
                print visit
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
            print 'No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit
            print image
            continue

        if min_visit2 is not None and min_db > 90:
            print 'More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit + ' on ' + image.ScanDate
            print 'Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                min_db) + ' days)'
            print 'Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                min_db2) + ' days)'

            # If image is too close to the date between two visits we prefer the earlier visit
            if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                    > datetime.strptime(image.ScanDate, "%Y-%m-%d")
                    > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                dif = days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                if abs((dif / 2.0) - min_db) < 30:
                    min_visit = min_visit2

            print 'We prefer ' + min_visit.VISCODE

        key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
        if key_min_visit not in visits.keys():
            visits[key_min_visit] = image.Visit
        elif visits[key_min_visit] != image.Visit:
            print 'Multiple visits for one timepoint!'
            print subject
            print key_min_visit
            print visits[key_min_visit]
            print image.Visit

    return visits
