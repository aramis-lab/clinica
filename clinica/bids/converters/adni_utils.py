# Paths computation methods
def compute_t1_paths(source_dir, clinical_dir, dest_dir, subjs_list):
    """

    :param source_dir:
    :param clinical_dir:
    :param dest_dir:
    :param subjs_list:
    :return:
    """

    import pandas as pd
    from os import path, walk

    t1_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                 'Study_ID', 'Field_Strength', 'Series_ID', 'Original']

    t1_df = pd.DataFrame(columns=t1_col_df)
    adni_merge_path = path.join(clinical_dir, 'ADNIMERGE.csv')
    adni_screening_path = path.join(clinical_dir, 'ADNI_ScreeningList_8_22_12.csv')
    ida_meta_path = path.join(clinical_dir, 'IDA_MR_METADATA_Listing.csv')
    mprage_meta_path = path.join(clinical_dir, 'MPRAGEMETA.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
    ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
    mprage_meta = pd.io.parsers.read_csv(mprage_meta_path, sep=',')

    for subj in subjs_list:
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        # Sort the values by examination date
        adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

        mprage_meta_subj = mprage_meta[mprage_meta.SubjectID == subj]
        mprage_meta_subj = mprage_meta_subj.sort_values('ScanDate')

        ida_meta_subj = ida_meta[ida_meta.Subject == subj]

        mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
        visits = self.visits_to_timepoints_t1(subj, mprage_meta_subj, adnimerge_subj)

        keys = visits.keys()
        keys.sort()
        for visit_info in visits.keys():
            if visit_info[1] == 'ADNI1':
                image_dict = self.adni1_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj,
                                              ida_meta_subj)
            elif visit_info[1] == 'ADNIGO':
                image_dict = self.adnigo_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj,
                                               ida_meta_subj,
                                               visit_info[2])
            else:  # ADNI2
                image_dict = self.adni2_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj_orig)

            if image_dict is None:
                image_dict = {'Subject_ID': subj,
                              'VISCODE': visit_info[0],
                              'Visit': visits[visit_info],
                              'Sequence': '',
                              'Scan_Date': '',
                              'Study_ID': '',
                              'Series_ID': '',
                              'Field_Strength': '',
                              'Original': True}

            row_to_append = pd.DataFrame(image_dict, index=['i', ])
            t1_df = t1_df.append(row_to_append, ignore_index=True)

    images = t1_df
    is_dicom = []
    nifti_paths = []
    count = 0

    for row in images.iterrows():

        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)

        count += 1
        # print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
        #     count) + ' / ' + str(total)

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

    # Drop all the lines that have the Path section empty
    images = images.drop(images[images.Path == ''].index)

    # Store the paths inside a file called conversion_info inside the input directory
    images.to_csv(path.join(dest_dir, 'conversion_info', 't1_paths.tsv'), sep='\t', index=False)

    return images


def compute_fdg_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list):
    import pandas as pd
    import os
    from os import walk, path
    from numpy import nan, argsort

    pet_fdg_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']

    pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
    petqc_path = path.join(clinical_dir, 'PETQC.csv')
    pet_meta_list_path = path.join(clinical_dir, 'PET_META_LIST.csv')
    petqc = pd.io.parsers.read_csv(petqc_path, sep=',')
    pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')

    for subj in subjs_list:
        pet_qc_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subj[-4:]))]
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]
        if subject_pet_meta.shape[0] < 1:
            # print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
            continue
        for visit in list(pet_qc_subj.VISCODE2.unique()):
            pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE2 == visit]
            if pet_qc_visit.shape[0] > 1:
                normal_images = []
                normal_meta = []
                for row in pet_qc_visit.iterrows():
                    image = row[1]
                    pet_meta_image = \
                        subject_pet_meta[(subject_pet_meta['Image ID'] == int(image.LONIUID[1:]))].iloc[0]
                    if pet_meta_image.Sequence.lower().find('early') < 0:
                        normal_images.append(image)
                        normal_meta.append(pet_meta_image)
                if len(normal_images) == 0:
                    print 'No regular FDG-PET image: Subject - ' + subj + ' for visit ' + visit
                    continue
                if len(normal_images) == 1:
                    qc_visit = normal_images[0]
                else:
                    qc_visit = None
                    index = argsort([x['Series ID'] for x in normal_meta])
                    for i in index[::-1]:
                        coreg_avg = subject_pet_meta[(subject_pet_meta['Sequence'] == 'Co-registered, Averaged')
                                                     & (
                                                         subject_pet_meta['Series ID'] == normal_meta[i][
                                                             'Series ID'])]
                        if coreg_avg.shape[0] > 0:
                            qc_visit = normal_images[i]
                            break
                    if qc_visit is None:
                        qc_visit = normal_images[index[len(index) - 1]]
            else:
                qc_visit = pet_qc_visit.iloc[0]
            int_image_id = int(qc_visit.LONIUID[1:])
            original_pet_meta = subject_pet_meta[
                (subject_pet_meta['Orig/Proc'] == 'Original') & (subject_pet_meta['Image ID'] == int_image_id)]
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & (subject_pet_meta.Sequence.map(
                    lambda x: (x.lower().find('fdg') > -1)))
                                                     & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]
                if original_pet_meta.shape[0] < 1:
                    print 'NO Screening: Subject - ' + subj + ' for visit ' + qc_visit.VISCODE2
                    continue
            original_image = original_pet_meta.iloc[0]
            averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == 'Co-registered, Averaged') & (
                subject_pet_meta['Series ID'] == original_image['Series ID'])]
            if averaged_pet_meta.shape[0] < 1:
                sel_image = original_image
                original = True
            else:
                sel_image = averaged_pet_meta.iloc[0]
                original = False
            visit = sel_image.Visit
            sequence = sel_image.Sequence
            sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                '_').replace(
                ')', '_').replace(':', '_')
            date = sel_image['Scan Date']
            study_id = sel_image['Study ID']
            series_id = sel_image['Series ID']
            image_id = sel_image['Image ID']

            row_to_append = pd.DataFrame(
                [[subj, qc_visit.VISCODE2, str(visit), sequence, date, str(study_id), str(series_id), str(image_id),
                  original]],
                columns=pet_fdg_col)
            pet_fdg_df = pet_fdg_df.append(row_to_append, ignore_index=True)

    images = pet_fdg_df
    count = 0
    total = images.shape[0]
    image_folders = []
    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)
        count += 1
        print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
            count) + ' / ' + str(total)
        image_path = ''
        for (dirpath, dirnames, filenames) in walk(seq_path):
            found = False
            for d in dirnames:
                if d == 'I' + str(image.Image_ID):
                    image_path = path.join(dirpath, d)
                    found = True
                    break
            if found:
                break
        if image.Original:
            for (dirpath, dirnames, filenames) in walk(image_path):
                for f in filenames:
                    if f.endswith(".nii"):
                        image_path = path.join(dirpath, f)
                        break
        image_folders.append(image_path)
        if image_path == '':
            print 'Not found ' + str(image.Subject_ID)
    images.loc[:, 'Path'] = pd.Series(image_folders, index=images.index)

    fdg_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(fdg_csv_path):
        os.mkdir(fdg_csv_path)
    images.to_csv(path.join(fdg_csv_path, 'fdg_pet_paths.tsv'), sep='\t', index=False)

    return images


def compute_av45_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list):
    import pandas as pd
    from numpy import nan
    from os import walk, path
    import os

    pet_av45_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                    'Study_ID', 'Series_ID', 'Image_ID', 'Original']

    pet_av45_df = pd.DataFrame(columns=pet_av45_col)

    pet_visits_corr = {'ADNI Baseline': 'bl', 'ADNI2 Baseline-New Pt': 'bl'}
    av45qc_path = path.join(clinical_dir, 'AV45QC.csv')
    av45qc = pd.io.parsers.read_csv(av45qc_path, sep=',')
    baseline = av45qc[(av45qc.VISCODE2 == 'bl') & (av45qc.PASS == 1) & av45qc.RID.isin(
        [int(s[-4:]) for s in subjs_list])]

    pet_meta_list_path = path.join(clinical_dir, 'PET_META_LIST.csv')
    pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')
    pet_meta_list = pet_meta_list[pet_meta_list.Visit.map(lambda x: x.find('Baseline') > -1)]

    for subject in subjs_list:
        pet_qc_subj = av45qc[(av45qc.PASS == 1) & (av45qc.RID == int(subject[-4:]))]
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subject]

        if subject_pet_meta.shape[0] < 1:
            # print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
            continue

        for visit in list(pet_qc_subj.VISCODE2.unique()):
            pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE2 == visit]
            if pet_qc_visit.shape[0] > 1:
                normal_images = []
                normal_meta = []
                for row in pet_qc_visit.iterrows():
                    image = row[1]
                    pet_meta_image = \
                        subject_pet_meta[(subject_pet_meta['Image ID'] == int(image.LONIUID[1:]))].iloc[0]
                    if pet_meta_image.Sequence.lower().find('early') < 0:
                        normal_images.append(image)
                        normal_meta.append(pet_meta_image)

                if len(normal_images) == 0:
                    print 'No regular AV45-PET image: Subject - ' + subject + ' for visit ' + visit
                    continue
                if len(normal_images) == 1:
                    qc_visit = normal_images[0]
                else:
                    qc_visit = None
                    index = nan.argsort([x['Series ID'] for x in normal_meta])
                    for i in index[::-1]:
                        coreg_avg = subject_pet_meta[
                            (subject_pet_meta['Sequence'] == 'AV45 Co-registered, Averaged')
                            & (subject_pet_meta['Series ID'] == normal_meta[i]['Series ID'])]
                        if coreg_avg.shape[0] > 0:
                            qc_visit = normal_images[i]
                            break
                    if qc_visit is None:
                        qc_visit = normal_images[index[len(index) - 1]]
            else:
                qc_visit = pet_qc_visit.iloc[0]
            int_image_id = int(qc_visit.LONIUID[1:])
            original_pet_meta = subject_pet_meta[
                (subject_pet_meta['Orig/Proc'] == 'Original') & (subject_pet_meta['Image ID'] == int_image_id)]
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & (subject_pet_meta.Sequence.map(
                    lambda x: (x.lower().find('av45') > -1)))
                                                     & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]
                if original_pet_meta.shape[0] < 1:
                    print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
                    continue
            original_image = original_pet_meta.iloc[0]
            averaged_pet_meta = subject_pet_meta[
                (subject_pet_meta['Sequence'] == 'AV45 Co-registered, Averaged') & (
                    subject_pet_meta['Series ID'] == original_image['Series ID'])]
            if averaged_pet_meta.shape[0] < 1:
                sel_image = original_image
                original = True
            else:
                sel_image = averaged_pet_meta.iloc[0]
                original = False
            visit = sel_image.Visit
            sequence = sel_image.Sequence
            sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                '_').replace(
                ')', '_').replace(':', '_')
            viscode = qc_visit['VISCODE2']
            date = sel_image['Scan Date']
            study_id = sel_image['Study ID']
            series_id = sel_image['Series ID']
            image_id = sel_image['Image ID']

            row_to_append = pd.DataFrame(
                [[subject, viscode, str(visit), sequence, date, str(study_id), str(series_id), str(image_id),
                  original]],
                columns=pet_av45_col)
            pet_av45_df = pet_av45_df.append(row_to_append, ignore_index=True)

    subjects = pet_av45_df
    count = 0
    total = subjects.shape[0]

    image_folders = []
    for row in subjects.iterrows():

        subject = row[1]
        seq_path = path.join(source_dir, str(subject.Subject_ID), subject.Sequence)

        count += 1
        print 'Processing subject ' + str(subject.Subject_ID) + ', ' + str(count) + ' / ' + str(total)

        image_path = ''
        for (dirpath, dirnames, filenames) in walk(seq_path):
            found = False
            for d in dirnames:
                if d == 'I' + str(subject.Image_ID):
                    image_path = path.join(dirpath, d)
                    found = True
                    break
            if found:
                break

        if subject.Original:
            for (dirpath, dirnames, filenames) in walk(image_path):
                for f in filenames:
                    if f.endswith(".nii"):
                        image_path = path.join(dirpath, f)
                        break

        image_folders.append(image_path)

        if image_path == '':
            print 'Not found ' + str(subject.Subject_ID)

    subjects.loc[:, 'Path'] = pd.Series(image_folders, index=subjects.index)

    av45_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(av45_csv_path):
        os.mkdir(av45_csv_path)
    subjects.to_csv(path.join(av45_csv_path, 'av45_pet_paths.tsv'), sep='\t', index=False)

    return subjects


def compute_fmri_path( source_dir, clinical_dir, dest_dir, subjs_list):
    '''
    Compute the paths to fmri images.

    The fmri images to convert into BIDS are chosen in the following way:
        - Extract the list of subjects from MAYOADIRL_MRI_FMRI_09_15_16.csv
        - Select the only the scans that came from PHILIPS machine (field Scanner from IDA_MR_Metadata_Listing.csv)
        - Discard all the subjects with column  series_quality = 4  (4 means that the scan is not usable) in MAYOADIRL_MRI_IMAGEQC_12_08_15.csv

    In case of multiple scans for the same session, same date the one to convert is chosen with the following criteria:
        - Check if in the file MAYOADIRL_MRI_IMAGEQC_12_08_15.csv there is a single scan with the field series_selected == 1
        - If yes choose the one with series_selected == 1
        - If no choose the scan with the best quality


    :param source_dir: path to the ADNI image folder
    :param clinical_dir: path to the directory with all the clinical data od ADNI
    :param dest_dir: path to the output_folder
    :param subjs_list: subjects list
    :return: pandas Dataframe containing the path for each fmri
    '''
    import os
    from os import path
    from os import walk
    import pandas as pd
    import logging

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
        print subj
        fmri_subjs_info = mayo_mri_fmri[(mayo_mri_fmri.RID == int(subj[-4:]))]
        # Extract visits available
        visits_list = fmri_subjs_info['VISCODE2'].tolist()
        # Removing duplicates
        visits_list = list(set(visits_list))

        if len(visits_list) != 0:
            for viscode in visits_list:
                scan_date = ''
                sequence = ''
                loni_uid = ''
                visit = ''
                mag_strenght = ''
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
                            imageuid_to_select = self.select_image_qc(fmri_imageuid, images_qc)

                        fmri_subj = fmri_subj[fmri_subj['IMAGEUID'] == int(imageuid_to_select)].iloc[0]
                    else:
                        fmri_subj = fmri_subj.iloc[0]

                    fmri_imageuid = fmri_subj['IMAGEUID']

                    # Discard scans made with non Philips scanner and with a bad quality
                    fmri_metadata = ida_mr_metadata[ida_mr_metadata['IMAGEUID'] == fmri_imageuid]

                    if not fmri_metadata.empty:
                        fmri_metadata = fmri_metadata.iloc[0]

                        if not 'Philips' in fmri_metadata['Scanner']:
                            print 'No Philips scanner for ', subj, 'visit', viscode, '. Skipped.'
                            continue

                        elif 4 in mayo_mri_imageqc[mayo_mri_imageqc['loni_image'] == 'I' + str(fmri_imageuid)][
                            'series_quality'].values:
                            print 'Bad scan quality for ', subj, 'visit', viscode, '. Skipped.'
                            continue

                        scan_date = fmri_subj.SCANDATE
                        sequence = self.replace_sequence_chars(fmri_subj.SERDESC)
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
                                        print 'Path not existing for subject:', subj, 'visit:', visit
                                    found = True
                                    break
                            if found:
                                break

                        # The session scmri correspond to the baseline
                        if viscode == 'scmri':
                            viscode = 'bl'
                    else:
                        print 'Missing visit, sequence, scan date and loniuid for ', subj, 'visit', visit
                        continue

                    row_to_append = pd.DataFrame(
                        [[subj, str(viscode), visit, str(fmri_imageuid), sequence, scan_date, str(loni_uid),
                          scanner, mag_strenght, image_path]], columns=fmri_col)

                    fmri_df = fmri_df.append(row_to_append, ignore_index=True)
                else:
                    logging.info('Missing fMRI for ', subj, 'visit', visit)

    fmri_df.to_csv(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t', index=False)
    return fmri_df


def compute_dti_paths(adni_dir, csv_dir, dest_dir, subjs_list):
    import pandas as pd
    from os import path, walk

    dti_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                  'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner', 'Enhanced']

    dti_df = pd.DataFrame(columns=dti_col_df)

    adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
    ida_meta_path = path.join(csv_dir, 'IDA_MR_METADATA_Listing.csv')
    mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

    adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
    ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
    ida_meta = ida_meta[ida_meta.Sequence.map(lambda x: x.lower().find('dti') > -1)]

    mri_qc = pd.io.parsers.read_csv(mri_qc_path, sep=',')
    mri_qc = mri_qc[mri_qc.series_type == 'DTI']

    for subj in subjs_list:

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
    images.to_csv(path.join(dest_dir, 'conversion_info', 'dwi_paths.tsv'), sep='\t', index=False)
    return images


def visits_to_timepoints_dti(subject, ida_meta_subj, adnimerge_subj):
    from datetime import datetime

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
                dif = self.days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
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


def dti_image(subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):
    sel_image = select_image_qc(list(ida_meta_scans.IMAGEUID), mri_qc_subj)
    if sel_image is None:
        return None

    sel_scan = ida_meta_scans[ida_meta_scans.IMAGEUID == sel_image].iloc[0]

    sequence = sel_scan.Sequence
    sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                        '_').replace(
        ')', '_').replace(':', '_')

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

        images_best_qc = images_not_rejected[images_not_rejected.series_quality == best_q]
        if images_best_qc.shape[0] == 1:
            selected_image = images_best_qc.iloc[0].loni_image[1:]
        else:
            selected_image = min(int_ids)

    return int(selected_image)


def days_between(d1, d2):
    from datetime import datetime
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


# Images correction methods
def center_nifti_origin(input_image, output_image):
    import nibabel as nib
    import numpy as np

    img = nib.load(input_image)
    canonical_img = nib.as_closest_canonical(img)
    hd = canonical_img.header
    if hd['quatern_b'] != 0 or hd['quatern_c'] != 0 or hd['quatern_d'] != 0:
        print 'Warning: Not all values in quatern are equal to zero'
    qform = np.zeros((4, 4))
    for i in range(1, 4):
        qform[i - 1, i - 1] = hd['pixdim'][i]
        qform[i - 1, 3] = -1.0 * hd['pixdim'][i] * hd['dim'][i] / 2.0
    new_img = nib.Nifti1Image(canonical_img.get_data(caching='unchanged'), qform)
    nib.save(new_img, output_image)


# Auxiliary methods
def check_two_dcm_folder(dicom_path, bids_folder, image_uid):
    '''
    Check if a folder contains more than one DICOM and if yes, copy the DICOM related to image id passed as parameter into
    a temporary folder called tmp_dicom_folder.


    :param dicom_path: path to the DICOM folder
    :param bids_folder: path to the BIDS folder where the dataset will be stored
    :param image_uid: image id of the fmri
    :return: the path to the original DICOM folder or the path to a temporary folder called tmp_dicom_folder where only
     the DICOM to convert is copied
    '''
    from glob import glob
    from os import path
    from shutil import copy
    import shutil
    import os

    temp_folder_name = 'tmp_dcm_folder'
    dest_path = path.join(bids_folder, temp_folder_name)

    # Check if there is more than one xml file inside the folder
    xml_list = glob(path.join(dicom_path,'*.xml*'))

    if len(xml_list) > 1:
        # Remove the precedent tmp_dcm_folder if is existing
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        os.mkdir(dest_path)
        dmc_to_conv = glob(path.join(dicom_path,'*'+image_uid+'.dcm*'))
        for d in dmc_to_conv:
            copy(d, dest_path)
        return dest_path
    else:
        return dicom_path




