
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


def replace_sequence_chars(sequence_name):
    import re
    return re.sub('[ /;*():]', '_', sequence_name)


def fill_zeros(s, length):
    return ('0' * (length - len(str(s)))) + str(s)


def days_between(d1, d2):
    '''
    Calculate the days between two dates

    :param d1: date 1
    :param d2: date 2
    :return: number of days between date 2 and date 1
    '''
    from datetime import datetime
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


def viscode_to_session(viscode):
    '''
    Replace the session label 'bl' with 'M00' or capitalize the session name passed
    as input.

    :param viscode: session name
    :return: M00 if is the baseline session or the original session name capitalized
    '''
    if viscode == 'bl':
        return 'M00'
    else:
        return viscode.capitalize()


def center_nifti_origin(input_image, output_image):
    '''
     Put the origin of the coordinate system at the center of the image.

    :param input_image: path to the input image
    :param output_image: path to the output image (where the result will be stored)
    :return:
    '''
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






