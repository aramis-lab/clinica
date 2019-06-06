# coding: utf-8

"""
 Module for converting FDG PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_fdg_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Args:
        source_dir:
        csv_dir:
        dest_dir:
        subjs_list:

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

    cprint('Calculating paths of FDG PET images. Output will be stored in ' + path.join(dest_dir, 'conversion_info') + '.')
    images = compute_fdg_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of FDG PET images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 'fdg')
    cprint(Fore.GREEN + 'FDG PET conversion done.' + Fore.RESET)


def compute_fdg_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """

    Args:
        source_dir:
        csv_dir:
        dest_dir:
        subjs_list:

    Returns:

    """
    import pandas as pd
    import os
    import operator
    from os import walk, path
    from numpy import argsort
    # from clinica.iotools.converters.adni_utils import replace_sequence_chars
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars
    from clinica.utils.stream import cprint
    from functools import reduce

    pet_fdg_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']

    pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
    petqc_path = path.join(csv_dir, 'PETQC.csv')
    pet_meta_list_path = path.join(csv_dir, 'PET_META_LIST.csv')
    petqc = pd.io.parsers.read_csv(petqc_path, sep=',')
    pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')

    for subj in subjs_list:
        pet_qc_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subj[-4:]))]
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]
        if subject_pet_meta.shape[0] < 1:
            # TODO Log somewhere subjects with problems
            cprint('NO Screening: Subject - ' + subj)
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
                    # TODO Log somewhere subjects with problems
                    cprint('No regular FDG-PET image: Subject - ' + subj + ' for visit ' + visit)
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
                (subject_pet_meta['Orig/Proc'] == 'Original') & (subject_pet_meta['Image ID'] == int_image_id)
                & (subject_pet_meta.Sequence.map(lambda s: (s.lower().find('early') < 0)))]
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[
                        (subject_pet_meta['Orig/Proc'] == 'Original')
                        & (subject_pet_meta.Sequence.map(
                            lambda x: (x.lower().find('fdg') > -1))
                           )
                        & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]
                if original_pet_meta.shape[0] < 1:
                    # TODO Log somewhere subjects with problems
                    cprint('NO Screening: Subject - ' + subj + ' for visit ' + qc_visit.VISCODE2)
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
            sequence = replace_sequence_chars(sel_image.Sequence)
            date = sel_image['Scan Date']
            study_id = sel_image['Study ID']
            series_id = sel_image['Series ID']
            image_id = sel_image['Image ID']

            row_to_append = pd.DataFrame(
                [[subj, qc_visit.VISCODE2, str(visit), sequence, date, str(study_id), str(series_id), str(image_id),
                  original]],
                columns=pet_fdg_col)
            pet_fdg_df = pet_fdg_df.append(row_to_append, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # NONAME.nii
                         ('037_S_1421', 'm36'),
                         ('037_S_1078', 'm36'),
                         # Empty folders
                         ('941_S_1195', 'm48'),
                         ('005_S_0223', 'm12')]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((pet_fdg_df.Subject_ID == conv_error[0])
                             & (pet_fdg_df.VISCODE == conv_error[1]))

    indices_to_remove = pet_fdg_df.index[reduce(operator.or_, error_indices, False)]
    pet_fdg_df.drop(indices_to_remove, inplace=True)

    images = pet_fdg_df
    # count = 0
    # total = images.shape[0]
    is_dicom = []
    image_folders = []
    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)
        # count += 1
        # print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
        #     count) + ' / ' + str(total)
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

        dicom = True
        for (dirpath, dirnames, filenames) in walk(image_path):
            for f in filenames:
                if f.endswith(".nii"):
                    dicom = False
                    image_path = path.join(dirpath, f)
                    break

        is_dicom.append(dicom)
        image_folders.append(image_path)
        if image_path == '':
            cprint('No FDG image path found for subject ' + str(image.Subject_ID))

    images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
    images.loc[:, 'Path'] = pd.Series(image_folders, index=images.index)

    fdg_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(fdg_csv_path):
        os.mkdir(fdg_csv_path)
    images.to_csv(path.join(fdg_csv_path, 'fdg_pet_paths.tsv'), sep='\t', index=False)

    return images
