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
    Convert FDG PET images of ADNI into BIDS format

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    """

    import pandas as pd
    from os import path
    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids.adni_utils import t1_pet_paths_to_bids
    from colorama import Fore

    if subjs_list is None:
        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path, sep=',', low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint('Calculating paths of FDG PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_fdg_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of FDG PET images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 'fdg')
    cprint(Fore.GREEN + 'FDG PET conversion done.' + Fore.RESET)


def compute_fdg_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to the FDG PET images and store them in a tsv file

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns:
        images: a dataframe with all the paths to the PET images that will be converted into BIDS

    """

    import pandas as pd
    import os
    import operator
    from os import path
    from functools import reduce
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars, find_image_path
    from clinica.utils.stream import cprint

    pet_fdg_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']
    pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)

    # Loading needed .csv files
    petqc = pd.read_csv(path.join(csv_dir, 'PETQC.csv'), sep=',', low_memory=False)
    petqc3 = pd.read_csv(path.join(csv_dir, 'PETC3.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.shape[0] < 1:
            # TODO Log somewhere subjects without FDG-PET images metadata
            continue

        # QC for FDG PET images for ADNI 1, GO and 2
        pet_qc_1go2_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subj[-4:]))]

        # QC for FDG PET images for ADNI 3
        pet_qc3_subj = petqc3[(petqc3.SCANQLTY == 1) & (petqc3.RID == int(subj[-4:]))]
        pet_qc3_subj.insert(0, 'EXAMDATE', pet_qc3_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        pet_qc_subj = pd.concat([pet_qc_1go2_subj, pet_qc3_subj], axis=0, ignore_index=True, sort=False)

        for visit in list(pet_qc_subj.VISCODE2.unique()):
            pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE2 == visit]

            # If there are several scans for a timepoint we keep image acquired last (higher LONIUID)
            pet_qc_visit = pet_qc_visit.sort_values("LONIUID", ascending=False)

            qc_visit = pet_qc_visit.iloc[0]

            # Corresponding LONI image ID for original scan in PET Meta List
            int_image_id = int(qc_visit.LONIUID[1:])

            original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                 & (subject_pet_meta['Image ID'] == int_image_id)
                                                 & (subject_pet_meta.Sequence.map(lambda s: (s.lower().find('early')
                                                                                             < 0)))]
            # If no corresponding FDG PET metadata for an original image,
            # take scan at the same date containing FDG in sequence name
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & (subject_pet_meta.Sequence.map(lambda x: (x.lower().find('fdg')
                                                                                                 > -1)))
                                                     & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]

                if original_pet_meta.shape[0] < 1:
                    # TODO Log somewhere QC visits without image metadata
                    cprint('No FDG-PET images metadata for subject - ' + subj + ' for visit ' + qc_visit.VISCODE2)
                    continue

            original_image = original_pet_meta.iloc[0]

            # Co-registered and Averaged image with the same Series ID of the original image
            averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == 'Co-registered, Averaged')
                                                 & (subject_pet_meta['Series ID'] == original_image['Series ID'])]

            # If an explicit Co-registered, Averaged image does not exist,
            # the original image is already in that preprocessing stage.

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
                [[qc_visit.Phase, subj, qc_visit.VISCODE2, str(visit), sequence, date, str(study_id), str(series_id),
                  str(image_id), original]],
                columns=pet_fdg_col)
            pet_fdg_df = pet_fdg_df.append(row_to_append, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # NONAME.nii
                         ('031_S_0294', 'bl'),
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

    images = find_image_path(pet_fdg_df, source_dir, 'FDG', 'I', 'Image_ID')

    fdg_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(fdg_csv_path):
        os.mkdir(fdg_csv_path)
    images.to_csv(path.join(fdg_csv_path, 'fdg_pet_paths.tsv'), sep='\t', index=False)

    return images


def check_exceptions(bids_dir):
    """
    Function for verifying that proper images are in bids folder
    """

    from os import path
    import pandas as pd
    from glob import glob

    fdg_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'fdg_pet_paths.tsv'), sep='\t')
    fdg_paths = fdg_paths[fdg_paths.Path.notnull()]
    fdg_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in fdg_paths.Subject_ID.to_list()]
    fdg_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in fdg_paths.VISCODE.to_list()]

    count = 0
    for r in fdg_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'pet')
        image_pattern = path.join(image_dir, '%s_%s_*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)
        if not files_list:
            print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 1:
            print("Too many images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)

    print(count)
