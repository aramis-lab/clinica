# coding: utf-8

"""
 Module for converting AV45 and Florbetaben PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_av45_fbb_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert AV-45 and Florbetaben PET images of ADNI into BIDS format

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

    cprint('Calculating paths of AV45 and Florbetaben PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_av45_fbb_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of AV45 and Florbetaben PET images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 'av45_fbb')
    cprint(Fore.GREEN + 'AV45 and Florbetaben PET conversion done.' + Fore.RESET)


def compute_av45_fbb_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to the AV45 and Florbetaben PET images and store them in a tsv file

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
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars, find_image_path
    from clinica.utils.stream import cprint
    from functools import reduce

    pet_amyloid_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID', 'Series_ID',
                       'Image_ID', 'Original', 'Tracer']
    pet_amyloid_df = pd.DataFrame(columns=pet_amyloid_col)

    # Loading needed .csv files
    av45qc = pd.read_csv(path.join(csv_dir, 'AV45QC.csv'), sep=',', low_memory=False)
    amyqc = pd.read_csv(path.join(csv_dir, 'AMYQC.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.shape[0] < 1:
            # TODO Log somewhere subjects without Amyloid PET images
            # cprint('No Amyloid PET images metadata for subject - ' + subj)
            continue

        # QC for AV45 PET images for ADNI 1, GO and 2
        av45_qc_subj = av45qc[(av45qc.PASS == 1) & (av45qc.RID == int(subj[-4:]))]

        # QC for Amyloid PET images for ADNI 3
        amy_qc_subj = amyqc[(amyqc.SCANQLTY == 1) & (amyqc.RID == int(subj[-4:]))]
        amy_qc_subj.insert(0, 'EXAMDATE', amy_qc_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        amyloid_qc_subj = pd.concat([av45_qc_subj, amy_qc_subj], axis=0, ignore_index=True, sort=False)

        for visit in list(amyloid_qc_subj.VISCODE2.unique()):
            amyloid_qc_visit = amyloid_qc_subj[amyloid_qc_subj.VISCODE2 == visit]

            # TODO Check
            # If there are several scans for a timepoint we keep image acquired last (higher LONIUID)
            amyloid_qc_visit = amyloid_qc_visit.sort_values("LONIUID", ascending=False)

            qc_visit = amyloid_qc_visit.iloc[0]

            # Corresponding LONI image ID for original scan in PET Meta List
            int_image_id = int(qc_visit.LONIUID[1:])

            original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                 & (subject_pet_meta['Image ID'] == int_image_id)
                                                 & (subject_pet_meta.Sequence.map(
                                                        lambda s: (s.lower().find('early') < 0)))]

            # If no corresponding Amyloid PET metadata for an original image,
            # take scan at the same date
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)
                                                     & (subject_pet_meta.Sequence.map(
                                                            lambda s: (s.lower().find('early') < 0)))]

                if original_pet_meta.shape[0] < 1:
                    # TODO Log somewhere QC visits without image metadata
                    cprint('No Amyloid PET images metadata for subject - ' + subj + ' for visit ' + qc_visit.VISCODE2)
                    continue

            original_image = original_pet_meta.iloc[0]

            # To determine type of amyloid PET tracer we find the
            # Coreg, Avg, Std Img and Vox Siz, Uniform Resolution image
            # with the same Series ID of the original image
            final_pet_meta = subject_pet_meta[(subject_pet_meta.Sequence.map(
                lambda x: (x.find('Coreg, Avg, Std Img and Vox Siz, Uniform Resolution') > 0)))
                                              & (subject_pet_meta['Series ID'] == original_image['Series ID'])]

            if final_pet_meta.shape[0] < 1:
                final_pet_meta = subject_pet_meta[(subject_pet_meta.Sequence.map(
                    lambda x: (x.find('Coreg, Avg, Std Img and Vox Siz, Uniform Resolution') > 0)))
                                                     & (subject_pet_meta['Scan Date'] == original_image['Scan Date'])]
                if final_pet_meta.shape[0] < 1:
                    # TODO Log
                    cprint('No "Coreg, Avg, Std Img and Vox Siz, Uniform Resolution" Amyloid PET image metadata for subject'
                           ' ' + subj + ' for visit ' + qc_visit.VISCODE2)
                    continue

            processed_sequence = final_pet_meta.iloc[0].Sequence

            if processed_sequence.startswith('AV45'):
                tracer = 'AV45'
            elif processed_sequence.startswith('FBB'):
                tracer = 'FBB'
            else:
                # TODO Log
                cprint('Unknown tracer for Amyloid PET image metadata for subject ' + subj +
                       ' for visit ' + qc_visit.VISCODE2)
                continue

            # Co-registered and Averaged image with the same Series ID of the original image
            averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == '%s Co-registered, Averaged' % tracer)
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
                  str(image_id), original, tracer]],
                columns=pet_amyloid_col)
            pet_amyloid_df = pet_amyloid_df.append(row_to_append, ignore_index=True)

    # TODO check for new exceptions in ADNI3
    # Exceptions
    # ==========
    conversion_errors = [  # Eq_1
                         ('128_S_2220', 'm48')]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((pet_amyloid_df.Subject_ID == conv_error[0])
                             & (pet_amyloid_df.VISCODE == conv_error[1]))

    if error_indices:
        indices_to_remove = pet_amyloid_df.index[reduce(operator.or_, error_indices, False)]
        pet_amyloid_df.drop(indices_to_remove, inplace=True)

    # DONE - Make a function reusing this code for different modalities
    # TODO check if it works properly

    images = find_image_path(pet_amyloid_df, source_dir, 'Amyloid', 'I', 'Image_ID')

    amyloid_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(amyloid_csv_path):
        os.mkdir(amyloid_csv_path)
    images.to_csv(path.join(amyloid_csv_path, 'amyloid_pet_paths.tsv'), sep='\t', index=False)

    return images


def check_exceptions(bids_dir):
    from os import path
    import pandas as pd
    from glob import glob

    amyloid_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'amyloid_pet_paths.tsv'), sep='\t')

    av45_paths = amyloid_paths[amyloid_paths.Tracer == 'AV45']

    av45_paths = av45_paths[av45_paths.Path.notnull()]

    av45_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in av45_paths.Subject_ID.to_list()]
    av45_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in av45_paths.VISCODE.to_list()]

    count = 0

    for r in av45_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'pet')
        image_pattern = path.join(image_dir, '%s_%s_*av45*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)
        if not files_list:
            print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 1:
            print("Too many images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)

    print(count)

    fbb_paths = amyloid_paths[amyloid_paths.Tracer == 'FBB']

    fbb_paths = fbb_paths[fbb_paths.Path.notnull()]

    fbb_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in fbb_paths.Subject_ID.to_list()]
    fbb_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in fbb_paths.VISCODE.to_list()]

    count = 0

    for r in fbb_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'pet')
        image_pattern = path.join(image_dir, '%s_%s_*fbb*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)
        if not files_list:
            print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 1:
            print("Too many images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)

    print(count)
