# coding: utf-8

"""
 Module for converting PIB PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_pib_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert PIB PET images of ADNI into BIDS format

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

    cprint('Calculating paths of PIB PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_pib_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of PIB PET images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 'pib')
    cprint(Fore.GREEN + 'PIB PET conversion done.' + Fore.RESET)


def compute_pib_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to the PIB PET images and store them in a tsv file

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

    pet_pib_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']
    pet_pib_df = pd.DataFrame(columns=pet_pib_col)

    # Loading needed .csv files
    pibqc = pd.read_csv(path.join(csv_dir, 'PIBQC.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.shape[0] < 1:
            # TODO Log somewhere subjects without PIB PET metadata
            continue

        # QC for PIB PET images
        pet_qc_subj = pibqc[(pibqc.PASS == 1) & (pibqc.RID == int(subj[-4:]))]

        for visit in list(pet_qc_subj.VISCODE.unique()):
            pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE == visit]

            qc_visit = pet_qc_visit.iloc[0]

            # Corresponding LONI image ID for original scan in PET Meta List
            int_image_id = int(qc_visit.LONIUID[1:])

            original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                 & (subject_pet_meta['Image ID'] == int_image_id)
                                                 & (subject_pet_meta.Sequence.map(lambda s: (s.lower().find('early')
                                                                                             < 0)))]
            # If no corresponding PIB PET metadata for an original image,
            # take scan at the same date containing PIB in sequence name
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & (subject_pet_meta.Sequence.map(lambda x: (x.lower().find('pib')
                                                                                                 > -1)))
                                                     & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]

                if original_pet_meta.shape[0] < 1:
                    # TODO Log somewhere QC visits without
                    cprint('No PIB-PET images metadata for subject - ' + subj + ' for visit ' + qc_visit.VISCODE)
                    continue

            original_image = original_pet_meta.iloc[0]

            # Co-registered and Averaged image with the same Series ID of the original image
            averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == 'PIB Co-registered, Averaged')
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
                [['ADNI1', subj, qc_visit.VISCODE, str(visit), sequence, date, str(study_id), str(series_id),
                  str(image_id), original]],
                columns=pet_pib_col)
            pet_pib_df = pet_pib_df.append(row_to_append, ignore_index=True)

    # TODO check for exceptions
    # Exceptions
    # ==========
    conversion_errors = []

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((pet_pib_df.Subject_ID == conv_error[0])
                             & (pet_pib_df.VISCODE == conv_error[1]))

    if error_indices:
        indices_to_remove = pet_pib_df.index[reduce(operator.or_, error_indices, False)]
        pet_pib_df.drop(indices_to_remove, inplace=True)

    # DONE - Make a function reusing this code for different modalities
    # TODO check if it works properly

    images = find_image_path(pet_pib_df, source_dir, 'PIB', 'I', 'Image_ID')

    pib_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(pib_csv_path):
        os.mkdir(pib_csv_path)
    images.to_csv(path.join(pib_csv_path, 'pib_pet_paths.tsv'), sep='\t', index=False)

    return images


def check_exceptions(bids_dir):
    from os import path
    import pandas as pd
    from glob import glob

    pib_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'pib_pet_paths.tsv'), sep='\t')

    pib_paths = pib_paths[pib_paths.Path.notnull()]

    pib_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in pib_paths.Subject_ID.to_list()]
    pib_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in pib_paths.VISCODE.to_list()]

    count = 0

    for r in pib_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'pet')
        image_pattern = path.join(image_dir, '%s_%s_*pib*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)
        if not files_list:
            print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 1:
            print("Too many images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)

    print(count)
