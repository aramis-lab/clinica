# coding: utf-8

"""
 Module for converting Tau PET of ADNI
"""

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def convert_adni_tau_pet(source_dir, csv_dir, dest_dir, subjs_list=None):
    """
    Convert Tau PET images of ADNI into BIDS format

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

    cprint('Calculating paths of TAU PET images. Output will be stored in %s.' % path.join(dest_dir, 'conversion_info'))
    images = compute_tau_pet_paths(source_dir, csv_dir, dest_dir, subjs_list)
    cprint('Paths of TAU PET images found. Exporting images into BIDS ...')
    t1_pet_paths_to_bids(images, dest_dir, 'tau')
    cprint(Fore.GREEN + 'TAU PET conversion done.' + Fore.RESET)


def compute_tau_pet_paths(source_dir, csv_dir, dest_dir, subjs_list):
    """
    Compute the paths to Tau PET images

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        subjs_list: subjects list

    Returns: pandas Dataframe containing the path for each Tau PET image

    """

    import pandas as pd
    import os
    import operator
    from os import path
    from functools import reduce
    from clinica.iotools.converters.adni_to_bids.adni_utils import replace_sequence_chars, find_image_path
    from clinica.utils.stream import cprint

    pet_tau_col = ['Phase', 'Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                   'Series_ID', 'Image_ID', 'Original']
    pet_tau_df = pd.DataFrame(columns=pet_tau_col)

    # Loading needed .csv files
    tauqc = pd.read_csv(path.join(csv_dir, 'TAUQC.csv'), sep=',', low_memory=False)
    tauqc3 = pd.read_csv(path.join(csv_dir, 'TAUQC3.csv'), sep=',', low_memory=False)
    pet_meta_list = pd.read_csv(path.join(csv_dir, 'PET_META_LIST.csv'), sep=',', low_memory=False)

    for subj in subjs_list:

        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]

        if subject_pet_meta.shape[0] < 1:
            # TODO Log somewhere subjects without TAU PET images metadata
            continue

        # QC for TAU PET images for ADNI 2
        tau_qc2_subj = tauqc[(tauqc.SCANQLTY == 1) & (tauqc.RID == int(subj[-4:]))]

        # QC for TAU PET images for ADNI 3
        tau_qc3_subj = tauqc3[(tauqc3.SCANQLTY == 1) & (tauqc3.RID == int(subj[-4:]))]

        # Concatenating visits in both QC files
        tau_qc_subj = pd.concat([tau_qc2_subj, tau_qc3_subj], axis=0, ignore_index=True, sort=False)

        for visit in list(tau_qc_subj.VISCODE2.unique()):
            # TODO Infer visit from ADNIMERGE visits
            if str(visit) == 'nan':
                continue

            pet_qc_visit = tau_qc_subj[tau_qc_subj.VISCODE2 == visit]

            # If there are several scans for a timepoint we keep image acquired last (higher LONIUID)
            pet_qc_visit = pet_qc_visit.sort_values("LONIUID", ascending=False)

            qc_visit = pet_qc_visit.iloc[0]

            # Corresponding LONI image ID for original scan in PET Meta List
            int_image_id = int(qc_visit.LONIUID[1:])

            original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                 & (subject_pet_meta['Image ID'] == int_image_id)]

            # If no corresponding TAU PET metadata for an original image,
            # take scan at the same date containing TAU or AV 45 in sequence name
            if original_pet_meta.shape[0] < 1:
                original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                     & subject_pet_meta.Sequence.map(
                    lambda x: ((x.lower().find('tau') > -1) |
                               (x.lower().find('av-1451') > -1) |
                               (x.lower().find('av1451') > -1)))
                                                        & (subject_pet_meta['Scan Date'] == qc_visit.SCANDATE)]

                if original_pet_meta.shape[0] < 1:
                    # TODO Log somewhere QC visits without image metadata
                    cprint('No TAU-PET images metadata for subject - ' + subj + ' for visit ' + qc_visit.VISCODE2)
                    continue

            original_image = original_pet_meta.iloc[0]

            # Co-registered and Averaged image with the same Series ID of the original image
            averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == 'AV1451 Co-registered, Averaged')
                                                 & (subject_pet_meta['Series ID'] == original_image['Series ID'])]

            # If an explicit AV1451 Co-registered, Averaged image does not exist,
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
                columns=pet_tau_col)
            pet_tau_df = pet_tau_df.append(row_to_append, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # Multiple output images
                         ('098_S_4275', 'm84')]

    error_indices = []
    for conv_error in conversion_errors:
        error_indices.append((pet_tau_df.Subject_ID == conv_error[0])
                             & (pet_tau_df.VISCODE == conv_error[1]))

    if error_indices:
        indices_to_remove = pet_tau_df.index[reduce(operator.or_, error_indices, False)]
        pet_tau_df.drop(indices_to_remove, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(pet_tau_df, source_dir, 'TAU', 'I', 'Image_ID')

    tau_csv_path = path.join(dest_dir, 'conversion_info')
    if not os.path.exists(tau_csv_path):
        os.mkdir(tau_csv_path)
    images.to_csv(path.join(tau_csv_path, 'tau_pet_paths.tsv'), sep='\t', index=False)

    return images


def check_exceptions(bids_dir):
    from os import path
    import pandas as pd
    from glob import glob

    tau_paths = pd.read_csv(path.join(bids_dir, 'conversion_info', 'tau_pet_paths.tsv'), sep='\t')

    tau_paths = tau_paths[tau_paths.Path.notnull()]

    tau_paths['BIDS_SubjID'] = ['sub-ADNI' + s.replace('_', '') for s in tau_paths.Subject_ID.to_list()]
    tau_paths['BIDS_Session'] = ['ses-' + s.replace('bl', 'm00').upper() for s in tau_paths.VISCODE.to_list()]

    count = 0

    for r in tau_paths.iterrows():
        image = r[1]
        image_dir = path.join(bids_dir, image.BIDS_SubjID, image.BIDS_Session, 'pet')
        image_pattern = path.join(image_dir, '%s_%s_*tau*' % (image.BIDS_SubjID, image.BIDS_Session))
        files_list = glob(image_pattern)
        if not files_list:
            print("No images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            count += 1

        elif len(files_list) > 1:
            print("Too many images for subject %s in session %s" % (image.BIDS_SubjID, image.BIDS_Session))
            print(files_list)

    print(count)
