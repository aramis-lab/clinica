# coding: utf-8

"""
Methods used by BIDS converters
"""
__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Completed"


# -- Methods for the clinical data --
# @ToDo:test this function
def create_participants_df(study_name, clinical_spec_path, clinical_data_dir, bids_ids, delete_non_bids_info=True):
    """
    Create the file participants.tsv

    Args:
        study_name: name of the study (Ex. ADNI)
        clinical_spec_path: path to the clinical file
        clinical_data_dir: path to the directory where the clinical data are
        stored
        bids_ids: list of bids ids
        delete_non_bids_info: if True delete all the rows of the subjects that
        are not available in the BIDS dataset

    Returns: a pandas dataframe that contains the participants data
    """
    import pandas as pd
    import os
    from os import path
    import numpy as np
    from clinica.utils.stream import cprint

    fields_bids = ['participant_id']
    prev_location = ''
    index_to_drop = []
    subjects_to_drop = []
    location_name = study_name + ' location'

    # Load the data from the clincal specification file
    participants_specs = pd.read_excel(clinical_spec_path, sheetname='participant.tsv')
    participant_fields_db = participants_specs[study_name]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs['BIDS CLINICA']

    # Extract the list of the available BIDS fields for the dataset
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])

    # Init the dataframe that will be saved in the file participants.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

    for i in range(0, len(participant_fields_db)):
        # If a field not empty is found
        if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            tmp = field_location[i].split('/')
            location = tmp[0]
            # If a sheet is available
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ''
            # Check if the file to open for a certain field is the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

                if file_ext == '.xlsx':
                    file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
                elif file_ext == '.csv':
                    file_to_read = pd.read_csv(file_to_read_path)
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                value_to_read = file_to_read[participant_fields_db[i]]
                if participant_fields_bids[i] == 'alternative_id_1' and\
                        (value_to_read.dtype == np.float64 or value_to_read.dtype == np.int64):
                    if not pd.isnull(file_to_read.get_value(j, participant_fields_db[i])):
                        value_to_append = str(file_to_read.get_value(j, participant_fields_db[i])).rstrip('.0')
                    else:
                        value_to_append = np.NaN
                else:
                    value_to_append = file_to_read.get_value(j, participant_fields_db[i])
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    if study_name == 'ADNI' or study_name == 'AIBL':
        # ADNImerge contains one row for each visits so there are duplicates
        participant_df = participant_df.drop_duplicates(subset=['alternative_id_1'], keep='first')
    participant_df.reset_index(inplace=True, drop=True)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        if study_name == 'OASIS':
            value = (participant_df['alternative_id_1'][i].split("_"))[1]
        else:
            value = remove_space_and_symbols(participant_df['alternative_id_1'][i])

        bids_id = [s for s in bids_ids if value in s]

        if len(bids_id) == 0:
            index_to_drop.append(i)
            subjects_to_drop.append(value)
        else:
            participant_df.set_value(i, 'participant_id', bids_id[0])

    if len(subjects_to_drop) > 0:
        cprint('The following subjects of ADNIMERGE were not found in your BIDS folder :\n'
               + ', '.join(subjects_to_drop))
    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.drop(index_to_drop)

    return participant_df


def create_sessions_dict(clinical_data_dir, study_name, clinical_spec_path, bids_ids, name_column_ids, subj_to_remove=[]):
    """

    Extract the information regarding the sessions and store them in a dictionary (session M00 only)

    Args:
        clinical_data_dir: path to the input folder
        study_name: name of the study (Ex: ADNI)
        clinical_spec_path:  path to the clinical file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject ids are stored
        subj_to_remove: subjects to remove

    """
    import pandas as pd
    from os import path
    import numpy as np
    import os

    # Load data
    location = study_name + ' location'
    sessions = pd.read_excel(clinical_spec_path, sheetname='sessions.tsv')
    sessions_fields = sessions[study_name]
    field_location = sessions[location]
    sessions_fields_bids = sessions['BIDS CLINICA']
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    sessions_df = pd.DataFrame(columns=fields_bids)

    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i].split('/')
            location = tmp[0]
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ''

            file_to_read_path = path.join(clinical_data_dir, location)
            file_ext = os.path.splitext(location)[1]

            if file_ext == '.xlsx':
                file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
            elif file_ext == '.csv':
                file_to_read = pd.read_csv(file_to_read_path)

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]
                # Extracts the subject ids columns from the dataframe
                subj_id = row[name_column_ids]
                if hasattr(subj_id, 'dtype'):
                    if subj_id.dtype == np.int64:
                        subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)
                if study_name == 'OASIS':
                    subj_id_alpha = str(subj_id[0:3] + 'IS' + subj_id[3] + subj_id[5:9])

                # Extract the corresponding BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if len(subj_bids) == 0:
                    # If the subject is not an excluded one
                    if subj_id not in subj_to_remove:
                        print(sessions_fields[i] + ' for ' + subj_id + ' not found in the BIDS converted.')
                else:
                    subj_bids = subj_bids[0]
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]
                    if subj_bids in sessions_dict:
                        (sessions_dict[subj_bids]['M00']).update({sessions_fields_bids[i]: row[sessions_fields[i]]})
                    else:
                        sessions_dict.update({subj_bids: {
                            'M00': {'session_id': 'ses-M00', sessions_fields_bids[i]: row[sessions_fields[i]]}}})

    return sessions_dict


def create_scans_dict(clinical_data_dir, study_name, clinic_specs_path, bids_ids, name_column_ids):
    """

    Args:
        clinical_data_dir: path to the directory where the clinical data are stored
        study_name: name of the study (Ex ADNI)
        clinic_specs_path: path to the clinical specification file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject id is contained

    Returns:

    """
    import pandas as pd
    from os import path
    scans_dict = {}
    prev_file = ''
    prev_sheet = ''

    if study_name not in get_supported_dataset():
        raise Exception('Dataset not supported. Supported datasets are:', get_supported_dataset())

    # Init the dictionary with the subject ids
    for bids_id in bids_ids:
        scans_dict.update({bids_id: {'T1/DWI/fMRI/FMAP': {}, 'FDG': {}}})

    scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
    fields_dataset = []
    fields_location = []
    fields_bids = []
    fields_mod = []

    # Extract the fields available and the corresponding bids name, location and type
    for i in range(0, len(scans_specs[study_name])):
        field = scans_specs[study_name][i]
        if not pd.isnull(field):
            fields_dataset.append(field)
            fields_bids.append(scans_specs['BIDS CLINICA'][i])
            fields_location.append(scans_specs[study_name + ' location'][i])
            fields_mod.append(scans_specs['Modalities related'][i])

    # For each field available extract the original name, extract from the file all the values and fill a data structure
    for i in range(0, len(fields_dataset)):
        # Location is composed by file/sheet
        location = fields_location[i].split('/')
        file_name = location[0]
        if len(location) > 1:
            sheet = location[1]
        else:
            sheet = ''
        # Check if the file to read is already opened
        if file_name == prev_file and sheet == prev_sheet:
            pass
        else:
            file_to_read_path = path.join(clinical_data_dir, file_name)
            if sheet != '':
                file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
            else:
                file_to_read = pd.read_excel(file_to_read_path)
            prev_file = file_name
            prev_sheet = sheet

        for bids_id in bids_ids:
            original_id = bids_id.replace('sub-'+study_name, '')
            row_to_extract = file_to_read[file_to_read[name_column_ids] == int(original_id)].index.tolist()
            if len(row_to_extract) > 0:
                row_to_extract = row_to_extract[0]
                # Fill the dictionary with all the information
                (scans_dict[bids_id][fields_mod[i]]).update(
                    {fields_bids[i]: file_to_read.iloc[row_to_extract][fields_dataset[i]]})
            else:
                print(" Scans information for " + bids_id + " not found.")

    return scans_dict


def write_sessions_tsv(bids_dir, sessions_dict):
    """
    Write the content of the function create scans dict in several tsv files following the BIDS specification

    Args:
        bids_dir: path to the bids directory
        sessions_dict: output of the function create_scans_dict

    Returns:

    """
    import os
    import pandas as pd
    from os import path
    from glob import glob

    bids_paths = glob(path.join(bids_dir, 'sub-*'))

    for sp in bids_paths:
        bids_id = sp.split(os.sep)[-1]

        if bids_id in sessions_dict:
            session_df = pd.DataFrame(sessions_dict[bids_id]['M00'], index=['i', ])
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False, encoding='utf8')
        else:
            print("No session data available for " + sp)
            session_df = pd.DataFrame(columns=['session_id'])
            session_df['session_id'] = pd.Series('M00')
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False, encoding='utf8')


def write_scans_tsv(bids_dir, bids_ids, scans_dict):
    """

    Write the scans dict into tsv files

    Args:
        bids_dir:  path to the BIDS directory
        bids_ids: list of bids ids
        scans_dict:  the output of the function create_scans_dict

    """
    import pandas as pd
    from os import path
    import os
    from glob import glob

    for bids_id in bids_ids:
        scans_df = pd.DataFrame()
        bids_id = bids_id.split(os.sep)[-1]
        # Create the file
        tsv_name = bids_id + "_ses-M00_scans.tsv"
        # If the file already exists, remove it
        if os.path.exists(path.join(bids_dir, bids_id, 'ses-M00', tsv_name)):
            os.remove(path.join(bids_dir, bids_id, 'ses-M00', tsv_name))

        mod_available = glob(path.join(bids_dir, bids_id, 'ses-M00', '*'))
        for mod in mod_available:
            mod_name = os.path.basename(mod)
            files = glob(path.join(mod, '*'))
            for file in files:
                file_name = os.path.basename(file)
                if mod_name == "anat" or mod_name == "dwi" or mod_name == "func":
                    f_type = 'T1/DWI/fMRI/FMAP'
                elif mod_name == 'pet':
                    f_type = 'FDG'

                row_to_append = pd.DataFrame(scans_dict[bids_id][f_type], index=[0])
                # Insert the column filename as first value
                row_to_append.insert(0, 'filename', path.join(mod_name, file_name))
                scans_df = scans_df.append(row_to_append)

            scans_df.to_csv(path.join(bids_dir, bids_id, 'ses-M00', tsv_name), sep='\t', index=False, encoding='utf8')


# -- Other methods --
def contain_dicom(folder_path):
    """
     Check if a folder contains DICOM images

    Args:
        folder_path: path to the folder

    Returns: True if dicom files are found inside the folder, False otherwise

    """
    from glob import glob
    from os import path

    dcm_files = glob(path.join(folder_path, '*.dcm'))
    if len(dcm_files) > 0:
        return True

    return False


def get_supported_dataset():
    """
    Return the list of supported datasets

    Returns: a list of supported datasets

    """
    return ['ADNI', 'CLINAD', 'PREVDEMALS', 'INSIGHT', 'OASIS', 'AIBL']


def dcm_to_nii(input_path, output_path, bids_name):
    """
    Convert DICOM to NIFTI trying to use several converters

    Args:
        input_path: path to the input folder with DICOM images
        output_path: path to the output folder
        bids_name: bids name to give to the converterd file

    """
    import os
    from os import path
    import subprocess
    from clinica.utils.stream import cprint
    from colorama import Fore

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if 'bold' in bids_name:
        # generation of the json file
        cmd = 'dcm2niix -b y -z y -o ' + output_path + ' -f ' + bids_name + ' ' + input_path
    else:
        cmd = 'dcm2niix -b n -z y -o ' + output_path + ' -f ' + bids_name + ' ' + input_path

    subprocess.run(cmd,
                   shell=True,
                   stderr=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL)

    # If dcm2niix didn't work use dcm2nii
    if not os.path.exists(path.join(output_path, bids_name + '.nii.gz')):
        cprint('\tConversion with dcm2niix failed, trying with dcm2nii')
        cmd = 'dcm2nii -a n -d n -e n -i y -g n -p n -m n -r n -x n -o ' + output_path + ' ' + input_path
        subprocess.run(cmd,
                       shell=True,
                       stderr=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL)

    # If the conversion failed with both tools
    if not os.path.exists(path.join(output_path, bids_name + '.nii.gz')):
        cprint(Fore.RED + 'Conversion of the dicom failed for ' + input_path + Fore.RESET)


def get_bids_subjs_list(bids_path):
    """

    Given a BIDS compliant dataset, returns the list of all the subjects available

    Args:/
        bids_path: path to the BIDS folder

    """
    import os
    from os import path

    return [d for d in os.listdir(bids_path) if os.path.isdir(path.join(bids_path, d))]


def get_bids_subjs_paths(bids_path):
    """

    Given a BIDS compliant dataset, returns the list of all paths to the subjects folders

    Args:
        bids_path: path to the BIDS folder
    """

    import os
    from os import path

    return [path.join(bids_path, d) for d in os.listdir(bids_path) if os.path.isdir(path.join(bids_path, d))]


def compute_new_subjects(original_ids, bids_ids):
    """
    Check for new subject to convert

    This function checks for news subjects to convert to the BIDS version i.e. subjects contained in the unorganised
    version that are not available in the bids version.


    Args:
        original_ids: list of all the ids of the unorganized folder.
        bids_ids: list of all the BIDS ids contained inside the bids converted version of the dataset

    Returns: a list containing the original_ids of the subjects that are not available in the bids converted version

    """
    to_return = []
    original_ids = remove_space_and_symbols(original_ids)

    for s in original_ids:
        if not any(s in id for id in bids_ids):
            to_return.append(s)

    return to_return


def remove_space_and_symbols(data):
    '''
    Remove spaces and  - _ from a list (or a single) of strings
    :param data: list of strings or a single string to clean
    :return:
        data: list of strings or a string without space and symbols _ and -
    '''

    import re

    if type(data) is list:
        for i in range(0, len(data)):
            data[i] = re.sub('[-_ ]', '', data[i])
    else:
        data = re.sub('[-_ ]', '', data)

    return data


def get_ext(file_path):

    import os

    root, ext = os.path.splitext(file_path)
    if ext in '.gz':
        file_ext = os.path.splitext(root)[1] + ext
    else:
        file_ext = ext
    return file_ext


def compress_nii(file_path):
    '''
    Compress nii files.

    :param file_path: path to the file to convert
    '''
    from os import remove
    import gzip
    import shutil

    with open(file_path, 'rb') as f_in:
        with gzip.open(file_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    remove(file_path)


def remove_rescan(list_path):
    """
    Remove all the folders containing the keyword 'rescan' from
    a given list of folders.

    Args:
        list_path (str): the list of the files to analize.

    Returns:
        The list of non rescanned folders.

    """

    import logging

    # Extract all the folders without the substring 'rescan'
    no_resc_lst = [s for s in list_path if 'rescan' not in s]
    if len(no_resc_lst) != len(list_path):
        for r_file in list(set(list_path)-set(no_resc_lst)):
            logging.warning('Rescan found '+r_file+' Ignored.')

    return no_resc_lst


def has_fixed(subj_id, ses, mod, special_list):
    """
    Check if a certain subject has some fixed folder name to convert for a certain modality.

    Args:
        subj_id: the subject to check
        ses: the session
        mod: the modality
        special_list: json file that contains the list of subjects and modality to convert

    Returns:
         False: if the fixed modality to convert is not available for the specified subject
         file_to_convert: name of the folder to choosefor the conversion
    """

    if subj_id in special_list:
        if mod in special_list[subj_id]:
            if ses in special_list[subj_id][mod]:
                file_to_convert = special_list[subj_id][mod]
                return file_to_convert
            else:
                return False
    return False


def choose_correction(dir, to_consider, mod):
    """
    Decides what is the best file type to choose depending on a priority list.

    The selection criteria is based on the input list toConsider that contains
    the files to consider in order of priority.

    Args:
        dir: the directory containing the files.
        to_consider: the list of files to consider in order of priority.
        mod: the modality

    Returns:
        The best file (according to the list to_consider) for the given modality or:
        -1 if if the modality is missing.
        0 if none of the desided corrections is available..

    """

    from os import path
    from glob import glob
    import os

    # Extract all the files available for a certain modality
    correction_list = remove_rescan(glob(path.join(dir, '*' + mod + '*')))
    if len(correction_list) == 0:
        return -1
    if len(correction_list) == 1:
        return correction_list[0].split(os.sep)[-1]
    else:
        for i in range(0, len(to_consider)):
            if any(to_consider[i] in c for c in correction_list):
                return to_consider[i]
        return 0


def get_bids_suff(mod):
    """
    Returns the BIDS suffix for a certain modality.

    Args:
        mod: modality.
    Returns:
        The suffix used in the BIDS standard for a certain modality.
    """

    bids_suff = {
        'T1': '_T1w',
        'T2': '_T2w',
        'Flair': '_FLAIR',
        'SingleMapPh': '_phasediff',
        'MultiMapPh': '_phase',
        'Map': '_magnitude',
        'fMRI': '_bold',
        'dwi': '_dwi',
        'pet': '_pet'
    }

    return bids_suff[mod]


def create_path_corr_file(out_dir):
    import os
    from os import path

    path_corr = path.join(out_dir, 'conversion_info', 'filename_correspondance.tsv')
    if os.path.exists(path_corr):
        os.remove(path_corr)

    f = open(path_corr, 'a')
    return f


def convert_T1(t1_path, output_path, t1_bids_name):
    """
    Convert into the BIDS specification a T1 image.

    Args:
        t1_path: the path of the T1 images to convert.
        output_path: output folder path.
        t1_bids_name: name to give to the file.

    """
    from os import path
    from shutil import copy
    import os

    bids_name = path.join(output_path, t1_bids_name + get_bids_suff('T1'))

    if contain_dicom(t1_path):
        print('DICOM found for t1 in ' + t1_path)
        dcm_to_nii(t1_path, output_path, t1_bids_name + get_bids_suff('T1'))
    else:
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        file_ext = get_ext(t1_path)
        bids_file = bids_name + file_ext
        print(bids_file)
        copy(t1_path, bids_file)
        # If  the original image is not compress, compress it
        if file_ext == '.nii':
            compress_nii(path.join(bids_file))


def convert_pet(folder_input, folder_output, pet_name, bids_name, task_name, acquisition=''):
    """
    Convert PET to BIDS

    Args:
        folder_input: path to the input folder
        folder_output: path to the BIDS folder
        pet_name: original name
        bids_name: bids name
        task_name: task name
        acquisition: acquisition name

    """

    from os import path
    from shutil import copy
    import os

    if acquisition != '':
        pet_bids_name = bids_name+'_task-'+task_name+'_acq-'+acquisition
    else:
        pet_bids_name = bids_name+'_task-'+task_name

    if contain_dicom(folder_input):
        print('DICOM found for PET')
        dcm_to_nii(folder_input, folder_output, pet_bids_name, 'pet')
    else:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)

        # If a prefixed pet is chosen for the conversion use it, otherwise extracts the pet file available into the folder
        if pet_name != '':
            pet_path = path.join(folder_input, pet_name)
        else:
            print('WARNING: feature to be implemented')

        file_ext = get_ext(pet_path)
        copy(pet_path, path.join(folder_output, pet_bids_name + get_bids_suff('pet') + file_ext))
        # If  the original image is not compress, compress it
        if file_ext == '.nii':
            compress_nii(path.join(folder_output, pet_bids_name + get_bids_suff('pet') + file_ext))


def convert_fieldmap(folder_input, folder_output, name, fixed_file=[False, False]):
    """
    Extract and convert into the BIDS specification fieldmap data.

    Args:
        folder_input: folder containing the fieldmap data.
        folder_output: folder where store the converted files.
        name: name of converted file.
        files_to_skip:

    Returns:
         -1 if the modality is not available.
         0 if the magnitude or the phase is missing (information incomplete).
    """

    from os import path
    from glob import glob
    from shutil import copy
    import nibabel as nib
    import os

    mag_missing = map_ph_missing = False

    # Check if there is a map or mapPh fixed to convert
    if fixed_file[0] is False and fixed_file[1] is False:
        map = remove_rescan(glob(path.join(folder_input, "*MAP_*", '*.nii.gz')))
        map_ph = remove_rescan(glob(path.join(folder_input, "*MAPph_*", '*.nii.gz')))
    elif fixed_file[0] is False and fixed_file[1] is not False:
        map = remove_rescan(glob(path.join(folder_input, "*MAP_*", '*.nii.gz')))
        map_ph = glob(path.join(folder_input, fixed_file[1], '*.nii.gz'))
    elif fixed_file[0] is not False and fixed_file[1] is False:
        map = glob(path.join(folder_input, fixed_file[0], '*.nii.gz'))
        map_ph = remove_rescan(glob(path.join(folder_input, "*MAPph_*", '*.nii.gz')))
    else:
        map = glob(path.join(folder_input, fixed_file[0], '*.nii.gz'))
        map_ph = glob(path.join(folder_input, fixed_file[1], '*.nii.gz'))

    files_to_skip = []
    if len(map) == 0:
        mag_missing = True
    if len(map_ph) == 0:
        map_ph_missing = True
    # If the information regarding the Fieldmap data are complete
    if len(map) > 0 and len(map_ph) > 0:
        map_ph_name = map_ph[0].split(os.sep)[-1]
        map_name = map[0].split(os.sep)[-1]

        # toSolve: there are some files that produce an error when loaded with Nibabel
        if (map_ph_name not in files_to_skip) and (map_name not in files_to_skip):
            # Open the files with Nibabel
            map_nib = nib.load(map[0])
            map_ph_nib = nib.load(map_ph[0])
            dim_map = (map_nib.header['dim'])[4]
            dim_map_ph = (map_ph_nib.header['dim'])[4]
            os.mkdir(path.join(folder_output))
            # Case 1: one phase difference image and at least one magnitude image
            if dim_map_ph == 1 and dim_map > 0:
                copy(map_ph[0], path.join(folder_output, name + get_bids_suff('SingleMapPh') + '.nii.gz'))
                os.system('fslsplit ' + map[0] + ' ' + path.join(folder_output, name + get_bids_suff('Map')))
                mag_list = glob(path.join(folder_output, name+get_bids_suff('Map')+'*'))

                for i in range(0, len(mag_list)):
                    old_mag_name = mag_list[i].split(os.sep)[-1]
                    # Remove the extension and the number sequence of fslsplit
                    new_mag_name = old_mag_name[:-11] + str(i+1)
                    os.rename(mag_list[i], path.join(folder_output, new_mag_name + ".nii.gz"))
            # Case 2: two phase images and two magnitude images'
            elif dim_map_ph == 2 and dim_map == 2:
                bids_name_ph = name + get_bids_suff('MultiMapPh')
                bids_name_mag = name + get_bids_suff('Map')
                os.system('fslsplit ' + map_ph[0] + ' ' + path.join(folder_output, bids_name_ph))
                os.system('fslsplit ' + map[0] + ' ' + path.join(folder_output, bids_name_mag))
                # Fslsplit produces files with a name not comply with BIDS specification, a rename is needed
                os.rename(path.join(folder_output, bids_name_mag + '0000.nii.gz'),
                          path.join(folder_output, bids_name_mag+'1.nii.gz'))
                os.rename(path.join(folder_output, bids_name_mag + '0001.nii.gz'),
                          path.join(folder_output, bids_name_mag + '2.nii.gz'))
                os.rename(path.join(folder_output, bids_name_ph + '0000.nii.gz'),
                          path.join(folder_output, bids_name_ph + '1.nii.gz'))
                os.rename(path.join(folder_output, bids_name_ph + '0001.nii.gz'),
                          path.join(folder_output, bids_name_ph + '2.nii.gz'))
    # The modalities is missing or incomplete
    else:
        if mag_missing and map_ph_missing:
            return -1
        else:
            return 0


def convert_flair(folder_input, folder_output, name, fixed_file=False):
    """
    Extracts and converts T2Flair data.

    Args:
        folder_input: folder containing T2Flair.
        folder_output: output folder path.
        name: name to give to the output file.

    Returns:
        -1 if no T2FLAIR is found in the input folder.
    """

    from os import path
    import logging
    from glob import glob
    from shutil import copy
    import os

    # If a given T2FLAIR is given for the conversion use it
    if fixed_file:
        fixed_flair_path = glob(path.join(folder_input, fixed_file, '*'))[0]
        copy(fixed_flair_path, path.join(folder_output, name + (get_bids_suff('Flair')) + '.nii.gz'))
    else:
        list_path = glob(path.join(folder_input, '*T2FLAIR*'))
        flair_lst = remove_rescan(list_path)
        if len(flair_lst) == 1:
            if not os.path.exists(folder_output):
                os.mkdir(folder_output)
            flair_path = glob(path.join(flair_lst[0], '*.nii.gz*'))[0]
            copy(flair_path, path.join(folder_output, name + (get_bids_suff('Flair')) + '.nii.gz'))
        elif len(flair_lst) == 0:
                return -1
        elif len(flair_lst) > 1:
                logging.warning('Multiple FLAIR found, computation aborted.')
                raise('Aborted')


def convert_fmri(folder_input, folder_output, name, fixed_fmri=False, task_name='rest'):
    """
    Extracts and converts into the BIDS specification fmri data.

    Args:
        folder_input: the folder containing the fmri to convert.
        folder_output: the output folder.
        name:

    Returns:
        -1 if no fMRI file is found in the input folder.
    """

    from os import path
    import logging
    from glob import glob
    from shutil import copy
    import os
    from clinica.utils.stream import cprint
    from colorama import Fore

    bids_name = name + '_task-' + task_name + get_bids_suff('fMRI')

    if contain_dicom(folder_input):
        dcm_to_nii(folder_input, folder_output, bids_name)
    else:
        if fixed_fmri is not False:
            fmri_lst = glob(path.join(folder_input, fixed_fmri))
        else:
            fmri_lst = remove_rescan(glob(path.join(folder_input, '*fMRI*')))
        if len(fmri_lst) > 0:
            if not os.path.exists(folder_output):
                os.mkdir(folder_output)
            fmri_file_path = glob(path.join(fmri_lst[0], '*.nii*'))[0]
            file_ext = get_ext(fmri_file_path)
            copy(fmri_file_path, path.join(folder_output, bids_name + file_ext))
            if file_ext == '.nii':
                cprint(Fore.YELLOW + 'Non compressed file found: ' + fmri_file_path + Fore.RESET)
                compress_nii(path.join(folder_output, bids_name + file_ext))
        else:
            logging.info('Non fMRI found for ' + folder_input)
            return -1


def merge_DTI(folder_input, folder_output, name, fixed_dti_list=False):
    """
    Merge all the DTI files of a given subject.

    For the merging only DTI folders containing all .nii.gz, .bval and .bvec are considered,
    otherwise the folder is ignored.

    Args:
        folder_input (str) : the folder containing all the DTI to merge.
        folder_output (str) : the folder where store the merged file.
        name (str): name to give to the merged file.

    Returns:
        -1 if the input folder doesn't contain any DTI folder.
        The list of incomplete DTI folders if there is some folders without bvec/bval/nii
    """

    from os import path
    from glob import glob
    import fileinput
    import os

    img = []
    bval = []
    bvec = []
    dti_list = []
    # The bvec output file has three rows that are the concatenation of all the merged bvec files
    lines_out_bvec = ['', '', '']
    # The bval output file has only one row that is the concatenation of all the merged bval files
    lines_out_bval = ['']

    if fixed_dti_list is not False:
        for dti in fixed_dti_list:
            dti_list.append(path.join(folder_input, dti))
    else:
        dti_list = remove_rescan(glob(path.join(folder_input, '*DTI*')))
    incomp_folders = []
    nr_dti = len(dti_list)
    if nr_dti == 0:
        return -1
    else:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)
        for folder in dti_list:
            if len(glob(path.join(folder, '*.bval'))) != 0 and len(glob(path.join(folder, '*.bvec'))) != 0:
                img.append(glob(path.join(folder, '*.nii*'))[0])
                bval.append(glob(path.join(folder, '*.bval'))[0])
                bvec.append(glob(path.join(folder, '*.bvec'))[0])
            else:
                incomp_folders.append(folder)

        # if it has been found at least a DTI folder complete with bvec, bval and nii.gz
        if len(img) > 0:
            file_suff = get_bids_suff('dwi')
            fin = fileinput.input(bval)
            # merge all the .nii.gz file with fslmerge
            os.system('fslmerge -t ' + path.join(folder_output, name + file_suff + '.nii.gz') + ' ' + " ".join(img))
            # merge all the .bval files
            fout = open(path.join(folder_output, name+file_suff+'.bval'), 'w')
            # Concatenate bval files
            for line in fin:
                if fileinput.isfirstline():
                    line_no = 0
                lines_out_bval[line_no] = lines_out_bval[line_no]+" "+line.rstrip()
                line_no += 1
            for i in range(0, len(lines_out_bval)):
                lines_out_bval[i] = lines_out_bval[i].lstrip()

                fout.write(lines_out_bval[i]+"\n")

            # Concatenate bvec files
            fin = fileinput.input(bvec)
            fout = open(path.join(folder_output, name + file_suff + '.bvec'), 'w')
            for line in fin:
                if fileinput.isfirstline():
                    line_no = 0

                lines_out_bvec[line_no] = lines_out_bvec[line_no]+" "+line.rstrip()
                line_no += 1
            for i in range(0, len(lines_out_bvec)):
                lines_out_bvec[i] = lines_out_bvec[i].lstrip()
                fout.write(lines_out_bvec[i] + "\n")

        if len(incomp_folders) > 0:
            return incomp_folders


def concatenate_bvec_bval(files_list, output_file, type):
    """
    Concatenate bvec and baval files

    Args:
        files_list: list of file
        output_file: path to the output file
        type: bvec of bval

    """
    import fileinput
    if type == 'bval':
        lines_out = ['']
    else:
        lines_out = ['', '', '']

    for line in files_list:
        if fileinput.isfirstline():
            line_no = 0
        lines_out[line_no] = lines_out[line_no] + " " + line.rstrip()
        line_no += 1
    for i in range(0, len(lines_out)):
        lines_out[i] = lines_out[i].lstrip()

        output_file.write(lines_out[i] + "\n")


def merge_noddi_dti(folder_input, folder_output, name):
    """
       Merge NODDI dti CATI organised following these rules:

       - DTI1, DTI3 and DTI5 are always Anterior to Posterior(AP/j-)
       - DTI2, DTI4, DTI6 are always Posterior to Anterior(PA/j)

    Args:
        folder_input: path to the folder where are the DTI
        folder_output: path to the BIDS folder
        name: BIDS id of the file
    """
    from glob import glob
    from os import path
    import fileinput
    import os

    dti_pa_nii = []
    dti_ap_nii = []
    dti_pa_bvec = []
    dti_ap_bvec = []
    dti_pa_bval = []
    dti_ap_bval = []

    dti_list = remove_rescan(glob(path.join(folder_input, '*DTI*')))

    if len(dti_list) != 6:
        raise Exception('Number of DTI found different from 6')

    if not os.path.exists(folder_output):
        os.mkdir(folder_output)

    out_file_name_pa = path.join(folder_output, name + '_seq-mshellPA' + get_bids_suff('dwi'))
    out_file_name_ap = path.join(folder_output, name + '_seq-mshellAP' + get_bids_suff('dwi'))

    dti_ap_paths = glob(path.join(folder_input, '*DTI[1,3,5]*'))
    dti_pa_paths = glob(path.join(folder_input, '*DTI[2,4,6]*'))

    for f in dti_pa_paths:
        dti_pa_nii.append(glob(path.join(f, '*.nii*'))[0])
        dti_pa_bvec.append(glob(path.join(f, '*.bvec'))[0])
        dti_pa_bval.append(glob(path.join(f, '*.bval'))[0])

    for f in dti_ap_paths:
        dti_ap_nii.append(glob(path.join(f, '*.nii*'))[0])
        dti_ap_bvec.append(glob(path.join(f, '*.bvec'))[0])
        dti_ap_bval.append(glob(path.join(f, '*.bval'))[0])

    # Merge all the nii files
    os.system('fslmerge -t ' + (out_file_name_pa + '.nii.gz') + ' ' + " ".join(dti_pa_nii))
    os.system('fslmerge -t ' + (out_file_name_ap + '.nii.gz') + ' ' + " ".join(dti_ap_nii))

    # Concatenate all the bvec/bval files
    concatenate_bvec_bval(fileinput.input(dti_pa_bval), open(path.join(out_file_name_pa + '.bval'), 'w'), 'bval')
    concatenate_bvec_bval(fileinput.input(dti_pa_bvec), open(path.join(out_file_name_pa + '.bvec'), 'w'), 'bvec')
    concatenate_bvec_bval(fileinput.input(dti_ap_bval), open(path.join(out_file_name_ap + '.bval'), 'w'), 'bval')
    concatenate_bvec_bval(fileinput.input(dti_ap_bvec), open(path.join(out_file_name_ap + '.bvec'), 'w'), 'bvec')
