# coding: utf-8

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def replace_sequence_chars(sequence_name):
    """
    Replace some special character with the sequence name given in input

    Args:
        sequence_name: sequence to replace

    Returns: the new string

    """
    import re
    return re.sub('[ /;*():]', '_', sequence_name)


def fill_zeros(s, length):
    return ('0' * (length - len(str(s)))) + str(s)


def days_between(d1, d2):
    """
    Calculate the days between two dates

    Args:
        d1: date 1
        d2: date 2

    Returns: number of days between date 2 and date 1

    """

    from datetime import datetime
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


def viscode_to_session(viscode):
    """
    Replace the session label 'bl' with 'M00' or capitalize the session name passed
    as input

    Args:
        viscode: session name

    Returns: M00 if is the baseline session or the original session name capitalized

    """

    if viscode == 'bl':
        return 'M00'
    else:
        return viscode.capitalize()


def center_nifti_origin(input_image, output_image):
    """

    Put the origin of the coordinate system at the center of the image

    Args:
        input_image: path to the input image
        output_image: path to the output image (where the result will be stored)

    Returns:

    """

    import nibabel as nib
    import numpy as np

    img = nib.load(input_image)
    canonical_img = nib.as_closest_canonical(img)
    hd = canonical_img.header
    # if hd['quatern_b'] != 0 or hd['quatern_c'] != 0 or hd['quatern_d'] != 0:
    #    print('Warning: Not all values in quatern are equal to zero')
    qform = np.zeros((4, 4))
    for i in range(1, 4):
        qform[i - 1, i - 1] = hd['pixdim'][i]
        qform[i - 1, 3] = -1.0 * hd['pixdim'][i] * hd['dim'][i] / 2.0
    new_img = nib.Nifti1Image(canonical_img.get_data(caching='unchanged'), qform)
    nib.save(new_img, output_image)


def remove_space_and_symbols(data):
    """
    Remove spaces and  - _ from a list (or a single) of strings

    Args:
        data: list of strings or a single string to clean

    Returns:
         data: list of strings or a string without space and symbols _ and -

    """
    import re

    if type(data) is list:
        for i in range(0, len(data)):
            data[i] = re.sub('[-_ ]', '', data[i])
    else:
        data = re.sub('[-_ ]', '', data)

    return data


def check_two_dcm_folder(dicom_path, bids_folder, image_uid):
    """
    Check if a folder contains more than one DICOM and if yes, copy the DICOM related to image id passed as parameter
    into a temporary folder called tmp_dicom_folder

    Args:
        dicom_path: path to the DICOM folder
        bids_folder: path to the BIDS folder where the dataset will be stored
        image_uid: image id of the fmri

    Returns: the path to the original DICOM folder or the path to a temporary folder called tmp_dicom_folder where only
     the DICOM to convert is copied

    """

    from glob import glob
    from os import path
    from shutil import copy
    import shutil
    import os

    temp_folder_name = 'tmp_dcm_folder_' + str(image_uid).strip(' ')
    dest_path = path.join(bids_folder, temp_folder_name)

    # Check if there is more than one xml file inside the folder
    xml_list = glob(path.join(dicom_path, '*.xml*'))

    if len(xml_list) > 1:
        # Remove the precedent tmp_dcm_folder if is existing
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        os.mkdir(dest_path)
        dmc_to_conv = glob(path.join(dicom_path, '*'+str(image_uid)+'.dcm*'))
        for d in dmc_to_conv:
            copy(d, dest_path)
        return dest_path
    else:
        return dicom_path


def remove_tmp_dmc_folder(bids_dir, image_id):
    """
    Remove the folder tmp_dmc_folder created by the method check_two_dcm_folder (if existing)

    Args:
        bids_dir: path to the BIDS directory

    Returns:

    """
    from os.path import exists, join
    from shutil import rmtree

    tmp_dcm_folder_path = join(bids_dir, 'tmp_dcm_folder_' + str(image_id).strip(' '))
    if exists(tmp_dcm_folder_path):
        rmtree(tmp_dcm_folder_path)


def check_bids_t1(bids_path, container='anat', extension='_T1w.nii.gz', subjects=None):

    import os

    if subjects is None:
        subjects = next(os.walk(bids_path))[1]

    errors = []
    print(len(subjects))
    for subject in subjects:
        sessions = next(os.walk(os.path.join(bids_path, subject)))[1]
        for session in sessions:
            image_name = subject + '_' + session + extension
            image_dir = os.path.join(bids_path, subject, session, container)
            if not os.path.isdir(image_dir):
                print('No directory: ' + image_dir)
            else:
                files = next(os.walk(image_dir))[2]
                if not files:
                    errors.append('Subject ' + subject + ' for session ' + session + ' folder is empty')
                for f in files:
                    if f != image_name:
                        errors.append('Subject ' + subject + ' for session ' + session + ' folder contains: ' + f)
    return errors


def check_bids_dwi(bids_path, container='dwi', extension=('_acq-axial_dwi.bvec', '_acq-axial_dwi.bval', '_acq-axial_dwi.nii.gz'), subjects=None):

    import os

    if subjects is None:
        subjects = next(os.walk(bids_path))[1]

    errors = []
    print(len(subjects))
    for subject in subjects:
        sessions = next(os.walk(os.path.join(bids_path, subject)))[1]
        for session in sessions:

            image_names = [subject + '_' + session + ext for ext in extension]
            image_names.sort()

            image_dir = os.path.join(bids_path, subject, session, container)

            if not os.path.isdir(image_dir):
                print('No directory: ' + image_dir)
            else:
                files = next(os.walk(image_dir))[2]
                files.sort()

                if not files:
                    errors.append('Subject ' + subject + ' for session ' + session + ' folder is empty')

                if image_names != files:
                    errors.append('Subject ' + subject + ' for session ' + session + ' folder contains: \n' + str(files))
    return errors


def is_nan(value):
    """
    Check if a value is nan

    Args:
        value: value that need to be checked

    Returns: true if is nan, false otherwise

    """

    from math import isnan
    import numpy as np

    if type(value) != float and type(value) != np.float64:
        return False

    if isnan(value):
        return True
    else:
        return False


def remove_fields_duplicated(bids_fields):
    seen = set()
    seen_add = seen.add
    return [x for x in bids_fields if not (x in seen or seen_add(x))]


def convert_diagnosis_code(diagnosis_code):
    """
    Convert the numeric field DXCURREN and DXCHANGE contained in DXSUM_PDXCONV_ADNIALL.csv into a code that identify the
    diagnosis

    Args:
        diagnosis_code: a string that represents a number between 1 and 9

    Returns:a code that identify a diagnosis

    """

    # Legenda
    diagnosis = {1: 'CN', 2: 'MCI', 3: 'AD', 4: 'MCI', 5: 'AD', 6: 'AD', 7: 'CN', 8: 'MCI', 9: 'CN'}

    if is_nan(diagnosis_code):
        return diagnosis_code
    else:
        return diagnosis[int(diagnosis_code)]


def write_adni_sessions_tsv(sessions_dict, fields_bids, bids_subjs_paths):
    """
    Write the result of method create_session_dict into several tsv files

    Args:
        sessions_dict: dictonary coming from the method create_sessions_dict
        fields_bids: fields bids to convert
        bids_subjs_paths: a list with the path to all bids subjects

    Returns:

    """
    from os import path
    import os
    import pandas as pd
    import numpy as np

    columns_order = remove_fields_duplicated(fields_bids)

    columns_order.insert(0, 'session_id')
    for sp in bids_subjs_paths:

        if not path.exists(sp):
            os.makedirs(sp)

        bids_id = sp.split(os.sep)[-1]

        fields_bids = list(set(fields_bids))
        sessions_df = pd.DataFrame(columns=fields_bids)

        if bids_id in sessions_dict:
            sess_aval = sessions_dict[bids_id].keys()
            for ses in sess_aval:
                sessions_df = sessions_df.append(pd.DataFrame(sessions_dict[bids_id][ses], index=['i', ]))

            sessions_df = sessions_df[columns_order]

            sessions_df[[
                'adas_Q1', 'adas_Q2', 'adas_Q3', 'adas_Q4',
                'adas_Q5', 'adas_Q6', 'adas_Q7', 'adas_Q8',
                'adas_Q9', 'adas_Q10', 'adas_Q11', 'adas_Q12',
                'adas_Q13'
                ]] = sessions_df[[
                    'adas_Q1', 'adas_Q2', 'adas_Q3', 'adas_Q4',
                    'adas_Q5', 'adas_Q6', 'adas_Q7', 'adas_Q8',
                    'adas_Q9', 'adas_Q10', 'adas_Q11', 'adas_Q12',
                    'adas_Q13']].apply(pd.to_numeric)

            sessions_df['adas_memory'] = (sessions_df['adas_Q1']
                                          + sessions_df['adas_Q4']
                                          + sessions_df['adas_Q7']
                                          + sessions_df['adas_Q8']
                                          + sessions_df['adas_Q9'])  # / 45
            sessions_df['adas_language'] = (sessions_df['adas_Q2']
                                            + sessions_df['adas_Q5']
                                            + sessions_df['adas_Q10']
                                            + sessions_df['adas_Q11']
                                            + sessions_df['adas_Q12'])  # / 25
            sessions_df['adas_praxis'] = (sessions_df['adas_Q3']
                                          + sessions_df['adas_Q6'])  # / 10
            sessions_df['adas_concentration'] = sessions_df['adas_Q13']  # / 5

            list_diagnosis_nan = np.where(pd.isnull(sessions_df['diagnosis']))
            diagnosis_change = {1: 'CN', 2: 'MCI', 3: 'AD'}

            for j in list_diagnosis_nan[0]:
                if not is_nan(sessions_df['adni_diagnosis_change'].iloc[j]) and int(sessions_df['adni_diagnosis_change'].iloc[j]) < 4:
                    sessions_df['diagnosis'].iloc[j] = diagnosis_change[int(sessions_df['adni_diagnosis_change'].iloc[j])]

            sessions_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False, encoding='utf-8')


def update_sessions_dict(sessions_dict, subj_bids, visit_id, field_value, bids_field_name):
    """
    Update the sessions dictionary for the bids subject specified by subj_bids
    created by the method create_adni_sessions_dict

    Args:
        sessions_dict: the session_dict created by the method
        create_adni_sessions_dict subj_bids: bids is of the subject for which
        the information need to be updated visit_id: session name (ex m54)
        field_value: value of the field extracted from the original clinical
        data bids_field_name: BIDS name of the field to update (Ex: diagnosis
        or examination date)

    Returns:
        session_dict: the dictonary updated

    """

    if visit_id == 'sc' or visit_id == 'uns1':
        return sessions_dict

    visit_id = viscode_to_session(visit_id)

    if bids_field_name == 'diagnosis':
        field_value = convert_diagnosis_code(field_value)

    # If the dictionary already contain the subject add or update information
    # regarding a specific session, otherwise create the entry
    if subj_bids in sessions_dict:
        sess_available = sessions_dict[subj_bids].keys()

        if visit_id in sess_available:
            # If a value is already contained, update it only if the previous value is nan
            if bids_field_name in sessions_dict[subj_bids][visit_id]:

                if is_nan(sessions_dict[subj_bids][visit_id][bids_field_name]):
                    sessions_dict[subj_bids][visit_id].update(
                        {bids_field_name: field_value})
            else:
                sessions_dict[subj_bids][visit_id].update(
                    {bids_field_name: field_value})
        else:
            sessions_dict[subj_bids].update({visit_id: {'session_id': 'ses-' + visit_id,
                                                        bids_field_name: field_value}})
    else:
        sessions_dict.update({subj_bids: {visit_id: {'session_id': 'ses-' + visit_id,
                                                     bids_field_name: field_value}}})

    return sessions_dict


def create_adni_sessions_dict(bids_ids, clinic_specs_path, clinical_data_dir, bids_subjs_paths):
    """
    Extract all the data required for the sessions files and organize them in a
    dictionary

    Args: bids_ids: clinic_specs_path: path to the specifications for
    converting the clinical data clinical_data_dir: path to the clinical data
    folder bids_subjs_paths: a list with the path to all the BIDS subjects

    Returns:

    """
    import pandas as pd
    from os import path
    import clinica.iotools.bids_utils as bids
    from colorama import Fore
    from clinica.utils.stream import cprint

    # Load data
    sessions = pd.read_excel(clinic_specs_path, sheetname='sessions.tsv')
    sessions_fields = sessions['ADNI']
    field_location = sessions['ADNI location']
    sessions_fields_bids = sessions['BIDS CLINICA']
    fields_dataset = []
    previous_location = ''
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    sessions_df = pd.DataFrame(columns=fields_bids)

    for i in range(0, len(field_location)):
        # If the i-th field is available
        if (not pd.isnull(field_location[i])) and path.exists(path.join(clinical_data_dir, field_location[i].split('/')[0])):
            # Load the file
            tmp = field_location[i].split('/')
            location = tmp[0]
            if location != previous_location:
                previous_location = location
                file_to_read_path = path.join(clinical_data_dir, location)
                cprint('\tReading clinical data file : ' + location)
                file_to_read = pd.read_csv(file_to_read_path, dtype=str)

                for r in range(0, len(file_to_read.values)):
                    row = file_to_read.iloc[r]

                    # Depending of the file that needs to be open, identify and
                    # do needed preprocessing on the column that contains the
                    # subjects ids
                    if location == 'ADNIMERGE.csv':
                        id_ref = 'PTID'
                        # What was the meaning of this line ?
                        # subj_id = row[id_ref.decode('utf-8')]
                        subj_id = row[id_ref]
                        subj_id = bids.remove_space_and_symbols(subj_id)
                    else:
                        id_ref = 'RID'
                        rid = str(row[id_ref])

                        # Fill the rid with the needed number of zero
                        if 4 - len(rid) > 0:
                            zeros_to_add = 4 - len(rid)
                            subj_id = '0' * zeros_to_add + rid
                        else:
                            subj_id = rid

                    # Extract the BIDS subject id related with the original
                    # subject id
                    subj_bids = [s for s in bids_ids if subj_id in s]

                    if len(subj_bids) == 0:
                        pass
                    elif len(subj_bids) > 1:
                        raise('Error: multiple subjects found for the same RID')
                    else:
                        subj_bids = subj_bids[0]
                        for i in range(0, len(sessions_fields)):
                            # If the i-th field is available
                            if not pd.isnull(sessions_fields[i]):
                                # Extract only the fields related to the current file opened
                                if location in field_location[i]:
                                    if location == 'ADAS_ADNIGO2.csv' or location == 'DXSUM_PDXCONV_ADNIALL.csv' or location == 'CDR.csv' or location == 'NEUROBAT.csv':
                                        if type(row['VISCODE2']) == float:
                                            continue
                                        visit_id = row['VISCODE2']
                                        # Convert sc to bl
                                        if visit_id == 'sc':
                                            visit_id = 'bl'
                                    else:
                                        visit_id = row['VISCODE']
                                    try:
                                        field_value = row[sessions_fields[i]]
                                        bids_field_name = sessions_fields_bids[i]
                                        sessions_dict = update_sessions_dict(sessions_dict, subj_bids, visit_id, field_value, bids_field_name)
                                    except KeyError:
                                        pass
                                        # cprint('Field value ' + Fore.RED + sessions_fields[i] + Fore.RESET + ' could not be added to sessions.tsv')
            else:
                continue

    # Write the sessions dictionary created in several tsv files
    write_adni_sessions_tsv(sessions_dict, fields_bids, bids_subjs_paths)


def create_adni_scans_files(clinic_specs_path, bids_subjs_paths, bids_ids):
    """

    Create scans.tsv files for ADNI

    Args:
        clinic_specs_path: path to the clinical file
        bids_subjs_paths: list of bids subject paths
        bids_ids: list of bids ids

    Returns:

    """
    from glob import glob
    import os
    import pandas as pd
    from os import path
    from os.path import normpath

    scans_dict = {}

    for bids_id in bids_ids:
        scans_dict.update({bids_id: {'T1/DWI/fMRI': {}, 'FDG': {}}})

    scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
    scans_fields_db = scans_specs['ADNI']
    scans_fields_bids = scans_specs['BIDS CLINICA']
    scans_fields_mod = scans_specs['Modalities related']
    fields_bids = ['filename']

    for i in range(0, len(scans_fields_db)):
        if not pd.isnull(scans_fields_db[i]):
            fields_bids.append(scans_fields_bids[i])

    scans_df = pd.DataFrame(columns=(fields_bids))

    for bids_subj_path in bids_subjs_paths:
        # Create the file
        bids_id = os.path.basename(normpath(bids_subj_path))

        sessions_paths = glob(path.join(bids_subj_path, 'ses-*'))
        for session_path in sessions_paths:
            session_name = session_path.split(os.sep)[-1]
            tsv_name = bids_id + '_' + session_name + "_scans.tsv"

            # If the file already exists, remove it
            if os.path.exists(path.join(session_path, tsv_name)):
                os.remove(path.join(session_path, tsv_name))

            scans_tsv = open(path.join(session_path, tsv_name), 'a')
            scans_df.to_csv(scans_tsv, sep='\t', index=False, encoding='utf-8')

            # Extract modalities available for each subject
            mod_available = glob(path.join(session_path, '*'))
            for mod in mod_available:
                mod_name = os.path.basename(mod)
                files = glob(path.join(mod, '*'))
                for file in files:
                    file_name = os.path.basename(file)
                    if mod == "anat" or mod == "dwi" or mod == "func":
                        type_mod = 'T1/DWI/fMRI'
                    else:
                        type_mod = 'FDG'

                    scans_df['filename'] = pd.Series(path.join(mod_name, file_name))
                    scans_df.to_csv(scans_tsv, header=False, sep='\t', index=False, encoding='utf-8')

            scans_df = pd.DataFrame(columns=(fields_bids))


def t1_pet_paths_to_bids(images, bids_dir, modality, mod_to_update=False):
    from multiprocessing import Pool, cpu_count, Value
    from clinica.iotools.converters.adni_to_bids import adni_utils
    from functools import partial

    if modality.lower() not in ['t1', 'av45', 'fdg']:
        # This should never be reached
        raise RuntimeError(modality.lower()
                           + ' is not supported for conversion in paths_to_bids')

    counter = None

    def init(args):
        ''' store the counter for later use '''
        global counter
        counter = args

    counter = Value('i', 0)
    total = images.shape[0]
    # Reshape inputs to give it as a list to the workers
    images_list = []
    for i in range(total):
        images_list.append(images.iloc[i])

    # Handle multiargument for create_file functions
    partial_create_file = partial(adni_utils.create_file,
                                  modality=modality,
                                  total=total,
                                  bids_dir=bids_dir,
                                  mod_to_update=mod_to_update)
    # intializer are used with the counter variable to keep track of how many
    # files have been processed
    poolrunner = Pool(max(cpu_count() - 1, 1), initializer=init, initargs=(counter,))
    output_file_treated = poolrunner.map(partial_create_file, images_list)
    del counter
    return output_file_treated


def create_file(image, modality, total, bids_dir, mod_to_update):
    import subprocess
    from colorama import Fore
    from clinica.utils.stream import cprint
    from clinica.iotools.converters.adni_to_bids import adni_utils
    from numpy import nan
    import os
    from os import path
    from glob import glob

    modality_specific = {'t1': {'output_path': 'anat',
                                'output_filename': '_T1w'},
                         'av45': {'output_path': 'pet',
                                  'output_filename': '_task-rest_acq-av45_pet'},
                         'fdg': {'output_path': 'pet',
                                 'output_filename': '_task-rest_acq-fdg_pet'}
                         }

    global counter

    with counter.get_lock():
        counter.value += 1

    subject = image.Subject_ID
    if image.Path == '':
        cprint(Fore.RED + '[' + modality.upper() + '] No path specified for '
               + image.Subject_ID + ' in session '
               + image.VISCODE + ' ' + str(counter.value)
               + ' / ' + str(total) + Fore.RESET)
        return nan
    cprint('[' + modality.upper() + '] Processing subject ' + str(subject)
           + ' - session ' + image.VISCODE + ', ' + str(counter.value)
           + ' / ' + str(total))
    session = adni_utils.viscode_to_session(image.VISCODE)
    image_path = image.Path
    image_id = image.Image_ID
    # If the original image is a DICOM, check if contains two DICOM
    # inside the same folder
    if image.Is_Dicom:
        image_path = adni_utils.check_two_dcm_folder(image_path,
                                                     bids_dir,
                                                     image_id)
    bids_subj = subject.replace('_', '')
    output_path = path.join(bids_dir, 'sub-ADNI' + bids_subj, 'ses-'
                            + session, modality_specific[modality]['output_path'])
    output_filename = 'sub-ADNI' + bids_subj + '_ses-' + session + modality_specific[modality]['output_filename']

    # If updated mode is selected, check if an old image is existing and remove it
    existing_im = glob(path.join(output_path, output_filename + '*'))
    if mod_to_update and len(existing_im) > 0:
        print('Removing the old image...')
        for im in existing_im:
            os.remove(im)

    try:
        os.makedirs(output_path)
    except OSError:
        # Folder already created with previous instance
        pass

    if image.Is_Dicom:
        command = 'dcm2niix -b n -z n -o ' + output_path + ' -f ' + output_filename + ' ' + image_path
        subprocess.run(command,
                       shell=True,
                       stderr=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL)
        nifti_file = path.join(output_path, output_filename + '.nii')
        output_image = nifti_file + '.gz'

        # Check if conversion worked (output file exists?)
        if not path.isfile(nifti_file):
            cprint('\tConversion with dcm2niix failed, trying with dcm2nii')
            command = 'dcm2nii -a n -d n -e n -i y -g n -p n -m n -r n -x n -o ' + output_path + ' ' + image_path
            subprocess.run(command,
                           shell=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            nifti_file = path.join(output_path,
                                   subject.replace('_', '') + '.nii')
            output_image = path.join(output_path,
                                     output_filename + '.nii.gz')

            if not path.isfile(nifti_file):
                cprint('DICOM to NIFTI conversion error for ' + image_path)
                return nan

        adni_utils.center_nifti_origin(nifti_file, output_image)
        os.remove(nifti_file)

    else:
        output_image = path.join(output_path, output_filename + '.nii.gz')
        adni_utils.center_nifti_origin(image_path, output_image)

    # Check if there is still the folder tmp_dcm_folder and remove it
    adni_utils.remove_tmp_dmc_folder(bids_dir, image_id)
    return output_image
