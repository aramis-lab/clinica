"""
Function used by BIDS converters.
"""

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"


# @ToDo:test this function
def create_participants_df(input_path,out_path, study_name, clinical_spec_path, bids_ids, delete_non_bids_info=True):
    """


    :param input_path: path to the original dataset
    :param out_path: path to the bids folder
    :param study_name: name of the study (Ex. ADNI)
    :param clinical_spec_path:
    :param bids_ids:
    :param delete_non_bids_info:
    :return: a pandas dataframe that contains the participants data
    """
    import pandas as pd
    import os
    from os import path
    import logging
    import numpy as np

    fields_bids = ['participant_id']
    fields_dataset = []
    prev_location = ''
    index_to_drop=[]

    location_name = study_name + ' location'

    participants_specs = pd.read_excel(clinical_spec_path , sheetname='participant.tsv')
    participant_fields_db = participants_specs[study_name]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs['BIDS CLINICA']

    # Extract the list of the available fields for the dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])
            fields_dataset.append(participant_fields_db[i])

    # Init the dataframe that will be saved in the file participant.tsv
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
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(input_path, 'clinicalData', location)

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
                if participant_fields_bids[i] == 'alternative_id_1' and\
                        (file_to_read[participant_fields_db[i]].dtype == np.float64 or file_to_read[participant_fields_db[i]].dtype == np.int64) :
                    if not pd.isnull(file_to_read.get_value(j, participant_fields_db[i])):
                        value_to_append = str(file_to_read.get_value(j, participant_fields_db[i])).rstrip('.0')
                    else:
                        value_to_append = np.NaN
                else:
                    value_to_append = file_to_read.get_value(j, participant_fields_db[i])
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    if study_name == 'ADNI':
        # ADNImerge contains one row for each visits so there are duplicates
        participant_df = participant_df.drop_duplicates(subset=['alternative_id_1'], keep='first')
    participant_df.reset_index(inplace=True, drop=True)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        value = remove_space_and_symbols(participant_df['alternative_id_1'][i])
        bids_id = [s for s in bids_ids if value in s]
        if len(bids_id) == 0:
            print "Subject " + value + " not found in the BIDS converted version of the dataset."
            logging.error("Subject " + value + " not found in the BIDS converted version of the dataset.")
            index_to_drop.append(i)
        else:
            participant_df['participant_id'][i] = bids_id[0]

    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info == True:
        participant_df = participant_df.drop(index_to_drop)

    return participant_df

def create_sessions_dict(input_path, study_name, clinical_spec_path, bids_ids, name_column_ids, subj_to_remove = []):
    """
    Extract the information regarding the sessions and store them in a dictionary (session M0 only)

    :param input_path:
    :param study_name: name of the study (Ex: ADNI)
    :param clinical_spec_path:
    :param bids_ids:
    :param name_column_ids:
    :param subj_to_remove:
    :return:
    """
    import pandas as pd
    from os import path
    import numpy as np

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
            sheet = tmp[1]
            file_to_read_path = path.join(input_path, 'clinicalData', location)
            file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]
                # Extracts the subject ids columns from the dataframe
                subj_id = row[name_column_ids.decode('utf-8')]
                if subj_id.dtype == np.int64:
                    subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)

                # Extract the corresponding BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if len(subj_bids) == 0:
                    # If the subject is not an excluded one
                    if not subj_id in subj_to_remove:
                        print sessions_fields[i] + ' for ' + subj_id + ' not found in the BIDS converted.'
                else:
                    subj_bids = subj_bids[0]
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]
                    if sessions_dict.has_key(subj_bids):
                        (sessions_dict[subj_bids]['M0']).update({sessions_fields_bids[i]: row[sessions_fields[i]]})
                    else:
                        sessions_dict.update({subj_bids: {
                            'M0': {'session_id': 'ses-M0', sessions_fields_bids[i]: row[sessions_fields[i]]}}})

    return sessions_dict

def write_sessions_tsv(out_path, bids_paths, sessions_dict, fields_bids, sessions_list = 'M0'):
    """

    :param out_path:
    :param bids_paths:
    :param sessions_dict:
    :param fields_bids:
    :param sessions_list:
    :return:
    """
    import os
    import pandas as pd
    from os import path

    for sp in bids_paths:
        sp = sp[:-1]
        bids_id = sp.split(os.sep)[-1]
        sessions_df = pd.DataFrame(columns=fields_bids)
        if sessions_dict.has_key(bids_id):
            session_df = pd.DataFrame(sessions_dict[bids_id]['M0'], index=['i', ])
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)
        else:
            print "No session data available for " + sp
            session_df = pd.DataFrame(columns=['session_id'])
            session_df['session_id'] = pd.Series('M0')
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)


def create_empty_sessions_tsv(): pass


def create_empy_scans_tsv(): pass


def get_bids_subjs_list(bids_path):
    """
    Given a BIDS compliant dataset, returns the list of all the subjects available

    :param bids_path: path to the BIDS folder
    :return:
    """
    import os
    from os import path

    return [d for d in os.listdir(bids_path) if os.path.isdir(path.join(bids_path, d))]


def get_bids_subjs_paths(bids_path):
    """
    Given a BIDS compliant dataset, returns the list of all paths to the subjects folders

    :param bids_path: path to the BIDS folder
    :return:
    """
    import os
    from os import path

    return [path.join(bids_path,d) for d in os.listdir(bids_path) if os.path.isdir(path.join(bids_path, d))]


def compute_new_subjects(original_ids, bids_ids):
    """
    Check for new subject to convert.

    This function checks for news subjects to convert to the BIDS version i.e. subjects contained in the unorganised
    version that are not available in the bids version.

    :param original_ids: list of all the ids of the unorganized folder.
    :param bids_ids: list of all the BIDS ids contained inside the bids converted version of the dataset
    :return: a list containing the original_ids of the subjects that are not available in the bids converted version
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
    Compress .nii file
    :param path:
    :return:
    '''

    from os import path
    import os
    import gzip

    f_in = open(file_path)
    f_out = gzip.open(path.join(file_path + '.gz'), 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()

    # Remove the original file
    os.remove(file_path)


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

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    #copy(t1_path, path.join(output_path, t1_bids_name + get_bids_suff('T1') + '.nii.gz'))
    file_ext = get_ext(t1_path)
    copy(t1_path, path.join(output_path, t1_bids_name + get_bids_suff('T1') + file_ext))
    # If  the original image is not compress, compress it
    if file_ext == '.nii':
        compress_nii(path.join(output_path, t1_bids_name + get_bids_suff('T1') + file_ext))


def convert_fieldmap(folder_input, folder_output, name, fixed_file=[False,False]):
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
        map = remove_rescan(glob(path.join(folder_input, "*MAP_*",'*.nii.gz')))
        map_ph = remove_rescan(glob(path.join(folder_input, "*MAPph_*", '*.nii.gz')))
    elif fixed_file[0] is False and fixed_file[1] is not False:
        map = remove_rescan(glob(path.join(folder_input, "*MAP_*",'*.nii.gz')))
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
                mag_list = glob(path.join(folder_output,name+get_bids_suff('Map')+'*'))

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
        if mag_missing == True and  map_ph_missing == True:
            return -1
        else:
            return 0


def convert_flair(folder_input, folder_output, name, fixed_file = False):
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
    if fixed_file!=False:
        fixed_flair_path = glob(path.join(folder_input,fixed_file,'*'))[0]
        copy(fixed_flair_path, path.join(folder_output, name + (get_bids_suff('Flair')) + '.nii.gz'))
    else:
        list_path = glob(path.join(folder_input,'*T2FLAIR*'))
        flair_lst = remove_rescan(list_path)
        if len(flair_lst) == 1:
            if not os.path.exists(folder_output):
                os.mkdir(folder_output)
            flair_path = glob(path.join(flair_lst[0], '*.nii.gz*'))[0]
            copy(flair_path, path.join(folder_output, name + (get_bids_suff('Flair')) + '.nii.gz'))
        elif len(flair_lst) == 0:
                return -1
        elif len(flair_lst)>1:
                logging.warning('Multiple FLAIR found, computation aborted.')
                raise


def convert_fmri(folder_input, folder_output, name, fixed_fmri=False):
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

    if fixed_fmri is not False:
        fmri_lst = glob(path.join(folder_input, fixed_fmri))
    else:
        fmri_lst = remove_rescan(glob(path.join(folder_input, '*fMRI*')))
    if len(fmri_lst) > 0:
        if not os.path.exists(folder_output):
            os.mkdir(folder_output)
        fmri_file_path = glob(path.join(fmri_lst[0], '*.nii*'))[0]
        file_ext = get_ext(fmri_file_path)
        copy(fmri_file_path, path.join(folder_output, name + '_task-rest' + get_bids_suff('fMRI') + file_ext))
        if file_ext == '.nii':
            logging.warning('Non compressed file found: '+fmri_file_path)
            compress_nii(path.join(folder_output, name + '_task-rest' + get_bids_suff('fMRI') + file_ext))
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
                img.append(glob(path.join(folder,'*.nii*'))[0])
                bval.append(glob(path.join(folder,'*.bval'))[0])
                bvec.append(glob(path.join(folder,'*.bvec'))[0])
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
                line_no +=1
            for i in range (0, len(lines_out_bval)):
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
