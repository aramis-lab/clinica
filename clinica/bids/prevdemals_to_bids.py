# coding=utf-8
"""
Covert the PREVDEMALS dataset into the BIDS specification.

ToFix:
    The conversion of the Clinical data doesn't work if the subjects list contained in each file to open is not exactly the same

@author: Sabrina Fontanella
"""

from os import path
from glob import glob
import os
from os.path import normpath
import logging
import bids_utils as bids
import pandas as pd
import pkg_resources as pkg
import sys
import re


def remove_space_and_symbols(data):
    '''
    Remove spaces and  - _ from a list (or a single) of strings
    :param data: list of strings or a single string to clean
    :return:
        data: list of strings or a string without space and symbols _ and -
    '''
    if type(data) is list:
        for i in range(0, len(data)):
            data[i] = re.sub('[-_ ]', '', data[i])
    else:
        data = re.sub('[-_ ]', '', data)

    return data


def remove_chars(data):
    '''
    Remove characters from a list of strings.

    :param data:
    :return:
    '''
    non_decimal = re.compile((r'[^\d.]+'))
    for i in range (0, len(data)):
        if type(data[i]) is unicode:
            to_clean = data[i]
            data[i] = non_decimal.sub('', to_clean)

    return data


def convert_clinical(input_path, out_path, bids_ids):
    pd.options.mode.chained_assignment = None
    reload(sys)
    sys.setdefaultencoding('utf-8')
    subj_to_remove = ['001-0029-CF', '001-0031-SB', '001-0032-SC', '003-0001-D-B', '004-0003-F-B']
    fields_bids_aval = ['participant_id']
    fields_aval = []
    prev_location = ""
    prev_sheet = ""
    prev_field = ""
    clinic_specs_path = pkg.resource_filename('clinica', 'bids/data/clinical_specifications.xlsx')

    genfi_subjs = glob((path.join(out_path, 'GENFI', 'sub-*')))
    icm_subjs =  glob((path.join(out_path, 'ICM', 'sub-*' )))

    # Extract all the bids subjects pats
    subjs_bids_path = glob(path.join(out_path,"GENFI", '*'))
    subjs_bids_path.append(glob(path.join(out_path,"ICM", '*')))

    # -- Creation of participants.tsv --
    logging.info("-- Creation of participants file --")
    print("-- Creation of participants file --")
    participants_specs = pd.read_excel(clinic_specs_path, sheetname='participant.tsv')

    # Extracts information regarding INSIGHT dataset
    participant_fields_db = participants_specs['PREVDEMALS']
    field_location = participants_specs['PREVDEMALS location']
    participant_fields_bids = participants_specs['BIDS CLINICA']
    field_type_bids = participants_specs['BIDS type']

    # Extract the list of the available fields for the CAPP dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids_aval.append(participant_fields_bids[i])
            fields_aval.append(participant_fields_db[i])

    # Init the dataframe that will be saved in the file participant.tsv
    participant_df = pd.DataFrame(columns=fields_bids_aval)

    for i in range(0, len(participant_fields_db)):
            # If a field not empty is found
            if not pd.isnull(participant_fields_db[i]):
                # Extract the file location of the field and read the value from the file
                tmp = field_location[i].split('/')
                location = tmp[0]
                sheet = tmp[1]
                # Check if the file to open for a certain field it's the same of the previous field
                if location == prev_location and sheet == prev_sheet:
                    pass
                else:
                    file_to_read_path = path.join(input_path, 'clinicalData', location)
                    if sheet in '1ere extract datas Paris_Mariem':
                        file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet, header=[1])
                    else:
                        file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
                    prev_location = location
                    prev_sheet = sheet

                field_col_values = []

                # It's a special field that need to be splitted
                if participant_fields_db[i] == prev_field:
                    pass
                else:
                    # For each field in fields_dataset extract all the column values
                    prev_field = participant_fields_db[i]
                    if participant_fields_db[i] in 'Transmitting parent phenotype':
                        # O = parent_stransmitting_fdt
                        # 1 = parent_transmitting_als
                        # 2 = relatives
                        field_col_values = [[], [], []]
                        for j in range(0, len(file_to_read)):
                            value_to_split = file_to_read.get_value(j, participant_fields_db[i])
                            # Is float when is nan
                            if type(value_to_split) != float:
                                if 'FTD' in value_to_split:
                                    field_col_values[0].append('Y')
                                else:
                                    field_col_values[0].append('N')
                                if "ALS" in value_to_split:
                                    field_col_values[1].append('Y')
                                else:
                                    if "SLA" in value_to_split:
                                        field_col_values[1].append('Y')
                                    else:
                                        field_col_values[1].append('N')
                                if 'mother' in value_to_split:
                                    field_col_values[2].append('mother')
                                else:
                                    field_col_values[2].append('father')
                            else:
                                field_col_values[0].append('')
                                field_col_values[1].append('')
                                field_col_values[2].append('')
                        participant_df['parent_transmitting_fdt'] = pd.Series(field_col_values[0])
                        participant_df['parent_transmitting_als'] = pd.Series(field_col_values[1])
                        participant_df['parent_transmitting'] = pd.Series(field_col_values[2])
                    elif participant_fields_db[i] == 'Family phenotype':
                        # familiar_phenotype_als
                        # familiar_phenotype_ftd
                        field_col_values = [[], []]

                        for j in range(0, len(file_to_read)):
                            value_to_split = file_to_read.get_value(j, participant_fields_db[i])
                            # Is float when is nan
                            if type(value_to_split) != float:
                                if 'FTD' in value_to_split:
                                    field_col_values[0].append('Y')
                                else:
                                    field_col_values[0].append('N')
                                if 'ALS' in value_to_split:
                                    field_col_values[1].append('Y')
                                else:
                                    field_col_values[1].append('N')
                            else:
                                field_col_values[0].append('')
                                field_col_values[1].append('')
                        participant_df['family_phenotype_ftd'] = pd.Series(field_col_values[0])
                        participant_df['family_phenotype_als'] = pd.Series(field_col_values[1])
                    elif participant_fields_db[i] == 'Genetic status (m=mutation, nm=no mutation)':
                        field_col_values = []
                        for j in range(0, len(file_to_read)):
                            genetic_status = file_to_read.get_value(j, participant_fields_db[i])
                            if type(genetic_status) != float:
                                if 'non' in genetic_status:
                                    field_col_values.append('N')
                                else:
                                    field_col_values.append('Y')
                            else:
                                field_col_values.append('')
                        participant_df['prevdemals_genetic_status_mutated'] = pd.Series(field_col_values)
                    elif participant_fields_db[i] == 'Status (A=patients, R=relatives)':
                        field_col_values = []
                        for j in range(0, len(file_to_read)):
                            prevdemals_status = file_to_read.get_value(j, participant_fields_db[i])
                            if type(prevdemals_status) != float:
                                if 'A' in prevdemals_status:
                                    field_col_values.append('P')
                                elif 'R' in prevdemals_status:
                                    field_col_values.append('R')
                            else:
                                field_col_values.append('')
                        participant_df['prevdemals_status'] = pd.Series(field_col_values)
                    else:
                        for j in range(0, len(file_to_read)):
                            field_col_values.append(file_to_read.get_value(j, participant_fields_db[i]))

                        if field_type_bids[i] == 'int' or field_type_bids[i] == 'float':
                            field_col_values = remove_chars(field_col_values)
                        # if participant_fields_bids[i] == 'aoo_mean_family' or participant_fields_bids[i] == 'expected_years_aoo_bl':
                        #     field_col_values = [round(float(num)) for num in field_col_values ]


                        # Add the extracted column to the participant_df
                        participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    # Adding participant_id column with BIDS ids and remove blank row
    index_to_remove = []
    for i in range(0, len(participant_df)):
        value = participant_df['alternative_id_1'][i]

        if pd.isnull(value):
            index_to_remove.append(i)
        else:
            # Removing all -
            value = value.replace('-', '')
            bids_id = [s for s in bids_ids if value in s]
            if len(bids_id) == 0:
                if not value in subj_to_remove:
                    print "Subject " + value + " not found in the BIDS converted version of PREVDEMALS."
                    logging.error("Subject " + value + " not found in the BIDS converted version of PREVDEMALS.")
                participant_df['participant_id'][i] = ''
            else:
                participant_df['participant_id'][i] = bids_id[0]

    # Removes all subjects discarded from the study
    for s in subj_to_remove:
        index_to_remove.append(participant_df[ participant_df['alternative_id_1'] == s ].index.tolist()[0])

    participant_df = participant_df.drop(index_to_remove)

    # Split the participant.tsv files in two files: one containing the GENFI subjects and the other the ICM subjects
    genfi_participant_df = pd.DataFrame(columns=fields_bids_aval)
    icm_participant_df = pd.DataFrame(columns=fields_bids_aval)
    for gsp in genfi_subjs:
        genfi_sub_name = gsp.split(os.sep)[-1]
        index_to_extract = participant_df[participant_df['participant_id'] == genfi_sub_name].index.tolist()
        if len(index_to_extract) == 0:
            print 'Subject '+genfi_sub_name+' not found in the list of subject available'
        else:
            genfi_participant_df = genfi_participant_df.append(participant_df.loc[index_to_extract[0]])

    for isp in icm_subjs:
        icm_sub_name = isp.split(os.sep)[-1]
        index_to_extract = participant_df[participant_df['participant_id'] == icm_sub_name].index.tolist()
        if len(index_to_extract) == 0:
            print 'Subject '+icm_sub_name+' not found in the list of subject available'
        else:
            icm_participant_df = icm_participant_df.append(participant_df.loc[index_to_extract[0]])


    icm_participant_df.to_csv(path.join(out_path,'ICM', 'participants.tsv'), sep='\t', index=False)
    genfi_participant_df.to_csv(path.join(out_path, 'GENFI', 'participants.tsv'), sep='\t', index=False)

    logging.info("Participants file created.\n")

    # -- Creation of sessions.tsv --
    logging.info("--Creation of sessions files. --")
    print("\nCreation of sessions files...")
    # Load data
    sessions = pd.read_excel(clinic_specs_path, sheetname='sessions.tsv')
    sessions_fields = sessions['PREVDEMALS']
    field_location = sessions['PREVDEMALS location']
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
            file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet, header=[1])

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]
                # Extracts the subject ids columns from the dataframe
                id_ref = 'Référence patient'
                subj_id = row[id_ref.decode('utf-8')]
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)
                # Extracts the correspondant BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if len(subj_bids) == 0:
                    # If the subject is not an exluded one
                    if not subj_id in subj_to_remove:
                        print sessions_fields[i]+' for '+subj_id+' not found in the BIDS converted.'
                else:
                    subj_bids = subj_bids[0]
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]
                    if sessions_dict.has_key(subj_bids):
                        (sessions_dict[subj_bids]['M00']).update({sessions_fields_bids[i]: row[sessions_fields[i]]})
                    else:
                        sessions_dict.update({subj_bids: {'M00': {'session_id': 'M00', sessions_fields_bids[i]: row[sessions_fields[i]]}}})

    subjs_bids_path = glob(path.join(out_path, 'GENFI','*/'))
    icm_subj = glob(path.join(out_path, 'ICM','*/'))
    for s in icm_subj:
        subjs_bids_path.append(s)
    for sp in subjs_bids_path:
        sp = sp[:-1]
        bids_id = sp.split(os.sep)[-1]
        sessions_df = pd.DataFrame(columns=fields_bids)
        if sessions_dict.has_key(bids_id):
            session_df = pd.DataFrame(sessions_dict[bids_id]['M00'], index=['i', ])
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index = False)
        else:
            logging.warning("No session data available for "+ sp)
            print "No session data available for " + sp
            session_df =  pd.DataFrame(columns=['session_id'])
            session_df['session_id'] = pd.Series('M00')
            session_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)


    logging.info("Sessions files created.")
    print("Sessions files created.")

    # -- Creation of *_scans.tsv --
    print 'Creation of scans files...'
    scans_dict = {}

    for bids_id in bids_ids:
        scans_dict.update({bids_id: {'T1/DWI/fMRI': {}, 'FDG': {}}})

    scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
    scans_fields_db = scans_specs['PREVDEMALS']
    field_location = scans_specs['PREVDEMALS location']
    scans_fields_bids = scans_specs['BIDS CLINICA']
    scans_fields_mod = scans_specs['Modalities related']
    fields_bids = ['filename']

    for i in range(0, len(scans_fields_db)):
        if not pd.isnull(scans_fields_db[i]):
            fields_bids.append(scans_fields_bids[i])

    scans_df = pd.DataFrame(columns=(fields_bids))

    for bids_subj_path in subjs_bids_path:
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
            scans_df.to_csv(scans_tsv, sep='\t', index=False)

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
                    scans_df.to_csv(scans_tsv, header=False, sep='\t', index=False)

            scans_df = pd.DataFrame(columns=(fields_bids))

    print '-- Scans files created for each subject. --'


def convert(source_dir, dest_dir, param=''):
    """
    Convert the PREVDEMALS dataset into the BIDS standard.

    Args:
        source_dir: directory of the input dataset
        dest_dir: output directory
    """

    t1_priority = ['3DT1_PN_noDIS', '3DT1_S']
    projects = {
        'GENFI': path.join(source_dir, 'PREV_DEMALS_GENFI', 'convertData', 'study'),
        'ICM': path.join(source_dir, 'PREV_DEMALS_ICM', 'convertData', 'study')
    }
    ses_available = ['M0', 'M24']
    #mmt = MissingModsTracker(['M0', 'M24'])
    if param!='-c':
        os.mkdir(dest_dir)
        os.mkdir(path.join(dest_dir, 'GENFI'))
        os.mkdir(path.join(dest_dir, 'ICM'))
    # summary_file_icm = open(path.join(dest_dir, 'ICM', 'conversion_summary.txt'), 'w')
    # summary_file_genfi = open(path.join(dest_dir, 'GENFI', 'conversion_summary.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'conversion_modalities.log'), format='%(asctime)s %(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)
    print "*******************************"
    print "PREVDEMALS to BIDS converter"
    print "*******************************"
    # Convert all the files contained in the two project folder PREV_DEMALS_GENFI and PREV_DEMALS_ICM
    if param != '-c':
        for proj in projects:
            cities_folder = glob(path.join(projects[proj], '*'))
            print

            pda_spath = []
            pda_ids = []
            bids_ids = []
            dest_dir_proj = path.join(dest_dir, proj)
            for cf in cities_folder:
                for subj_path in glob(path.join(cf, "*")):
                    subj_id = subj_path.split(os.sep)[-1]
                    pda_spath.append(subj_path)
                    pda_ids.append(subj_id)
                    bids_ids.append('sub-PREVDEMALS' + subj_id)
                    # Create the subject folder in the BIDS converted dataset
                    os.mkdir(path.join(dest_dir_proj, bids_ids[-1]))


            # For each subject extract the list of files and convert them into BIDS specification
            for subj_path in pda_spath:
                print "Converting:", subj_path
                logging.info('Converting:'+subj_path)
                # Extract the session(s) available
                sessions = glob(path.join(subj_path, "*"))
                # For each sub-session
                for ses_path in sessions:
                    ses = ses_path.split(os.sep)[-1]
                    if 'rescan' in ses:
                        logging.warning('Rescan of a session found: '+ses+'. Ignored.')
                        continue

                    # Extracting the index of the subject that is the same for the PREVDEMALS list and BIDS list
                    subj_index = pda_spath.index(subj_path)
                    os.mkdir(path.join(dest_dir_proj, bids_ids[subj_index], 'ses-' + ses))
                    ses_dir_bids = path.join(dest_dir_proj, bids_ids[subj_index], 'ses-' + ses)
                    bids_file_name = bids_ids[subj_index] + '_ses-' + ses
                    session_path = path.join(subj_path, ses)
                    mods_folder_path = path.join(session_path, 'NIFTI')


                    # Merge all the valid DTI folders for the same subject
                    if 'ICM' in projects[proj]:
                        logging.info('DTI ignored for PREV_DEMALS_ICM.')
                    else:
                        out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                        if out is not None:  # Missing or incomplete DTI
                            if out == -1:
                                logging.info('No DTI found for ' + mods_folder_path)
                                #mmt.add_missing_mod('DTI', ses)
                            else:
                                for e in out:
                                    logging.warning('.bvec or .bval not found for DTI folder ' + e + ' Skipped')
                                    #mmt.add_incomplete_mod('DTI', ses)

                    # Decide the best T1 to take and convert the file format into the BIDS standard
                    t1_selected = bids.choose_correction(mods_folder_path, t1_priority, 'T1')
                    if t1_selected != -1 and t1_selected != 0:
                        t1_sel_path = glob(path.join(mods_folder_path, '*' + t1_selected + '*', '*.nii*'))[0]
                        bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                    else:
                        if t1_selected == -1:
                            logging.info('No T1 found for ' + mods_folder_path)
                            #mmt.add_missing_mod('T1', ses)
                        else:
                            logging.warning('None of the desiderd T1 corrections is available for ' + mods_folder_path)

                    # Extract and convert the T2 FLAIR modality if is available
                    out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                    if out == -1:
                        logging.warning('No FLAIR found for '+ mods_folder_path)
                        #mmt.add_missing_mod('FLAIR', ses)

                    logging.info("Conversion for the subject terminated.\n")

            # Printing the statistics about missing modalities into a file
            # if proj == 'ICM':
            #     file_to_write = summary_file_icm
            # else:
            #     file_to_write = summary_file_genfi
            # print_statistics(file_to_write, len(pda_ids), ses_available, mmt)
            # file_to_write.close()
    else:
        print '** Conversion of clinical data only **'
        if not os.path.exists(dest_dir):
            print dest_dir, ' not found.'
            raise

    if os.path.exists(path.join(dest_dir, 'conversion_clinical.log')):
        os.remove(path.join(dest_dir, 'conversion_clinical.log'))

    logging.basicConfig(filename=path.join(dest_dir, 'conversion_clinical.log'),
                        format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG,
                        datefmt='%m/%d/%Y %I:%M')

    print 'Converting clinical data...'
    logging.info('Converting clinical data...')
    bids_ids = [d for d in os.listdir(path.join(dest_dir, 'GENFI')) if os.path.isdir(path.join(dest_dir, 'GENFI', d))]
    icm_subjs = ([d for d in os.listdir(path.join(dest_dir, 'ICM')) if os.path.isdir(path.join(dest_dir, 'ICM', d))])
    for s in icm_subjs:
        bids_ids.append(s)
    convert_clinical(source_dir, dest_dir, bids_ids)


