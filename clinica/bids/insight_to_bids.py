"""
Covert the INSIGHT dataset into the BIDS specification.

ToDo:
- Subjects missing problems and multiple t2flairs

@author: Sabrina Fontanella
"""


from os import path
from glob import glob
import os, sys
import logging
from converter_utils import MissingModsTracker, print_statistics
import bids_utils as bids
import pandas as pd


def convert_clinical(input_path, out_path, bids_ids, original_ids, orig_ids_alpha):
    '''
    Convert clinical data of INSIGHT (CATI organized) dataset (subject included only).

    Args:
        input_path:
        out_path:
        bids_ids:
        original_ids:
    :return:

    '''

    fields_bids = ['participant_id']
    fields_dataset = []
    clinic_specs_path = path.join(os.path.dirname(__file__), 'bids-utils-docs', 'clinical_specifications.xlsx')

    # -- Creation of participant.tsv --
    participants_specs = pd.read_excel(clinic_specs_path, sheetname='participant.tsv')
    # Extracts information regarding INSIGHT dataset
    participant_fields_db = participants_specs['INSIGHT']
    field_location = participants_specs['INSIGHT location']
    participant_fields_bids = participants_specs['BIDS CLINICA']

    # Extract the list of the available fields for the INSIGHT dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])
            fields_dataset.append(participant_fields_db[i])

   # Init the dataframe that will be saved in the file participant.tsv
    participant_df = pd.DataFrame(columns= fields_bids)

    for i in range(0, len(participant_fields_db)):
         # If a field not empty is found
         if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            file_to_read_path = path.join(input_path, 'clinicalData', field_location[i]+'.xls')
            file_to_read = pd.read_excel(file_to_read_path)

    field_col_values = []
    # For each field in fields_dataset extract all the column values
    for f in range(0, len(fields_dataset)):
        for i in range(0, len(file_to_read.values)):
            field_col_values.append(file_to_read.get_value(i, fields_dataset[f]))

        # Add the extracted column to the participant_df
        participant_df[fields_bids[f+1]] = pd.Series(field_col_values)
        field_col_values = []

   # Remove data from not included subjects
    index_to_delete = []
    for i in range(0, len(participant_df)):
        value = participant_df['alternative_id_1'][i]

        # Removing the first 0 and deleting the -
        value = value[1:].replace("-", "")
        # Check if the id is contained in bids_ids
        if not any(value in id for id in bids_ids):
            # Get the index of the element to delete
            # TOFix: it works but maybe there a way much more elegant!
            index_to_delete.append((participant_df[participant_df['alternative_id_1'] == participant_df['alternative_id_1'][i]].index.tolist())[0])
        else:
            bid = [s for s in bids_ids if value in s][0]
            participant_df['participant_id'][i] = bid

    # Remove all the informations regarding the subjects not included from the dataframe
    participant_df = participant_df.drop(index_to_delete)

    participant_df.to_csv(path.join(out_path,'participants.tsv'), sep='\t', index = False)


    # Creation of sub-<subject_label>_sessions.tsv
    sessions = pd.read_excel(clinic_specs_path, sheetname='sessions.tsv')
    sessions_fields_db = sessions['INSIGHT']
    field_location = sessions['INSIGHT location']
    sessions_fields_bids = sessions['BIDS CLINICA']
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields_db)):
        if not pd.isnull(sessions_fields_db[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields_db[i])

    sessions_df = pd.DataFrame(columns= fields_bids)

    for i in range(0, len(sessions_fields_db)):
        # If a field not empty is found
        if not pd.isnull(sessions_fields_db[i]):
            # Estrae il corrispettivo valore della locazione dove si trova il campo
            tmp = field_location[i].split('/')
            location = tmp[0]
            sheet = tmp[1]
            file_to_read_path = path.join(input_path, 'clinicalData', location + '.xlsx')
            file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]
                # Extract the subject id
                subj_id = row['Subject']
                # Extract the correspondant BIDS id and create the output file if doesn't exist
                subj = [s for s in orig_ids_alpha if s in subj_id]

                # If the ID is an included one
                if len(subj) != 0:
                    subj = subj[0]
                    subj_index = orig_ids_alpha.index(subj)
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields_db[i]]
                    # Create the files
                    if sessions_dict.has_key(bids_ids[subj_index]):
                        (sessions_dict[bids_ids[subj_index]]['M0']).update({sessions_fields_bids[i]: row[sessions_fields_db[i]]})
                    else:
                        sessions_dict.update({bids_ids[subj_index]: {'M0' : {'session_id': 'M0', sessions_fields_bids[i]: row[sessions_fields_db[i]]}}})

    # Create the several sessions.tsv file extracting data from the dictionary
    for bids_id in bids_ids:
        sessions_df = pd.DataFrame(columns=fields_bids)

        if sessions_dict.has_key(bids_id):
            session_df = pd.DataFrame(sessions_dict[bids_id]['M0'], index=['i', ])
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
            session_df.to_csv(path.join(out_path,bids_id, bids_id + '_sessions.tsv'), sep='\t', index = False)
        else:
            print bids_id, ' not found into clinical datas'


def converter(source_dir, dest_dir, mod=''):
    # List of the all the subjects to take in consideration for the conversion
    subjs_list = path.join(source_dir, '_INSIGHT-AD_IDs_INCLUDED__Status-07-06-2016.xlsx')

    data_dir = path.join(source_dir, 'convertData', 'study', 'Paris') # Folder where data are stored
    t1_priority = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    mmt = MissingModsTracker(['M0', 'M24'])
    ses_available = ['M0', 'M24']

    insight_ids = []
    original_insight_ids = []
    original_insight_ids_alpha = []
    memento_ids = []
    subjs_paths = []
    bids_ids = []
    subjs_included = pd.read_excel(subjs_list, sheetname='IDs INCLUDED')

    os.mkdir(dest_dir)
    summary = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'conversion'
                                                     '.log'), format='%(asctime)s %(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)

    print "*******************************"
    print "INSIGHT to BIDS converter"
    print "*******************************"


     # The last 7 rows from the excel files are not ids but other datas
    for i in range(0, len(subjs_included.index) - 7):
        original_insight_ids.append(str(subjs_included['INSIGHT IDs'][i]).replace(" ", ""))
        original_insight_ids_alpha.append(original_insight_ids[-1].replace("-", ""))
        insight_ids.append((str(subjs_included['INSIGHT IDs'][i]).replace("-", "")).replace(" ", ""))
        memento_ids.append(str(subjs_included['MEMENTO IDs'][i]).replace("-", ""))
        # The BIDS specification allows only numbers and letter inside the subject id
        bids_ids.append('sub-INSIGHT'+insight_ids[-1].replace("-", ""))
        os.mkdir(path.join(dest_dir, bids_ids[-1]))
        # In the INSIGHT cati organized the subject is is the insight id without the character -
        subj_path = glob(path.join(data_dir, '*' + insight_ids[-1] + '*'))[0]
        subjs_paths.append(subj_path)


    # For each subject extract the list of files and convert them into BIDS specification
    for subj_path in subjs_paths:
        if mod!='-c':
            print "Converting:", subj_path
        logging.info('Converting:' + subj_path)
        # Extract the session(s) available
        sessions = glob(path.join(subj_path, "*"))

        # Check if the subject has some missing session
        for s in ses_available:
            if not any(s in ses_found for ses_found in sessions ):
                mmt.incr_missing_session(s)

        # For each sub-session
        for ses_path in sessions:
            ses = ses_path.split(os.sep)[-1]

            if 'rescan' in ses:
                logging.warning('Rescan of a session found: ' + ses + 'for ' + subj_path + '. Ignored.')
                print 'Rescan of a session found: ' + ses + '. Ignored.'
                continue

            # Extracting the index of the subject
            subj_index = subjs_paths.index(subj_path)
            os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses))
            ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses)
            bids_file_name = bids_ids[subj_index] + '_ses-' + ses
            session_path = path.join(subj_path, ses)
            mods_folder_path = path.join(session_path, 'NIFTI')

            if mod != "-c":
                # Convert the fieldmap data
                out = bids.convert_fieldmap(mods_folder_path, path.join(ses_dir_bids, "fmap"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('Fieldmap', ses)
                    logging.warning("No Fieldmap found for "+mods_folder_path)
                elif out == 0:
                    logging.warning("Map or MapPh missing for "+mods_folder_path +': skipped.')
                    mmt.add_incomplete_mod('Fieldmap', ses)

                #Convert DTI data
                out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                if out != None:  # missing or incomplete DTI
                    if out == -1:
                        logging.info('No DTI found for ' + mods_folder_path)
                        mmt.add_missing_mod('DTI', ses)
                    else:
                        for e in out:
                            logging.warning('.bvec or .bval not found for DTI folder ' + e + ': skipped.')
                            mmt.add_incomplete_mod('DTI', ses)

                # Decide the best T1 to take and convert the file format into the BIDS standard
                t1_selected = bids.choose_correction(mods_folder_path, t1_priority, 'T1')
                if t1_selected != -1 and t1_selected != 0:
                    t1_sel_path = glob(path.join(mods_folder_path, '*' + t1_selected + '*', '*.nii*'))[0]
                    bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                else:
                    if t1_selected == -1:
                        logging.info('No T1 found for ' + mods_folder_path)
                        mmt.add_missing_mod('T1', ses)

                    else:
                        logging.warning('None of the desiderd T1 corrections is available for ' + mods_folder_path)

                # Convert T2FLAIR
                out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                if out == -1:
                    logging.warning('No FLAIR found for '+ mods_folder_path)
                    mmt.add_missing_mod('FLAIR', ses)

                # Convert fMRI
                out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name)
                if out == -1:
                        mmt.add_missing_mod('fMRI', ses)

        logging.info("Conversion for the subject terminated.\n")

    print 'Converting the clinical data...'
    convert_clinical(source_dir, dest_dir, bids_ids, original_insight_ids, original_insight_ids_alpha)
    # Printing the statistics about missing modalities into a file
    print_statistics(summary, len(subjs_paths), ses_available, mmt)
    summary.close()