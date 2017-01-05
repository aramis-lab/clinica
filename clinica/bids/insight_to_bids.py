"""
Covert the INSIGHT dataset into the BIDS specification.


@author: Sabrina Fontanella
"""

from os import path
from glob import glob
import os
import logging
import bids_utils as bids
import pandas as pd
from converter_utils import  remove_space_and_symbols

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"

def convert_clinical(input_path, out_path, bids_ids, original_ids, orig_ids_alpha, memento_ids):
    '''
    Convert clinical data of INSIGHT (CATI organized) dataset (subject included only).

    Args:
        input_path:
        out_path:
        bids_ids:
        original_ids:
    :return:

    '''

    memento_ids_alpha = remove_space_and_symbols(memento_ids)
    use_memento = False
    fields_bids = ['participant_id']
    fields_dataset = []
    clinic_specs_path = path.join(os.path.dirname(__file__), 'data', 'clinical_specifications.xlsx')

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
            # ToDO: it works but maybe there a way much more elegant!
            index_to_delete.append((participant_df[participant_df['alternative_id_1'] == participant_df['alternative_id_1'][i]].index.tolist())[0])
        else:
            bid = [s for s in bids_ids if value in s][0]
            participant_df['participant_id'][i] = bid

    # Remove all the informations regarding the subjects not included from the dataframe
    index_to_delete.append((participant_df[participant_df['alternative_id_2'] == '0040052_HEPH'].index.tolist())[0])
    participant_df = participant_df.drop(index_to_delete)
    # Remove the row with for subject 0040052_HEPH (non included in the list of available)
    participant_df.to_csv(path.join(out_path,'participants.tsv'), sep='\t', index = False)

    # --Creation sessions.tsv--
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
            tmp = field_location[i].split('/')
            location = tmp[0]
            sheet = tmp[1]
            file_to_read_path = path.join(input_path, 'clinicalData', location)
            file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]

                # Extract the subject id
                if 'Subject' in row:
                    subj_id = row['Subject']
                    use_memento = False
                if 'PAT_ID' in row:
                    subj_id = row['PAT_ID']
                    subj_id = remove_space_and_symbols(subj_id)
                    use_memento = True
                if 'subject' in row:
                    subj_id = row['subject']
                    use_memento = True


                # Extract the correspondant BIDS id and create the output file if doesn't exist
                if use_memento:
                    subj = [s for s in memento_ids_alpha if s in subj_id]
                else:
                    subj = [s for s in orig_ids_alpha if s in subj_id]


                # If the ID is an included one
                if len(subj) != 0:
                    subj = subj[0]
                    if use_memento:
                        subj_index = memento_ids_alpha.index(subj)
                    else:
                        subj_index = orig_ids_alpha.index(subj)
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields_db[i]]

                    # Create the files
                    if sessions_dict.has_key(bids_ids[subj_index]):
                        (sessions_dict[bids_ids[subj_index]]['M0']).update({sessions_fields_bids[i]: row[sessions_fields_db[i]]})
                    else:
                        sessions_dict.update({bids_ids[subj_index]: {'M0' : {'session_id': 'M0', sessions_fields_bids[i]: row[sessions_fields_db[i]]}}})



    # Create the several sessions.tsv file extracting data from the dictionary
    fields_bids.insert(0, 'session_id')
    cols = []
    for bids_id in bids_ids:
        sessions_df = pd.DataFrame(columns=fields_bids)

        if sessions_dict.has_key(bids_id):
            session_df = pd.DataFrame(sessions_dict[bids_id]['M0'], index=['i', ])
            cols = session_df.columns.tolist()
            cols.remove('session_id')
            cols.insert(0, 'session_id')
            session_df = session_df[cols]
            session_df.to_csv(path.join(out_path,bids_id, bids_id + '_sessions.tsv'), sep='\t', index = False)
        else:
            print bids_id, ' not found into clinical data'

    # -- Creation of scans.tsv --
    scans = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
    scans_fields_db = scans['INSIGHT']
    field_location = scans['INSIGHT location']
    scans_fields_bids = scans['BIDS CLINICA']
    fields_dataset = []
    fields_bids = ['filename']

    for i in range(0, len(scans_fields_db)):
        if not pd.isnull(scans_fields_db[i]):
            fields_bids.append(scans_fields_bids[i])
            fields_dataset.append(scans_fields_db[i])

    scans_df = pd.DataFrame(columns=(fields_bids))

    for bids_id in bids_ids:
        sessions_paths = glob(path.join(out_path, bids_id, 'ses-*'))

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




def convert(source_dir, dest_dir, mod=''):
    # List of the all the subjects to take in consideration for the conversion
    subjs_list = path.join(source_dir,'clinicalData', '_INSIGHT-AD_IDs_INCLUDED__Status-07-06-2016.xlsx')

    data_dir = path.join(source_dir, 'convertData', 'study', 'Paris') # Folder where data are stored
    t1_priority = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    #mmt = MissingModsTracker(['M0', 'M24'])
    ses_available = ['M0', 'M24']

    insight_ids = []
    original_insight_ids = []
    original_insight_ids_alpha = []
    memento_ids = []
    subjs_paths = []
    bids_ids = []
    subjs_included = pd.read_excel(subjs_list, sheetname='IDs INCLUDED')

    if mod!='-c':
        os.mkdir(dest_dir)
        #summary = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
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
        if mod!='-c':
            os.mkdir(path.join(dest_dir, bids_ids[-1]))
        # In the INSIGHT cati organized the subject is is the insight id without the character -
        subj_path = glob(path.join(data_dir, '*' + insight_ids[-1] + '*'))[0]
        subjs_paths.append(subj_path)


    # For each subject extract the list of files and convert them into BIDS specification
    if mod!='c':
        for subj_path in subjs_paths:
            logging.info('Converting:' + subj_path)
            # Extract the session(s) available
            sessions = glob(path.join(subj_path, "*"))

            # Check if the subject has some missing session
            #for s in ses_available:
                # if not any(s in ses_found for ses_found in sessions ):
                #     mmt.incr_missing_session(s)

            # For each sub-session
            for ses_path in sessions:
                ses = ses_path.split(os.sep)[-1]

                if 'rescan' in ses:
                    logging.warning('Rescan of a session found: ' + ses + 'for ' + subj_path + '. Ignored.')
                    print 'Rescan of a session found: ' + ses + '. Ignored.'
                    continue

                # Extracting the index of the subject
                subj_index = subjs_paths.index(subj_path)
                if mod!='-c':
                    os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses))
                ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses)
                bids_file_name = bids_ids[subj_index] + '_ses-' + ses
                session_path = path.join(subj_path, ses)
                mods_folder_path = path.join(session_path, 'NIFTI')

                if mod != "-c":
                    # Convert the fieldmap data
                    out = bids.convert_fieldmap(mods_folder_path, path.join(ses_dir_bids, "fmap"), bids_file_name)
                    if out == -1:

                        logging.warning("No Fieldmap found for "+mods_folder_path)
                    elif out == 0:
                        logging.warning("Map or MapPh missing for "+mods_folder_path +': skipped.')


                    #Convert DTI data
                    out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                    if out != None:  # missing or incomplete DTI
                        if out == -1:
                            logging.info('No DTI found for ' + mods_folder_path)
                        else:
                            for e in out:
                                logging.warning('.bvec or .bval not found for DTI folder ' + e + ': skipped.')

                    # Decide the best T1 to take and convert the file format into the BIDS standard
                    t1_selected = bids.choose_correction(mods_folder_path, t1_priority, 'T1')
                    if t1_selected != -1 and t1_selected != 0:
                        t1_sel_path = glob(path.join(mods_folder_path, '*' + t1_selected + '*', '*.nii*'))[0]
                        bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                    else:
                        if t1_selected == -1:
                            logging.info('No T1 found for ' + mods_folder_path)


                        else:
                            logging.warning('None of the desiderd T1 corrections is available for ' + mods_folder_path)

                    # Convert T2FLAIR
                    out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                    if out == -1:
                        logging.warning('No FLAIR found for '+ mods_folder_path)


                    # Convert fMRI
                    out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name)
                    if out == -1:
                        logging.warning('No fMRI found for ' + mods_folder_path)


            logging.info("Conversion for the subject terminated.\n")

    print 'Converting the clinical data...'
    convert_clinical(source_dir, dest_dir, bids_ids, original_insight_ids, original_insight_ids_alpha, memento_ids)


