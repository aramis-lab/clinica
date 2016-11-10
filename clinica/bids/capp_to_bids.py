"""
Covert the CAPP dataset into the BIDS specification.

ToDo:
- Integrate a list that indicates, in case of rescan, what is the best version to consider
- Provide support for the metadata

@author: Sabrina Fontanella
"""

from os import path
from glob import glob
import os
import logging
from converter_utils import MissingModsTracker, print_statistics
import bids_utils as bids
import pandas as pd
import json


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
         file_to_convert: name of the folder to choose for the conversion
    """
    if subj_id in special_list:
        if mod in special_list[subj_id]:
            if ses in special_list[subj_id][mod]:
                file_to_convert = special_list[subj_id][mod]
                return file_to_convert
            else:
                return False
    return False


def convert_clinical(input_path, out_path, bids_ids):
    pd.options.mode.chained_assignment = None
    fields_bids = ['participant_id']
    fields_dataset = []
    prev_location = ""
    prev_sheet = ""
    capp_origin_ids = []
    clinic_specs_path = path.join(os.path.dirname(__file__), 'bids-utils-docs', 'clinical_specifications.xlsx')

    # -- Creation of participant.tsv --
    participants_specs = pd.read_excel(clinic_specs_path, sheetname='participant.tsv')
    # Extracts information regarding INSIGHT dataset
    participant_fields_db = participants_specs['CAPP']
    field_location = participants_specs['CAPP location']
    participant_fields_bids = participants_specs['BIDS CLINICA']

    # Extract the list of the available fields for the CAPP dataset (and the corresponding BIDS version)
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
            sheet = tmp[1]
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_to_read_path = path.join(input_path, 'clinicalData', location + '.xlsx')
                file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, 114):
                field_col_values.append(file_to_read.get_value(j, participant_fields_db[i]))
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        value = participant_df['alternative_id_1'][i]
        # Removing _M0
        value = value[:-4]
        bids_id = [s for s in bids_ids if value in s]
        if len(bids_id)==0:
            print "Subject "+value+" not found in the BIDS converted version of the dataset."
            logging.error("Subject "+value+" not found in the BIDS converted version of the dataset.")
        else:
            participant_df['participant_id'][i] = bids_id[0]

    participant_df.to_csv(path.join(out_path, 'participants.tsv'), sep='\t', index=False)

    # -- Creation of *_scans.tsv --
    scans_dict = {}
    for bids_id in bids_ids:
        scans_dict.update({bids_id: {'T1/DWI/fMRI': {}, 'FDG': {}}})
    scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
    scans_fields_db = scans_specs['CAPP']
    field_location = scans_specs['CAPP location']
    scans_fields_bids = scans_specs['BIDS CLINICA']
    scans_fields_mod = scans_specs['Modalities related']


    # For each field available extract the original name, extract from the file all the values and fill a data structure
    for i in range(0, len(scans_fields_db)):
        if not pd.isnull(scans_fields_db[i]):
            tmp = field_location[i].split('/')
            location = tmp[0]
            sheet = tmp[1]
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_to_read_path = path.join(input_path, 'clinicalData', location + '.xlsx')
                file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
                prev_location = location
                prev_sheet = sheet
            for bids_id in bids_ids:
                cati_id = bids_id[8:]+'_M00'
                row_to_extract = file_to_read[file_to_read['Identifiant CATI'] == cati_id].index.tolist()
                if len(row_to_extract)>0:
                    row_to_extract = row_to_extract[0]
                    # fill the dictionary with all the information
                    (scans_dict[bids_id][scans_fields_mod[i]]).update({scans_fields_bids[i]: file_to_read.iloc[row_to_extract][scans_fields_db[i]]})

    scans_df = pd.DataFrame(columns=('filename', 'acq_time'))
    # Extract the information for each modality and create the several tsv files
    for bids_id in bids_ids:
        # Create the file
        tsv_name = bids_id+"_ses-M00_scans.tsv"
        if os.path.exists(path.join(out_path,bids_id,'ses-M00',tsv_name)):
            os.remove(path.join(out_path,bids_id,'ses-M00',tsv_name))

        scans_tsv = open(path.join(out_path,bids_id,'ses-M00', tsv_name), 'a')
        scans_df.to_csv(scans_tsv, sep='\t')
        # Extract modalities available for each subject
        mod_available = glob(path.join(out_path, bids_id, 'ses-M00', '*'))
        for mod in mod_available:
            mod_name =  os.path.basename(mod)
            files = glob(path.join(mod,'*'))
            for file in files:
                file_name = os.path.basename(file)
                if mod == "anat" or mod == "dwi" or mod == "func":
                    type = 'T1/DWI/fMRI'
                else:
                    type = 'FDG'

                if (scans_dict[bids_id][type]).has_key('acq_time'):
                    scans_df = pd.DataFrame(scans_dict[bids_id][type], index=['i', ])
                    scans_df['filename'] = path.join(mod_name,file_name)
                    cols = scans_df.columns.tolist()
                    cols = cols[-1:] + cols[:-1]
                    scans_df = scans_df[cols]
                    scans_df.to_csv(scans_tsv, header=False, sep='\t', index=False)


def convert(source_dir, dest_dir, param=''):
    """
    Convert the CAPP dataset into the BIDS standard.

    Args:
        source_dir: directory of the input dataset
        dest_dir: output directory
        param: parameters for the conversion
    """
    capp_spath = []
    capp_ids = []
    bids_ids = []
    special_subjs_path = path.join(os.path.dirname(__file__), 'bids-utils-docs', 'special_subjs_CAPP.json')
    special_subjs_file = open(special_subjs_path, 'r')
    special_subjs_json = json.load(special_subjs_file)
    t1_to_consider = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    ses_available = ["M00", "M18"]
    mmt = MissingModsTracker(ses_available)

    print "*******************************"
    print "CAPP to BIDS converter"
    print "*******************************"

    if param != '-c':
        os.mkdir(dest_dir)
        summary_file = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
        logging.basicConfig(filename=path.join(dest_dir, 'conversion_modalities.log'),
                            format='%(asctime)s %(levelname)s:%(message)s',
                            level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M')
        cities_folder = glob(path.join(source_dir, '*', '*', '*'))


        # Create the lists of bids_ids extracting the list of the subjects from the dataset
        for cf in cities_folder:
            for sub_path in glob(path.join(cf, "*")):
                sub_id = sub_path.split(os.sep)[-1]
                # The substring 'xxx' inside a folder subject name indicates an acquisition to ignore
                if 'xxx' in sub_id:
                    logging.warning('Anomalous subject folder: '+sub_id+'. Ignored.\n')
                else:
                    capp_spath.append(sub_path)
                    capp_ids.append(sub_id)
                    bids_ids.append('sub-CAPP'+sub_id)
                    os.mkdir(path.join(dest_dir, bids_ids[-1]))

        print "Number of subjects: ", len(capp_spath)
        # For each subject extract the list of files and convert them into BIDS specification
        for sub_path in capp_spath:
            if param!='-c':
                print "Converting:", sub_path
                logging.info("Converting:"+sub_path)
            sessions = glob(path.join(sub_path, "*"))

            # Check if the subject has some missing session
            for s in ses_available:
                if not any(s in ses_found for ses_found in sessions):
                    mmt.incr_missing_session(s)

            # The subject 09001PMC it's a particular in which the rescaned version of M18 session is chosed instead of the M18
            if sub_path.split(os.sep)[-1] in "09001PMC":
                logging.warning("Predifiened sessions to convert for subject 09001PMC: M00 and M18_rescaned.")
                print("WARNING: Predifiened sessions to convert for subject 09001PMC: M00 and M18_rescaned.")
                sessions = ["M00", "M18_rescaned"]

            # For each sub-session
            for ses_path in sessions:
                ses = ses_path.split(os.sep)[-1]

                if 'rescan' in ses:
                    logging.warning('Rescan of a session found: ' + ses + ' for ' + sub_path + '.')
                    print 'Rescan of a session found: ' + ses + 'for ' + sub_path + '.'
                    # Remove the word 'rescan' from the session
                    ses = ses[:-9]

                if 'discard' in ses:
                    continue

                if os.path.exists(path.join(sub_path,ses)):
                    # Extracting the index of the subject
                    subj_index = capp_spath.index(sub_path)
                    os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses))
                    ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses)
                    bids_file_name = bids_ids[subj_index]+'_ses-'+ses
                    session_path = path.join(sub_path, ses)
                    mods_folder_path = path.join(session_path, 'NIFTI')

                if param!='-c':
                    # Extract and convert fieldmap data
                    fixed_map_to_convert = has_fixed(capp_ids[subj_index], ses, 'map', special_subjs_json)
                    fixed_map_ph_to_convert = has_fixed(capp_ids[subj_index], ses, 'map_ph', special_subjs_json)
                    if fixed_map_to_convert != False or fixed_map_ph_to_convert!=False:
                        logging.info("Fixed fielmap to convert for " + mods_folder_path)
                        out = bids.convert_fieldmap(mods_folder_path, path.join(ses_dir_bids, 'fmap'), bids_file_name, [fixed_map_to_convert, fixed_map_ph_to_convert])
                    else:
                        out = bids.convert_fieldmap(mods_folder_path, path.join(ses_dir_bids, 'fmap'), bids_file_name)
                        if out == -1:
                            mmt.add_missing_mod('Fieldmap', ses)
                            logging.warning("No Fieldmap found for "+mods_folder_path)
                        elif out == 0:
                            logging.warning("Map or MapPh missing for " + mods_folder_path + ': skipped')
                            mmt.add_incomplete_mod('Fieldmap', ses)


                    # Decide which T1 will be considered for the conversion
                    fixed_file_to_convert = has_fixed(capp_ids[subj_index], ses, 't1', special_subjs_json)
                    if fixed_file_to_convert != False:
                        logging.info("Fixed T1 to convert for " + mods_folder_path + ": " + fixed_file_to_convert)
                        t1_sel_path = glob(path.join(mods_folder_path, fixed_file_to_convert, '*.nii*'))[0]
                        bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                    else:
                        t1_selected = bids.choose_correction(mods_folder_path, t1_to_consider, 'T1')
                        if t1_selected != -1 and t1_selected != 0:
                            t1_sel_path = glob(path.join(mods_folder_path, '*' + t1_selected + '*', '*.nii*'))[0]
                            bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                        else:
                            if t1_selected == -1:
                                logging.info('No T1 found for ' + mods_folder_path)
                                mmt.add_missing_mod('T1', ses)
                            else:
                                logging.warning('None of the desiderd T1 corrections is available for ' + mods_folder_path)

                    # Extract and convert the T2 FLAIR modality if is available
                    # Check if the subject contains particular specification for the T2flair to choose
                    fixed_file_to_convert = has_fixed(capp_ids[subj_index], ses, 'flair', special_subjs_json)
                    if fixed_file_to_convert!=False:
                            logging.info("Fixed FLAIR to convert for " + mods_folder_path + ": " + fixed_file_to_convert)
                            out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name, fixed_file_to_convert)
                    else :
                        out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name, False)
                    if out == -1:
                        mmt.add_missing_mod('FLAIR', ses)
                        logging.warning('No FLAIR found for '+ mods_folder_path)


                    # Merge and convert all valid DTI folders for the same subject
                    fixed_file_to_convert = has_fixed(capp_ids[subj_index], ses, 'dwi', special_subjs_json)
                    if fixed_file_to_convert != False:
                        out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name, fixed_file_to_convert)
                    else:
                        out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                    if out != None: #missing or incomplete DTI
                        if out == -1:
                            logging.info('No DTI found for ' + mods_folder_path)
                            mmt.add_missing_mod('DTI', ses)
                        else:
                            for e in out:
                                logging.warning('.bvec or .bval not found for DTI folder ' + e + ' Skipped')
                                mmt.add_incomplete_mod('DTI', ses)

                    # Extract and convert the fMRI
                    fixed_file_to_convert = has_fixed(capp_ids[subj_index], ses, 'fmri', special_subjs_json)
                    if fixed_file_to_convert != False:
                        out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name, fixed_file_to_convert)
                    else:
                        out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name)
                    if out == -1:
                        mmt.add_missing_mod('fMRI', ses)

            logging.info("Conversion for the subject terminated.\n")
    else:
        print '** Conversion of clinical data only **'
        logging.basicConfig(filename=path.join(dest_dir, 'conversion_clinical.log'),
                            format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M')
        if not os.path.exists(dest_dir):
            print dest_dir, ' not found.'
            raise

        bids_ids = [d for d in os.listdir(dest_dir) if os.path.isdir(path.join(dest_dir,d))]

    print 'Converting clinical data...'
    logging.info('Converting clinical data...')
    convert_clinical(source_dir, dest_dir, bids_ids)
    # Printing the statistics about missing modalities into a file
    if param != '-c':
        print_statistics(summary_file, len(capp_spath), ses_available, mmt)
        summary_file.close()