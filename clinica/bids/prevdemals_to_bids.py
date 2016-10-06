"""
Covert the PREVDEMALS dataset into the BIDS specification.

ToDo:
- Check if the subject list is the same for both projects
- Implement the log file

@author: Sabrina Fontanella
"""
from os import path
from glob import glob
import os
import logging
import bids_utils as bids
from shutil import copy
from converter_utils import MissingModsTracker, print_statistics

def convert(source_dir, dest_dir):
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
    mmt = MissingModsTracker(['M0', 'M24'])
    os.mkdir(dest_dir)
    os.mkdir(path.join(dest_dir, 'GENFI'))
    os.mkdir(path.join(dest_dir, 'ICM'))
    logging.basicConfig(filename=path.join(dest_dir, 'conversion.log'), format='%(asctime)s %(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)
    print "*******************************"
    print "PREVDEMALS to BIDS converter"
    print "*******************************"
    # Convert all the files contained in the two project folder PREV_DEMALS_GENFI and PREV_DEMALS_ICM
    for proj in projects:
        cities_folder = glob(path.join(projects[proj], '*'))
        pda_spath = []
        pda_ids = []
        bids_ids = []
        dest_dir_proj = path.join(dest_dir, proj)
        participants = open(path.join(dest_dir_proj, 'participants.tsv'), 'w')
        participants.write("------------------------------\n")
        participants.write("participant_id   BIDS_id\n")
        participants.write("------------------------------\n")
        for cf in cities_folder:
            for subj_path in glob(path.join(cf, "*")):
                subj_id = subj_path.split(os.sep)[-1]
                pda_spath.append(subj_path)
                pda_ids.append(subj_id)
                bids_ids.append('sub-PREVDELMALS' + subj_id)
                # Create the subject folder in the BIDS converted dataset
                os.mkdir(path.join(dest_dir_proj, bids_ids[-1]))

        for p_ids, b_ids in zip(pda_ids, bids_ids):
            participants.writelines(p_ids + "     " + b_ids + '\n')

        participants.close()

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
                    if out != None: #missing or incomplete DTI
                        if out == -1:
                            logging.info('No DTI found for ' + mods_folder_path)
                            mmt.add_missing_mod('DTI', ses)
                        else:
                            for e in out:
                                logging.warning('.bvec or .bval not found for DTI folder ' + e + ' Skipped')
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

                # Extract and convert the T2 FLAIR modality if is available
                out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('FLAIR', ses)

                logging.info("Conversion for the subject terminated.\n")
