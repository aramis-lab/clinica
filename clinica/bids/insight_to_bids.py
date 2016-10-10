"""
Covert the INSIGHT dataset into the BIDS specification.

ToDo:
- Subjects missing problems and multiple t2flairs

@author: Sabrina Fontanella
"""

from os import path
from os import path
from glob import glob
import os
from shutil import copy
import logging
import nibabel as nib
from converter_utils import MissingModsTracker,print_statistics
import bids_utils as bids
import pandas as pd
import numpy as np


def converter(source_dir, dest_dir):
    # List of the all the subjects to take in consideration for the conversion
    subjs_list = path.join(source_dir, '_INSIGHT-AD_IDs_INCLUDED__Status-07-06-2016.xlsx')

    source_dir = path.join(source_dir, 'convertData', 'study', 'Paris') # Folder where data are stored
    t1_priority = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    mmt = MissingModsTracker(['M0', 'M24'])
    ses_available = ['M0', 'M24']


    insight_ids = []
    memento_ids = []
    memento_subjs_paths = []
    subjs_paths = []
    bids_ids = []
    subjs_included = pd.read_excel(subjs_list, sheetname='IDs INCLUDED')

    os.mkdir(dest_dir)
    participants = open(path.join(dest_dir, 'participants.tsv'), 'w')
    summary = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'conversion'
                                                     '.log'), format='%(asctime)s %(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)

    print "*******************************"
    print "INSIGHT to BIDS converter"
    print "*******************************"

     # The last 7 rows from the excel files are not ids but other datas
    for i in range(0, len(subjs_included.index) - 7):
        insight_ids.append((str(subjs_included['INSIGHT IDs'][i]).replace("-", "")).replace(" ", ""))
        memento_ids.append(str(subjs_included['MEMENTO IDs'][i]).replace("-", ""))
        # The BIDS specification allows only numbers and letter inside the subject id
        bids_ids.append('sub-INSIGHT'+insight_ids[-1].replace("-", ""))
        os.mkdir(path.join(dest_dir, bids_ids[-1]))
        subj_path = glob(path.join(source_dir, '*' + insight_ids[-1] + '*'))[0]
        subjs_paths.append(subj_path)


    # For each subject extract the list of files and convert them into BIDS specification
    for subj_path in subjs_paths:
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
                logging.warning('Rescan of a session found: ' + ses + '. Ignored.')
                print 'Rescan of a session found: ' + ses + '. Ignored.'
                continue

            # Extracting the index of the subject
            subj_index = subjs_paths.index(subj_path)
            os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses))
            ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-' + ses)
            bids_file_name = bids_ids[subj_index] + '_ses-' + ses
            session_path = path.join(subj_path, ses)
            mods_folder_path = path.join(session_path, 'NIFTI')

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

    # Printing the statistics about missing modalities into a file
    print_statistics(summary, len(subjs_paths), ses_available, mmt)
    summary.close()
    participants.close()








