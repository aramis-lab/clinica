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


def convert(source_dir, dest_dir):
    """
    Convert the CAPP dataset into the BIDS standard.

    Args:
        source_dir: directory of the input dataset
        dest_dir: output directory
    """
    capp_spath = []
    capp_ids = []
    bids_ids = []
    t1_to_consider = ['3DT1_noPN_DIS', '3DT1_noSCIC_GW', '3DT1_noCLEAR_GEO', '3DT1_CLEAR_GEO', '3DT1_S']
    ses_available = ["M00", "M18"]
    map_problems = ["13001PBA20150623M18B0MAPph_S016.nii.gz", "13002PRJ20150922M18B0MAPph_S014.nii.gz",
                    "11001PGM20130704M00B0MAPph_S010.bval","07002PPP20150116M18B0MAPph_S009.nii.gz",
                    "07003PGM20141217M18B0MAPph_S010.nii.gz"]
    mmt = MissingModsTracker(ses_available)
    os.mkdir(dest_dir)
    summary_file = open(path.join(dest_dir, 'conversion_summary.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'conversion.log'), format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %I:%M' )
    participants = open(path.join(dest_dir, 'participants.tsv'), 'w')
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

    participants.write("------------------------------\n")
    participants.write("participant_id   BIDS_id\n")
    participants.write("------------------------------\n")

    for c_ids, b_ids in zip(capp_ids, bids_ids):
        participants.writelines(c_ids+"     "+b_ids+'\n')
    print "*******************************"
    print "CAPP to BIDS converter"
    print "*******************************"
    print "Number of subjects: ", len(capp_spath)
    # For each subject extract the list of files and convert them into BIDS specification
    for sub_path in capp_spath:
        print "Converting:", sub_path
        logging.info("Converting:"+sub_path)
        sessions = glob(path.join(sub_path, "*"))

        # Check if the subject has some missing session
        for s in ses_available:
            if not any(s in ses_found for ses_found in sessions):
                mmt.incr_missing_session(s)

        # For each sub-session
        for ses_path in sessions:
            ses = ses_path.split(os.sep)[-1]

            if 'rescan' in ses:
                logging.warning('Rescan of a session found: ' + ses + '. Ignored.')
                print 'Rescan of a session found: ' + ses + '. Ignored.'
                continue

            if 'discard' in ses:
                continue

            if os.path.exists(path.join(sub_path,ses)):
                # Extracting the index of the subject
                subj_index = capp_spath.index(sub_path)
                os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses))
                ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses)
                bids_file_name = bids_ids[subj_index]+'_ses-'+ses
                session_path = path.join(sub_path,ses)
                mods_folder_path = path.join(session_path, 'NIFTI')

                # Extract and convert fieldmap data
                out = bids.convert_fieldmap(mods_folder_path, path.join(ses_dir_bids, 'fmap'), bids_file_name, map_problems)
                if out == -1:
                    mmt.add_missing_mod('Fieldmap', ses)
                    logging.warning("No Fieldmap found for "+mods_folder_path)
                elif out == 0:
                    logging.warning("Map or MapPh missing for " + mods_folder_path + ': skipped')
                    mmt.add_incomplete_mod('Fieldmap', ses)


                # Decide which T1 will be considered for the conversion
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
                out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('FLAIR', ses)
                    logging.warning('No FLAIR found for '+ mods_folder_path)


                # Merge and convert all valid DTI folders for the same subject
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
                out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('fMRI', ses)

        logging.info("Conversion for the subject terminated.\n")

    # Printing the statistics about missing modalities into a file
    print_statistics(summary_file, len(capp_spath), ses_available, mmt)
    summary_file.close()