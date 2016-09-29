"""
Covert the CAPP dataset into the BIDS specification.

ToDo:
- Integrate a list that indicates, in case of rescan, what is the best version to consider
- Provide support for the metadata
- Improve the mergeDTI function merging the .bvec and .bval files with python instead that with bash
- Fix the code style
- Fix the problem with python logging and nipype

@author: Sabrina Fontanella
"""

from os import path
from glob import glob
import os
from shutil import copy
# import nipype
import logging
import nibabel as nib
from missing_mods_tracker import MissingModsTracker
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
    sessions = ["M00", "M18"]
    map_problems = ["13001PBA20150623M18B0MAPph_S016.nii.gz", "13002PRJ20150922M18B0MAPph_S014.nii.gz",
                    "11001PGM20130704M00B0MAPph_S010.bval","07002PPP20150116M18B0MAPph_S009.nii.gz",
                    "07003PGM20141217M18B0MAPph_S010.nii.gz"]
    bids_suff = {
        'T1': '_T1w',
        'T2': '_T2w',
        'Flair': '_FLAIR',
        'MapPh': '_phasediff',
        'Map': '_magnitude',
        'fMRI': '_bold',
        'dwi': '_dwi'
    }

    mmt = MissingModsTracker()
    os.mkdir(dest_dir)
    readme = open(path.join(dest_dir, 'README.txt'), 'w')
    logging.basicConfig(filename=path.join(dest_dir, 'logfile.log'), format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %I:%M' )
    participants = open(path.join(dest_dir,'participants.tsv'), 'w')
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
        # For each sub-session
        for ses in sessions:
            if os.path.exists(path.join(sub_path,ses)):
                # Extracting the index of the subject
                subj_index = capp_spath.index(sub_path)
                os.mkdir(path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses))
                ses_dir_bids = path.join(dest_dir, bids_ids[subj_index], 'ses-'+ses)
                bids_file_name = bids_ids[subj_index]+'_ses-'+ses
                session_path = path.join(sub_path,ses)
                mods_folder_path = path.join(session_path, 'NIFTI')

                # Convert the fieldmap data
                map = bids.remove_rescan(glob(path.join(sub_path, ses, "NIFTI", "*MAP_*", '*.nii.gz')))
                map_ph = bids.remove_rescan(glob(path.join(sub_path, ses, "NIFTI", "*MAPph_*", '*.nii.gz')))
                if len(map) == 0:
                    logging.warning('Missing magnitude image(s) for '+mods_folder_path)
                    mmt.add_missing_mod('Map')
                if len(map_ph) == 0:
                    logging.warning('Missing phase image(s) for '+mods_folder_path)
                    mmt.add_missing_mod('MapPh')
                # If the information regarding the Fieldmap data are complete
                if len(map) > 0 and len(map_ph) > 0:
                    map_ph_name = map_ph[0].split(os.sep)[-1]
                    map_name = map[0].split(os.sep)[-1]
                    # toSolve: there are some files that produce an error when loaded with Nibabel!
                    if (map_ph_name not in map_problems) and (map_name not in map_problems):
                        # Open the files with Nibabel
                        map_nib = nib.load(map[0])
                        map_ph_nib = nib.load(map_ph[0])
                        dim_map = (map_nib.header['dim'])[4]
                        dim_map_ph = (map_ph_nib.header['dim'])[4]
                        os.mkdir(path.join(ses_dir_bids, "fmap"))
                        # Case 1: one phase difference image and at least one magnitude image
                        if dim_map_ph ==1 and dim_map >0:
                            copy(map_ph[0], path.join(ses_dir_bids, "fmap", bids_file_name + bids_suff['MapPh'] + '.nii.gz'))
                            os.system('fslsplit ' + map[0]+ ' '+path.join(ses_dir_bids, "fmap/"+bids_file_name+bids_suff['Map']))
                        # Case 2: two phase images and two magnitude images'
                        elif dim_map_ph == 2 and dim_map == 2:
                            os.system('fslsplit ' + map_ph[0] + ' ' + path.join(ses_dir_bids, "fmap/"+bids_file_name+bids_suff['MapPh']))
                            os.system('fslsplit ' + map[0] + ' ' + path.join(ses_dir_bids,"fmap/" + bids_file_name + bids_suff['Map']))

                # Decide which T1 will be considered for the conversion
                t1_selected = bids.choose_correction(mods_folder_path, t1_to_consider, 'T1')
                if t1_selected != -1 and t1_selected != 0:
                    t1_sel_path = glob(path.join(mods_folder_path, '*' + t1_selected + '*', '*.nii*'))[0]
                    bids.convert_T1(t1_sel_path, path.join(ses_dir_bids, 'anat'), bids_file_name)
                else:
                    if t1_selected == -1:
                        logging.info('No T1 found for ' + mods_folder_path)
                        mmt.add_missing_mod('T1')
                    else:
                        logging.warning('None of the desiderd T1 corrections is available for ' + mods_folder_path)

                # Extract and convert the T2 FLAIR modality if is available
                out = bids.convert_flair(mods_folder_path, path.join(ses_dir_bids, "anat"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('FLAIR')

                # Merge and convert all valid DTI folders for the same subject
                out = bids.merge_DTI(mods_folder_path, path.join(ses_dir_bids, 'dwi'), bids_file_name)
                if out != None: #missing or incomplete DTI
                    if out == -1:
                        logging.info('No DTI found for ' + mods_folder_path)
                        mmt.add_missing_mod('DTI')
                    else:
                        for e in out:
                            logging.warning('.bvec or .bval not found for DTI folder ' + e + ' Skipped')

                # Extract and convert the fMRI
                out = bids.convert_fmri(mods_folder_path, path.join(ses_dir_bids, "func"), bids_file_name)
                if out == -1:
                    mmt.add_missing_mod('fMRI')

        logging.info("Conversion for the subject terminated.\n")

    # Printing the statistics about missing modalities into a file
    readme.write('Number of subjects converted: '+str(len(capp_spath))+'\n')
    readme.write('********************************\n')
    readme.write('Number of missing files per modalities (more details inside the logfile):\n')
    mmt.remove_mods_unused()
    missing_list = mmt.get_missing_list()
    for key in missing_list:
        readme.write(key+':'+str(missing_list[key])+"\n")
    readme.write('********************************\n')
    readme.close()