from os import path
from glob import glob
import pandas as pd
import os

def find_mod_available(dataset_dir):
    mods_dict = {}
    subjects_paths_lists = glob(path.join(dataset_dir, '*sub-*'))
    for sub_path in subjects_paths_lists:
        ses_paths = glob(path.join(sub_path, '*ses-*'))
        for session in ses_paths:
            mods_paths_folders = glob(path.join(session, '*/'))
            if ['func' in m for m in mods_paths_folders]:
                if not mods_dict.has_key('func'):
                    mods_dict.update({'func': 'func'})
            if ['dwi' in m for m in mods_paths_folders]:
                if not mods_dict.has_key('dwi'):
                    mods_dict.update({'dwi': 'dwi'})
            if ['anat' in m for m in mods_paths_folders]:
                anat_files_paths = glob(path.join(session, 'anat', '*'))

                for anat_file in anat_files_paths:
                    anat_name = anat_file.split(os.sep)[-1]

                    # Extract the name of the file without the extension
                    if '.nii.gz' in anat_name:
                        anat_name = anat_name.replace('.nii.gz', '')
                        anat_ext = 'nii.gz'
                    else:
                        anat_name = os.path.splitext(anat_name.split(os.sep)[-1])[0]
                        anat_ext = os.path.splitext(anat_name.split(os.sep)[-1])[1]

                    if anat_ext != 'json':
                        file_parts = anat_name.split("_")
                        anat_type = file_parts[len(file_parts) - 1]

                        if mods_dict.has_key('anat'):
                            if anat_type not in mods_dict['anat']:
                                anat_aval = mods_dict['anat']
                                anat_aval.append(anat_type)
                                mods_dict.update({'anat': anat_aval})
                        else:
                            mods_dict.update({'anat': [anat_type]})


    print mods_dict

def create_subs_sess_list(dataset_path, out_dir):
    file_name = out_dir.split(os.sep)[-1]
    if len(file_name) == 0:
        file_name = 'subjects_sessions_list.tsv'
    else:
        # Extract the path of the file
        out_dir = os.path.dirname(out_dir)

    if '.' not in file_name:
        file_name = file_name+'.tsv'

    subjs_sess_tsv = open(path.join(out_dir, file_name), 'w')
    subjs_sess_tsv.write('participant_id' + '\t' + 'session_id' + '\n')

    subjects_paths = glob(path.join(dataset_path, '*sub-*'))
    if len(subjects_paths) == 0:
        raise 'Dataset empty or not BIDS-compliant.'
    for sub_path in subjects_paths:
        subj_id = sub_path.split(os.sep)[-1]
        sess_list = glob(path.join(sub_path, '*ses-*'))
        for ses_path in sess_list:
            session_name = ses_path.split(os.sep)[-1]
            subjs_sess_tsv.write(subj_id+'\t'+session_name+'\n')

    subjs_sess_tsv.close()





