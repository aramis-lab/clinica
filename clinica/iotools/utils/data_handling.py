from os import path
from glob import glob
import pandas as pd
import os
from ..converter_utils import MissingModsTracker, print_statistics

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Sabrina Fontanella"]
__license__ = ""
__version__ = "0.1.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"


def create_merge_file(bids_dir, out_dir, true_false_mode = False):
    """
    Merge all the .tsv files containing clinical data of a BIDS compliant dataset and store
    the result inside a .tsv file

    :param bids_dir: path to the BIDS folder
    :param out_dir: path to the output foler
    :param true_false_mode: if True convert all the binary values to True/False

    """

    col_list = []
    scans_dict = {}

    if not os.path.isfile(path.join(bids_dir, 'participants.tsv')):
        raise 'participants.tsv not found'
    participants_df = pd.read_csv(path.join(bids_dir, 'participants.tsv'), sep='\t')
    subjs_paths = glob(path.join(bids_dir, '*sub-*'))
    subjs_paths.sort()

    out_file_name = out_dir.split(os.sep)[-1]
    if len(out_file_name) == 0 or out_dir == '.':
        out_file_name = 'merge_tsv.tsv'
    else:
        # Extract the path of the file
        out_dir = os.path.dirname(out_dir)

    if '.' not in out_file_name:
        out_file_name = out_file_name + '.tsv'
    else:
        extension = os.path.splitext(out_file_name)[1]
        if extension != '.tsv':
            raise 'Output file must be .tsv.'

    if out_dir == '.':
        out_dir = os.getcwd()

    for col in participants_df.columns.values:
        col_list.append(col)

    merged_df = pd.DataFrame(columns=col_list)

    for sub_path in subjs_paths:
        sub_name = sub_path.split(os.sep)[-1]
        # For each subject, extract the relative row from the dataframe
        row_participant = participants_df[participants_df['participant_id'] == sub_name]
        # Open the sessions file related to the subject
        sessions_df = pd.read_csv(path.join(sub_path, sub_name+'_sessions.tsv'), sep='\t')

        # For each session found extract the information contained in the scans files
        for line in range(0, len(sessions_df)):
            # Extract and convert to a dictonary information regarding the session
            row_sessions = sessions_df.iloc[line]
            row_session_df = pd.DataFrame([row_sessions])
            new_cols = [s for s in row_session_df.columns.values if s not in col_list]
            if len(new_cols)!=0:
                for i in range(0, len(new_cols)):
                    col_list.append(new_cols[i])

            session_id = row_sessions['session_id']
            if os.path.isfile(path.join(bids_dir, sub_name, 'ses-'+session_id, sub_name+'_'+'ses-'+session_id+'_scans.tsv')):
                scans_df = pd.read_csv(path.join(bids_dir, sub_name, 'ses-'+session_id, sub_name+'_'+'ses-'+session_id+'_scans.tsv'), sep='\t')
                for i in range(0, len(scans_df)):
                    for col in scans_df.columns.values:
                        if col == 'filename':
                            pass
                        else:
                            file_scan = scans_df.iloc[i]['filename']
                            file_name = file_scan.split('/')[1]
                            # Remove the extension .nii.gz
                            file_name = os.path.splitext(os.path.splitext(file_name)[0])[0]
                            file_parts = file_name.split('_')
                            last_pattern_index = len(file_parts) - 1
                            mod_type = file_parts[last_pattern_index]
                            value = scans_df.iloc[i][col]
                            new_col_name = col+'_'+mod_type
                            scans_dict.update({new_col_name:value})
                row_scans = pd.DataFrame(scans_dict, index=[0])
            else:
                row_scans = pd.DataFrame()

            new_cols = [s for s in row_scans.columns.values if s not in col_list]
            if len(new_cols)!=0:
                for i in range(0, len(new_cols)):
                    col_list.append(new_cols[i])

            row_to_append_df= pd.DataFrame(columns=row_participant.columns)
            for col in row_participant:
                row_to_append_df[col] = row_participant[col]

            # Append all the data inside session_df
            for col in row_session_df:
                row_to_append_df[col] = row_session_df[col].values[0]

            for col in row_scans:
                row_to_append_df[col] = row_scans[col].values[0]

            merged_df = merged_df.append(row_to_append_df)
        scans_dict = {}

    old_index = col_list.index('session_id')
    col_list.insert(1, col_list.pop(old_index))
    merged_df = merged_df[col_list]
    merged_df.to_csv(path.join(out_dir, out_file_name), sep='\t', index=False)

    # Call the script for computing the missing modalities and append the result to the merged file
    compute_missing_mods(bids_dir, out_dir, 'tmpG7VIY0')
    tmp_ses = glob(path.join(out_dir, 'tmpG7VIY0*'))
    for f in tmp_ses:
        # Skip the summary file
        if not 'summary' in f:
            # Load the file
            mss_df = pd.read_csv(f, sep='\t')
            f_name = f.split(os.sep)[-1]
            patterns = f_name.split('-')
            ses_id = patterns[len(patterns)-1]
            ses_id = ses_id.replace('.tsv', '')
            cols = mss_df.columns.values

            # If the file opened contains new columns, add them to the exstisting merged_df
            for col_name in cols:
                if not col_name in col_list:
                    merged_df[col_name] = 0

            for i in range(0, len(mss_df)):
                row = mss_df.iloc[i]
                subj_idx = merged_df[ (merged_df['participant_id'] == row['participant_id']) & (merged_df['session_id'] == ses_id)].index.tolist()


                if len(subj_idx)>1:
                    raise ValueError('Multiple row for the same visit in the merge-tsv file.')
                elif len(subj_idx) ==0:
                    print 'Warning: Found modalities missing information for the subject:' + row['participant_id']+ ' visit:' + ses_id + ' but the subject is not included in the column participant_id.'
                    continue
                else:
                    subj_idx = subj_idx[0]
                for col_name in cols:
                    if not col_name=='participant_id':
                        merged_df.iloc[subj_idx, merged_df.columns.get_loc(col_name)] = row[col_name]

    # Remove all the temporary files created
    for f in tmp_ses:
        os.remove(f)

    if true_false_mode:
        merged_df = merged_df.replace(['Y','N'], ['0', '1'])

    merged_df.to_csv(path.join(out_dir, out_file_name), sep='\t', index=False)


def find_mods_and_sess(dataset_dir):
    """
    Finds all the modalities and sessions available for a given BIDS dataset

    :param dataset_dir:
    :return: a dictionary that stores the sessions and modalities found and has the following structure.
    Example:
    {
        'sessions': ['ses-M00', 'ses-M18'],
        'fmap': ['fmap'],
        'anat': ['flair', 't1w'],
        'func': ['func_task-rest'],
        'dwi': ['dwi']
    }
    """

    mods_dict = {}
    mods_list = []
    subjects_paths_lists = glob(path.join(dataset_dir, '*sub-*'))

    for sub_path in subjects_paths_lists:
        ses_paths = glob(path.join(sub_path, '*ses-*'))
        for session in ses_paths:
            ses_name = session.split(os.sep)[-1]
            mods_aval = []
            if mods_dict.has_key('sessions'):
                if not ses_name in mods_dict['sessions']:
                    mods_dict['sessions'].append(ses_name)
            else:
                mods_dict.update({'sessions':[ses_name]})
            mods_paths_folders = glob(path.join(session, '*/'))

            for p in mods_paths_folders:
                p = p[:-1]
                mods_aval.append(p.split('/').pop())

            if 'func' in mods_aval:
                list_funcs_paths = glob(path.join(session, 'func', '*bold.nii.gz'))
                for func_path in list_funcs_paths:
                    func_name = func_path.split(os.sep)[-1]
                    func_name_tokens = func_name.split('_')
                    func_task = func_name_tokens[2]
                if mods_dict.has_key('func'):
                    if not 'func_'+func_task in mods_dict['func']:
                        mods_dict['func'].append('func_'+func_task)
                else:
                    mods_dict.update({'func': ['func_'+func_task]})

                if not 'func_'+func_task in mods_list:
                    mods_list.append('func_'+func_task)

            if 'dwi' in mods_aval:
                if not mods_dict.has_key('dwi'):
                    mods_dict.update({'dwi': ['dwi']})
                if not 'dwi' in mods_list:
                    mods_list.append('dwi')

            if 'fmap' in mods_aval:
                if not mods_dict.has_key('fmap'):
                    mods_dict.update({'fmap': ['fmap']})
                if not 'fmap' in mods_list:
                    mods_list.append('fmap')

            if 'pet' in mods_aval:
                if not mods_dict.has_key('pet'):
                    mods_dict.update({'pet': ['pet']})
                if not 'pet' in mods_list:
                    mods_list.append('pet')

            if 'anat' in mods_aval:
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
                        anat_type = str.lower(file_parts[len(file_parts) - 1])
                        if mods_dict.has_key('anat'):
                            if anat_type not in mods_dict['anat']:
                                anat_aval = mods_dict['anat']
                                anat_aval.append(anat_type)
                                mods_dict.update({'anat': anat_aval})
                        else:
                            mods_dict.update({'anat': [anat_type]})

                        if anat_type not in mods_list:
                            mods_list.append(anat_type)

    return mods_dict


def compute_missing_mods(in_dir, out_dir, output_prefix = ''):
    """
    Compute the list of missing modalities for each subject in a BIDS compliant dataset.

    :param in_dir: path to the BIDS directory
    :param out_dir: path to the output directory
    :param output_prefix: string that replace the default prefix ('missing_mods_') in the name of all the output files
    created
    """

    # Find all the modalities and sessions available for the input dataset
    mods_and_sess= find_mods_and_sess(in_dir)
    sessions_found = mods_and_sess['sessions']
    mods_and_sess.pop('sessions')
    mods_avail_dict = mods_and_sess
    mods_avail = [j for i in mods_avail_dict.values() for j in i]
    cols_dataframe = mods_avail[:]
    cols_dataframe.insert(0, 'participant_id')
    mmt = MissingModsTracker(sessions_found, mods_avail)

    if output_prefix == '':
        out_file_name = 'missing_mods_'
    else:
        out_file_name = output_prefix + '_'

    summary_file = open(path.join(out_dir, out_file_name + 'summary.txt'), 'w')
    missing_mods_df = pd.DataFrame(columns=cols_dataframe)
    row_to_append_df = pd.DataFrame(columns=cols_dataframe)
    subjects_paths_lists = glob(path.join(in_dir, '*sub-*'))
    subjects_paths_lists.sort()

    if len(subjects_paths_lists) == 0:
        raise "No subjects found or dataset not BIDS complaint."
    # Check the modalities available for each session
    for ses in sessions_found:
        mods_avail_bids = []
        for sub_path in subjects_paths_lists:
            mods_avail_bids = []
            subj_id = sub_path.split(os.sep)[-1]
            row_to_append_df['participant_id'] = pd.Series(subj_id)
            ses_path_avail = glob(path.join(sub_path, ses))
            if len(ses_path_avail)==0:
                mmt.increase_missing_ses(ses)
                for mod in mods_avail:
                    row_to_append_df[mod] = pd.Series('0')
            else:

                ses_path = ses_path_avail[0]
                mods_paths_folders = glob(path.join(ses_path, '*/'))

                for p in mods_paths_folders:
                    p = p[:-1]
                    mods_avail_bids.append(p.split('/').pop())

                # Check if a modality folder is available and if is empty
                if 'func' in mods_avail_bids:
                    # Extract all the task available
                    for m in mods_avail_dict['func']:
                        tokens = m.split('_')
                        task_name = tokens[1]
                        task_avaL_list = glob(path.join(ses_path,'func' ,'*'+task_name+'*'))

                        if len(task_avaL_list) == 0:
                            row_to_append_df[m] = pd.Series('0')
                        else:
                            row_to_append_df[m] = pd.Series('1')
                # If the folder is not available but the modality is in the list of the available one mark it as missing
                else:
                    if 'func' in mods_avail_dict:
                        for m in mods_avail_dict['func']:
                            row_to_append_df[m] = pd.Series('0')
                        mmt.add_missing_mod(ses, m)

                if 'dwi' in mods_avail_bids:
                    row_to_append_df['dwi'] = pd.Series('1')
                else:
                    if 'dwi' in mods_avail:
                        row_to_append_df['dwi'] = pd.Series('0')
                        mmt.add_missing_mod(ses, 'dwi')

                if 'anat' in mods_avail_bids:
                    for m in mods_avail_dict['anat']:
                        anat_aval_list = glob(path.join(ses_path, 'anat', '*.nii.gz'))
                        if len(anat_aval_list) > 0:
                            row_to_append_df[m] = pd.Series('1')
                        else:
                            row_to_append_df[m] = pd.Series('0')
                            mmt.add_missing_mod(ses, m)
                else:
                    if mods_avail_dict.has_key('anat'):
                        for m in mods_avail_dict['anat']:
                            row_to_append_df[m] = pd.Series('0')
                            mmt.add_missing_mod(ses, m)

                if 'fmap' in mods_avail_bids:
                    row_to_append_df['fmap'] = pd.Series('1')
                else:
                    if 'fmap' in mods_avail:
                        row_to_append_df['fmap'] = pd.Series('0')
                        mmt.add_missing_mod(ses, 'fmap')
                if 'pet' in mods_avail_bids:
                    row_to_append_df['pet'] = pd.Series('1')
                else:
                    if 'pet' in mods_avail:
                        row_to_append_df['pet'] = pd.Series('0')
                        mmt.add_missing_mod(ses, 'pet')


            missing_mods_df = missing_mods_df.append(row_to_append_df)
            row_to_append_df = pd.DataFrame(columns=cols_dataframe)

        missing_mods_df = missing_mods_df[cols_dataframe]
        missing_mods_df.to_csv(path.join(out_dir, out_file_name+ses+'.tsv'), sep='\t', index=False)
        missing_mods_df = pd.DataFrame(columns=cols_dataframe)

    print_statistics(summary_file, len(subjects_paths_lists), sessions_found, mmt )


def create_subs_sess_list(dataset_path, out_dir, file_name = ''):
    """
    Create the file subject_session_list.tsv that contains the list of the visits for each subject for a BIDS compliant
    dataset.


    :param dataset_path: path to the BIDS directory
    :param out_dir: path to the output directory
    :param file_name: name of the output file

    """

    if file_name == '':
        file_name = 'subjects_sessions_list.tsv'

    subjs_sess_tsv = open(path.join(out_dir, file_name), 'w')
    subjs_sess_tsv.write('participant_id' + '\t' + 'session_id' + '\n')
    subjects_paths = glob(path.join(dataset_path, '*sub-*'))
    # Sort the subjects list
    subjects_paths.sort()

    if len(subjects_paths) == 0:
        raise Exception('Dataset empty or not BIDS-compliant.')

    for sub_path in subjects_paths:
        subj_id = sub_path.split(os.sep)[-1]
        sess_list = glob(path.join(sub_path, '*ses-*'))

        for ses_path in sess_list:
            session_name = ses_path.split(os.sep)[-1]
            subjs_sess_tsv.write(subj_id+'\t'+session_name+'\n')

    subjs_sess_tsv.close()




