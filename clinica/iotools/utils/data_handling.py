# coding: utf8

"""
Data handling scripts
"""
__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Completed"


def create_merge_file(bids_dir, out_tsv, caps_dir=None, tsv_file=None, pipelines=None, **kwargs):
    """
    Merge all the .TSV files containing clinical data of a BIDS compliant dataset and store
    the result inside a .TSV file.

    Args:
        bids_dir: path to the BIDS folder
        out_tsv: path to the output tsv file
        caps_dir: path to the CAPS folder (optional)
        tsv_file: TSV file containing the subjects with their sessions (optional)
        pipelines: when adding CAPS information, indicates the pipelines that will be merged (optional)

    """
    from os import path
    from glob import glob
    import os
    import pandas as pd
    import numpy as np
    import warnings
    from .pipeline_handling import InitException, DatasetError
    from ...pipelines.engine import get_subject_session_list

    if caps_dir is not None:
        if not path.isdir(caps_dir):
            raise IOError('The path to the CAPS directory is wrong')

    col_list = []
    scans_dict = {}

    if not os.path.isfile(path.join(bids_dir, 'participants.tsv')):
        raise IOError('participants.tsv not found in the specified BIDS directory')
    participants_df = pd.read_csv(path.join(bids_dir, 'participants.tsv'), sep='\t')

    sessions, subjects = get_subject_session_list(bids_dir, ss_file=tsv_file, use_session_tsv=True)
    n_sessions = len(sessions)

    # Find what is dir and what is file_name
    if os.sep not in out_tsv:
        out_dir = os.getcwd()
        if out_tsv == '.':
            out_file_name = 'merge_tsv.tsv'
        else:
            out_file_name = out_tsv
    else:
        out_file_name = out_tsv.split(os.sep)[-1]
        out_dir = path.dirname(out_tsv)

    if len(out_file_name) == 0:
        out_file_name = 'merge_tsv.tsv'

    if '.' not in out_file_name:
        out_file_name = out_file_name + '.tsv'
    else:
        extension = os.path.splitext(out_file_name)[1]
        if extension != '.tsv':
            raise TypeError('Output file must be .tsv.')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for col in participants_df.columns.values:
        col_list.append(col)

    merged_df = pd.DataFrame(columns=col_list)

    # BIDS part
    i_subject = 0
    while i_subject < n_sessions:
        sub_path = path.join(bids_dir, subjects[i_subject])
        sub_name = sub_path.split(os.sep)[-1]
        # For each subject, extract the relative row from the dataframe
        row_participant = participants_df[participants_df['participant_id'] == sub_name]
        # Open the sessions file related to the subject
        sessions_df = pd.read_csv(path.join(sub_path, sub_name + '_sessions.tsv'), sep='\t')
        # Looking for the sessions corresponding to the subject
        loc_sessions = []
        i_session = i_subject
        while i_session < n_sessions and subjects[i_session] == subjects[i_subject]:
            loc_sessions.append(i_session)
            i_session += 1

        # For each session found
        # extract the information contained in the scans files
        # for line in range(0, len(sessions_df)):
        for i_session in loc_sessions:
            # Extract and convert to a dictionary information
            # regarding the session
            row_session_df = sessions_df[sessions_df.session_id == sessions[i_session]]
            row_session_df.reset_index(inplace=True, drop=True)
            if len(row_session_df) == 0:
                raise DatasetError(sessions_df.loc[0, 'session_id'] + ' / ' + sessions[i_session])

            new_cols = [s for s in row_session_df.columns.values if s not in col_list]
            if len(new_cols) != 0:
                for i in range(0, len(new_cols)):
                    col_list.append(new_cols[i])

            session_id = row_session_df.loc[0, 'session_id']
            if os.path.isfile(path.join(bids_dir, sub_name, 'ses-' + session_id,
                                        sub_name + '_' + 'ses-' + session_id + '_scans.tsv')):
                scans_df = pd.read_csv(path.join(bids_dir, sub_name, 'ses-' + session_id,
                                                 sub_name + '_' + 'ses-' + session_id + '_scans.tsv'), sep='\t')
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
                            new_col_name = col + '_' + mod_type
                            scans_dict.update({new_col_name: value})
                row_scans = pd.DataFrame(scans_dict, index=[0])
            else:
                row_scans = pd.DataFrame()

            new_cols = [s for s in row_scans.columns.values if s not in col_list]
            if len(new_cols) != 0:
                for i in range(0, len(new_cols)):
                    col_list.append(new_cols[i])

            row_to_append_df = pd.DataFrame(columns=row_participant.columns)
            for col in row_participant:
                row_to_append_df[col] = row_participant[col]

            # Append all the data inside session_df
            for col in row_session_df:
                row_to_append_df[col] = row_session_df[col].values[0]

            for col in row_scans:
                row_to_append_df[col] = row_scans[col].values[0]

            merged_df = merged_df.append(row_to_append_df)
        scans_dict = {}
        i_subject = loc_sessions[-1] + 1

    old_index = col_list.index('session_id')
    col_list.insert(1, col_list.pop(old_index))
    merged_df = merged_df[col_list]
    merged_df.to_csv(path.join(out_dir, out_file_name), sep='\t', index=False)

    len_BIDS = len(merged_df.columns)

    # # Call the script for computing the missing modalities
    # # and append the result to the merged file
    # compute_missing_mods(bids_dir, out_dir, 'tmpG7VIY0')
    # tmp_ses = glob(path.join(out_dir, 'tmpG7VIY0*'))
    # for f in tmp_ses:
    #     # Skip the summary file
    #     if 'summary' not in f:
    #         # Load the file
    #         mss_df = pd.read_csv(f, sep='\t')
    #         f_name = f.split(os.sep)[-1]
    #         patterns = f_name.split('-')
    #         ses_id = patterns[len(patterns) - 1]
    #         ses_id = ses_id.replace('.tsv', '')
    #         cols = mss_df.columns.values
    #
    #         # If the file opened contains new columns,
    #         # add them to the existing merged_df
    #         for col_name in cols:
    #             if col_name not in col_list:
    #                 merged_df[col_name] = 0
    #
    #         for i in range(0, len(mss_df)):
    #             row = mss_df.iloc[i]
    #             subj_idx = merged_df[(merged_df['participant_id'] == row['participant_id']) & (
    #                     merged_df['session_id'] == ses_id)].index.tolist()
    #
    #             if len(subj_idx) > 1:
    #                 raise ValueError('Multiple row for the same visit in the merge-tsv file.')
    #             elif len(subj_idx) == 0:
    #                 print 'Warning: Found modalities missing information for the subject:' + row[
    #                     'participant_id'] + ' visit:' + ses_id + ' but the subject is not included in the column participant_id.'
    #                 continue
    #             else:
    #                 subj_idx = subj_idx[0]
    #             for col_name in cols:
    #                 if not col_name == 'participant_id':
    #                     merged_df.iloc[subj_idx, merged_df.columns.get_loc(col_name)] = row[col_name]
    #
    # # Remove all the temporary files created
    # for f in tmp_ses:
    #     os.remove(f)
    #
    # if true_false_mode:
    #     merged_df = merged_df.replace(['Y', 'N'], ['0', '1'])

    merged_df = merged_df.reset_index(drop=True)

    # CAPS
    if caps_dir is not None:
        # Call the different pipelines
        from .pipeline_handling import t1_volume_pipeline, pet_volume_pipeline

        pipeline_options = {
            't1-volume': t1_volume_pipeline,
            'pet-volume': pet_volume_pipeline
        }
        columns_summary = ['pipeline_name', 'group_id', 'atlas_id', 'regions_number', 'first_column_name', 'last_column_name']
        merged_summary_df = pd.DataFrame(columns=columns_summary)
        if pipelines is None:
            for key, pipeline in pipeline_options.items():
                try:
                    merged_df, summary_df = pipeline(caps_dir, merged_df, **kwargs)
                    merged_summary_df = pd.concat([merged_summary_df, summary_df])

                except InitException:
                    warnings.warn('This pipeline was not initialized: ' + key)
        else:
            for pipeline in pipelines:
                merged_df, summary_df = pipeline_options[pipeline](caps_dir, merged_df, **kwargs)
                merged_summary_df = pd.concat([merged_summary_df, summary_df])

        n_atlas = len(merged_summary_df)
        index_column_df = pd.DataFrame(index=np.arange(n_atlas), columns=['first_column_index', 'last_column_index'])
        index_column_df.iat[0, 0] = len_BIDS
        index_column_df.iat[n_atlas - 1, 1] = np.shape(merged_df)[1] - 1
        for i in range(1, n_atlas):
            index_column_df.iat[i, 0] = index_column_df.iat[i-1, 0] + merged_summary_df.iat[i-1, 3]
            index_column_df.iat[i-1, 1] = index_column_df.iat[i, 0] - 1

        merged_summary_df.reset_index(inplace=True, drop=True)
        merged_summary_df = pd.concat([merged_summary_df, index_column_df], axis=1)
        summary_filename = out_file_name.split('.')[0] + '_summary.tsv'
        merged_summary_df.to_csv(path.join(out_dir, summary_filename), sep='\t', index=False)

    merged_df.to_csv(path.join(out_dir, out_file_name), sep='\t', index=False)


def find_mods_and_sess(bids_dir):
    """
    Find all the modalities and sessions available for a given BIDS dataset

    Args:
        bids_dir: path to the BIDS dataset

    Returns:
        mods_dict: a dictionary that stores the sessions and modalities found and has the following structure.
    Example:
    {
        'sessions': ['ses-M00', 'ses-M18'],
        'fmap': ['fmap'],
        'anat': ['flair', 't1w'],
        'func': ['func_task-rest'],
        'dwi': ['dwi']
    }

    """
    from glob import glob
    from os import path
    import os

    mods_dict = {}
    mods_list = []
    subjects_paths_lists = glob(path.join(bids_dir, '*sub-*'))

    for sub_path in subjects_paths_lists:
        ses_paths = glob(path.join(sub_path, '*ses-*'))
        for session in ses_paths:
            ses_name = session.split(os.sep)[-1]
            mods_avail = []
            if 'sessions' in mods_dict:
                if ses_name not in mods_dict['sessions']:
                    mods_dict['sessions'].append(ses_name)
            else:
                mods_dict.update({'sessions': [ses_name]})
            mods_paths_folders = glob(path.join(session, '*/'))

            for p in mods_paths_folders:
                p = p[:-1]
                mods_avail.append(p.split('/').pop())

            if 'func' in mods_avail:
                list_funcs_paths = glob(path.join(session, 'func', '*bold.nii.gz'))
                for func_path in list_funcs_paths:
                    func_name = func_path.split(os.sep)[-1]
                    func_name_tokens = func_name.split('_')
                    func_task = func_name_tokens[2]
                if 'func' in mods_dict:
                    if 'func_' + func_task not in mods_dict['func']:
                        mods_dict['func'].append('func_' + func_task)
                else:
                    mods_dict.update({'func': ['func_' + func_task]})

                if 'func_' + func_task not in mods_list:
                    mods_list.append('func_' + func_task)

            if 'dwi' in mods_avail:
                if 'dwi' not in mods_dict:
                    mods_dict.update({'dwi': ['dwi']})
                if 'dwi' not in mods_list:
                    mods_list.append('dwi')

            if 'fmap' in mods_avail:
                if 'fmap' not in mods_dict:
                    mods_dict.update({'fmap': ['fmap']})
                if 'fmap' not in mods_list:
                    mods_list.append('fmap')

            if 'pet' in mods_avail:
                if 'pet' not in mods_dict:
                    mods_dict.update({'pet': ['pet']})
                if 'pet' not in mods_list:
                    mods_list.append('pet')

            if 'anat' in mods_avail:
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
                        if 'anat' in mods_dict:
                            if anat_type not in mods_dict['anat']:
                                anat_aval = mods_dict['anat']
                                anat_aval.append(anat_type)
                                mods_dict.update({'anat': anat_aval})
                        else:
                            mods_dict.update({'anat': [anat_type]})

                        if anat_type not in mods_list:
                            mods_list.append(anat_type)

    return mods_dict


def compute_missing_mods(bids_dir, out_dir, output_prefix=''):
    """
    Compute the list of missing modalities for each subject in a BIDS compliant dataset

    Args:
        bids_dir: path to the BIDS directory
        out_dir: path to the output folder
        output_prefix: string that replace the default prefix ('missing_mods_') in the name of all the output files
    created
    """
    from ..converter_utils import MissingModsTracker, print_statistics
    import os
    from os import path
    import pandas as pd
    from glob import glob

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Find all the modalities and sessions available for the input dataset
    mods_and_sess = find_mods_and_sess(bids_dir)
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
    subjects_paths_lists = glob(path.join(bids_dir, '*sub-*'))
    subjects_paths_lists.sort()

    if len(subjects_paths_lists) == 0:
        raise IOError("No subjects found or dataset not BIDS complaint.")
    # Check the modalities available for each session
    for ses in sessions_found:
        mods_avail_bids = []
        for sub_path in subjects_paths_lists:
            mods_avail_bids = []
            subj_id = sub_path.split(os.sep)[-1]
            row_to_append_df['participant_id'] = pd.Series(subj_id)
            ses_path_avail = glob(path.join(sub_path, ses))
            if len(ses_path_avail) == 0:
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
                        task_avail_list = glob(path.join(
                            ses_path, 'func', '*' + task_name + '*')
                        )

                        if len(task_avail_list) == 0:
                            row_to_append_df[m] = pd.Series('0')
                        else:
                            row_to_append_df[m] = pd.Series('1')
                # If the folder is not available but the modality is
                # in the list of the available one mark it as missing
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
                    if 'anat' in mods_avail_dict:
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
        missing_mods_df.to_csv(path.join(out_dir, out_file_name + ses + '.tsv'), sep='\t', index=False,
                               encoding='utf-8')
        missing_mods_df = pd.DataFrame(columns=cols_dataframe)

    print_statistics(summary_file, len(subjects_paths_lists), sessions_found, mmt)


def create_subs_sess_list(input_dir, output_dir,
                          file_name=None, is_bids_dir=True, use_session_tsv=False):
    """
    Create the file subject_session_list.tsv that contains the list
    of the visits for each subject for a BIDS or CAPS compliant dataset.

    Args:
        input_dir (str): Path to the BIDS or CAPS directory.
        output_dir (str): Path to the output directory
        file_name: name of the output file
        is_bids_dir (boolean): Specify if input_dir is a BIDS directory or
            not (i.e. a CAPS directory)
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files
    """
    from os import path
    from glob import glob
    import pandas as pd
    import os

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if file_name is None:
        file_name = 'subjects_sessions_list.tsv'
    subjs_sess_tsv = open(path.join(output_dir, file_name), 'w')
    subjs_sess_tsv.write('participant_id' + '\t' + 'session_id' + '\n')

    if is_bids_dir:
        path_to_search = input_dir
    else:
        path_to_search = path.join(input_dir, 'subjects')
    subjects_paths = glob(path.join(path_to_search, '*sub-*'))

    # Sort the subjects list
    subjects_paths.sort()

    if len(subjects_paths) == 0:
        raise IOError('Dataset empty or not BIDS/CAPS compliant.')

    for sub_path in subjects_paths:
        subj_id = sub_path.split(os.sep)[-1]

        if use_session_tsv:
            session_df = pd.read_csv(path.join(sub_path, subj_id + '_sessions.tsv'), sep='\t')
            session_list = list(session_df['session_id'].values)
            for session in session_list:
                subjs_sess_tsv.write(subj_id + '\t' + session + '\n')

        else:
            sess_list = glob(path.join(sub_path, '*ses-*'))

            for ses_path in sess_list:
                session_name = ses_path.split(os.sep)[-1]
                subjs_sess_tsv.write(subj_id + '\t' + session_name + '\n')

    subjs_sess_tsv.close()
