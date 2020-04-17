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
    from clinica.utils.participant import get_subject_session_list

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
                        anat_aval_list = [elem for elem in anat_aval_list if m.lower() in elem.lower()]
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
            session_list = list(session_df['session_id'].to_numpy())
            for session in session_list:
                subjs_sess_tsv.write(subj_id + '\t' + session + '\n')

        else:
            sess_list = glob(path.join(sub_path, '*ses-*'))

            for ses_path in sess_list:
                session_name = ses_path.split(os.sep)[-1]
                subjs_sess_tsv.write(subj_id + '\t' + session_name + '\n')

    subjs_sess_tsv.close()


def center_nifti_origin(input_image, output_image):

    """

    Put the origin of the coordinate system at the center of the image

    Args:
        input_image: path to the input image
        output_image: path to the output image (where the result will be stored)

    Returns:
        path of the output image created
    """

    import nibabel as nib
    import numpy as np
    from colorama import Fore
    from nibabel.spatialimages import ImageFileError
    from os.path import isfile
    import os

    error_str = None
    try:
        img = nib.load(input_image)
    except FileNotFoundError:
        error_str = Fore.RED + '[Error] No such file ' + input_image + Fore.RESET
    except ImageFileError:
        error_str = Fore.RED + '[Error] File ' + input_image + ' could not be read.' + Fore.RESET
    except Exception as e:
        error_str = Fore.RED + '[Error] File ' + input_image + ' could not be loaded with nibabel: ' + \
                    str(e) + Fore.RESET

    if not error_str:
        try:
            canonical_img = nib.as_closest_canonical(img)
            hd = canonical_img.header

            qform = np.zeros((4, 4))
            for i in range(1, 4):
                qform[i - 1, i - 1] = hd['pixdim'][i]
                qform[i - 1, 3] = -1.0 * hd['pixdim'][i] * hd['dim'][i] / 2.0
            new_img = nib.Nifti1Image(canonical_img.get_data(caching='unchanged'), affine=qform, header=hd)

            # Without deleting already-existing file, nib.save causes a severe bug on Linux system
            if isfile(output_image):
                os.remove(output_image)

            nib.save(new_img, output_image)
            if not isfile(output_image):
                error_str = Fore.RED + '[Error] NIfTI file created but Clinica could not save it to ' \
                            + output_image + '. Please check that the output folder has the correct permissions.' \
                            + Fore.RESET
        except Exception as e:
            error_str = Fore.RED + '[Error] File ' + input_image + ' could not be processed with nibabel: ' + \
                        str(e) + Fore.RESET

    return output_image, error_str


def center_all_nifti(bids_dir, output_dir, modality, center_all_files=False):
    """
    Center all the NIfTI images of the input BIDS folder into the empty output_dir specified in argument.
    All the files from bids_dir are copied into output_dir, then all the NIfTI images we can found are replaced by their
    centered version if their center if off the origin by more than 50 mm.

    Args:
        bids_dir: (str) path to bids directory
        output_dir: (str) path to EMPTY output directory
        modality: (list of str) modalities to convert
        center_all_files: (bool) center only files that may cause problem for SPM if false. If true, center all NIfTI

    Returns:
        List of the centered files
    """
    from colorama import Fore
    from clinica.utils.inputs import check_bids_folder
    from clinica.utils.exceptions import ClinicaBIDSError
    from os.path import join, basename
    from glob import glob
    from os import listdir
    from os.path import isdir, isfile
    from shutil import copy2, copytree

    # output and input must be different, so that we do not mess with user's data
    if bids_dir == output_dir:
        raise ClinicaBIDSError(Fore.RED + '[Error] Input BIDS and output directories must be different' + Fore.RESET)

    assert isinstance(modality, list), 'modality arg must be a list of str'

    # check that input is a BIDS dir
    check_bids_folder(bids_dir)

    for f in listdir(bids_dir):
        if isdir(join(bids_dir, f)) and not isdir(join(output_dir, f)):
            copytree(join(bids_dir, f), join(output_dir, f))
        elif isfile(join(bids_dir, f)) and not isfile(join(output_dir, f)):
            copy2(join(bids_dir, f), output_dir)

    pattern = join(output_dir, '**/*.nii*')
    nifti_files = glob(pattern, recursive=True)

    # Now filter this list by elements in modality list
    #   For each file:
    #       if any modality name (lowercase) is found in the basename of the file:
    #           keep the file
    nifti_files_filtered = [f for f in nifti_files
                            if any(elem.lower() in basename(f).lower() for elem in modality)]

    # Remove those who are centered
    if not center_all_files:
        nifti_files_filtered = [file for file in nifti_files_filtered if not is_centered(file)]

    all_errors = []
    for f in nifti_files_filtered:
        print('Handling ' + f)
        _, current_error = center_nifti_origin(f, f)
        if current_error:
            all_errors.append(current_error)
    if len(all_errors) > 0:
        final_error_msg = Fore.RED + '[Error] Clinica encoutered ' + str(len(all_errors)) \
                          + ' error(s) while trying to center all NIfTI images.\n'
        for error in all_errors:
            final_error_msg += '\n' + error
        raise RuntimeError(final_error_msg)
    return nifti_files_filtered


def are_far_appart(file1, file2, threshold=80):
    """
    Tells if 2 files have a center located at more than a threshold distance
    Args:
        file1: (str) path to the first nifti file
        file2: (str) path to the second nifti file
        threshold: threshold to consider whether 2 files are too far appart

    Returns:
        True if distance between `file1` and `file2` is greter than `threshold`, False otherwise.
    """
    from os.path import isfile
    import numpy as np

    assert isfile(file1)
    assert isfile(file2)

    center1 = get_world_coordinate_of_center(file1)
    center2 = get_world_coordinate_of_center(file2)

    return np.linalg.norm(center2 - center1, ord=2) > threshold


def write_list_of_files(file_list, output_file):
    """
    Save `file_list` list of files into `output_file` text file.
    Args:
        file_list: (list of str) of path to files
        output_file: (str) path to the output txt file

    Returns:
        output_file
    """
    from os.path import isfile

    assert isinstance(file_list, list), 'First argument must be a list'
    assert isinstance(output_file, str), 'Second argument must be a str'
    if isfile(output_file):
        return None

    text_file = open(output_file, 'w+')
    for created_file in file_list:
        text_file.write(created_file + '\n')
    text_file.close()
    return output_file


def check_relative_volume_location_in_world_coordinate_system(label_1, nifti_list1,
                                                              label_2, nifti_list2,
                                                              bids_dir,
                                                              modality):
    """
    Check if the NIfTI file list nifti_list1 and nifti_list2 provided in argument are not too far apart (otherwise coreg
    in SPM may fail. Norm between center of volumes of 2 files must be less than 80 mm.

    Args:
        label_1: label of the first nifti_list1 files (used in potential warning message)
        nifti_list1: first set of files
        label_2: label of the second nifti_list
        nifti_list2: second set of files, must be same length as nifti_list1
        bids_dir: bids directory (used in potential warning message)
        modality: string that must be used in argument of: clinica iotools bids --modality MODALITY (used in potential
                warning message)
    Returns:
        Nothing
    """
    import numpy as np
    from colorama import Fore
    from os.path import abspath, basename
    from clinica.utils.stream import cprint
    import sys

    center_coordinate_1 = [get_world_coordinate_of_center(file) for file in nifti_list1]
    center_coordinate_2 = [get_world_coordinate_of_center(file) for file in nifti_list2]

    l2_norm = [np.linalg.norm(center_1 - center_2) for center_1, center_2 in zip(center_coordinate_1, center_coordinate_2)]
    pairs_with_problems = [i for i, norm in enumerate(l2_norm) if norm > 80]

    if len(pairs_with_problems) > 0:
        warning_message = (Fore.YELLOW + '[Warning] It appears that ' + str(len(pairs_with_problems)) + ' pairs of files'
                           + ' have an important relative offset. SPM coregistration has a high probability to fail '
                           + 'on these files:\n\n')

        # File column width : 3 spaces more than the longest string to display
        file1_width = max(3 + len(label_1),
                          3 + max(len(basename(file)) for file in [nifti_list1[k] for k in pairs_with_problems]))
        file2_width = max(3 + len(label_2),
                          3 + max(len(basename(file)) for file in [nifti_list2[k] for k in pairs_with_problems]))

        norm_width = len('Relative distance')

        warning_message += ('%-' + str(file1_width)
                            + 's%-' + str(file2_width)
                            + 's%-' + str(norm_width) + 's') % (label_1,
                                                                label_2,
                                                                'Relative distance')

        warning_message += '\n' + '-' * (file1_width + file2_width + norm_width) + '\n'
        for file1, file2, norm in zip([nifti_list1[k] for k in pairs_with_problems],
                                      [nifti_list2[k] for k in pairs_with_problems],
                                      [l2_norm[k] for k in pairs_with_problems]):
            # Nice formatting as array
            # % escape character
            # - aligned to the left, with the size of the column
            # s = string, f = float
            # . for precision with float
            # https://docs.python.org/2/library/stdtypes.html#string-formatting for more information
            warning_message += ('%-' + str(file1_width)
                                + 's%-' + str(file2_width)
                                + 's%-' + str(norm_width) + '.2f\n') % (str(basename(file1)),
                                                                        str(basename(file2)),
                                                                        norm)
        warning_message += '\nClinica provides a tool to counter this problem by replacing the center of the volume' \
                           + ' at the origin of the world coordinates.\nUse the following command line to correct the '\
                           + 'header of the faulty NIFTI volumes in a new folder:\n' + Fore.RESET \
                           + Fore.BLUE + '\nclinica iotools center-nifti ' + abspath(bids_dir) + ' ' \
                           + abspath(bids_dir) + '_centered --modality "' + modality + '"\n\n'  \
                           + Fore.YELLOW + 'You will find more information on the command by typing ' + Fore.BLUE \
                           + 'clinica iotools center-nifti' + Fore.YELLOW + ' in the console.\nDo you still want to ' \
                           + 'launch the pipeline now?' + Fore.RESET
        cprint(warning_message)
        while True:
            cprint('Your answer [yes/no]:')
            answer = input()
            if answer.lower() in ['yes', 'no']:
                break
            else:
                cprint(Fore.RED + 'You must answer yes or no' + Fore.RESET)
        if answer.lower() == 'no':
            cprint(Fore.RED + 'Clinica will now exit...' + Fore.RESET)
            sys.exit(0)


def check_volume_location_in_world_coordinate_system(nifti_list, bids_dir, modality='t1w'):
    """
    Check if the NIfTI file list nifti_list provided in argument are aproximately centered around the origin of the
    world coordinates. (Problem may arise with SPM segmentation

    If yes, we warn the user of this problem, and propose him to exit clinica in order for him to run:
        clinica iotools center-nifti ...
    or to continue with the execution of the pipeline

    Args:
        nifti_list: (list of str) list of path to nifti files
        bids_dir: (str) path to bids directory associated with this check (in order to propose directly the good
            command line for center-nifti tool)
        modality: (str) to propose directly the good command line option

    Returns:
        Nothing
    """
    from colorama import Fore
    from clinica.utils.stream import cprint
    from os.path import abspath, basename
    import numpy as np
    import sys

    list_non_centered_files = [file for file in nifti_list if not is_centered(file)]
    if len(list_non_centered_files) > 0:
        centers = [get_world_coordinate_of_center(file) for file in list_non_centered_files]
        l2_norm = [np.linalg.norm(center, ord=2) for center in centers]

        # File column width : 3 spaces more than the longest string to display
        file_width = 3 + max(len(basename(file)) for file in list_non_centered_files)
        # Center column width (with a fixed minimum size) : 3 spaces more than the longest string to display
        center_width = max(len('Coordinate of center') + 3,
                           3 + max(len(str(center)) for center in centers))

        warning_message = (Fore.YELLOW + '[Warning] It appears that ' + str(len(list_non_centered_files)) + ' files '
                           + 'have a center way out of the origin of the world coordinate system. SPM has a high prob'
                           + 'ability to fail on these files (for coregistration or segmentation):\n\n')
        warning_message += ('%-' + str(file_width) + 's%-' + str(center_width) + 's%-s') % ('File',
                                                                                            'Coordinate of center',
                                                                                            'Distance to origin')
        # 18 is the length of the string 'Distance to origin'
        warning_message += '\n' + '-' * (file_width + center_width + 18) + '\n'
        for file, center, l2 in zip(list_non_centered_files, centers, l2_norm):
            # Nice formatting as array
            # % escape character
            # - aligned to the left, with the size of the column
            # s = string, f = float
            # . for precision with float
            # https://docs.python.org/2/library/stdtypes.html#string-formatting for more information
            warning_message += ('%-' + str(file_width) + 's%-' + str(center_width) + 's%-25.2f\n') % (basename(file),
                                                                                                      str(center),
                                                                                                      l2)

        cmd_line = (Fore.BLUE
                    + '\nclinica iotools center-nifti '
                    + abspath(bids_dir) + ' '
                    + abspath(bids_dir) + '_centered'
                    + '--modality "' + modality + '"'
                    + '\n\n' + Fore.YELLOW)

        warning_message += ('\nIf you are trying to launch the t1-freesurfer pipeline, you can ignore this message '
                            + 'if you do not want to run the pet-surface pipeline afterward.')

        warning_message += '\nClinica provides a tool to counter this problem by replacing the center of the volume' \
                           + ' at the origin of the world coordinates.\nUse the following command line to correct the '\
                           + 'header of the faulty NIFTI volumes in a new folder:\n'\
                           + cmd_line \
                           + 'You will find more information on the command by typing ' + Fore.BLUE \
                           + 'clinica iotools center-nifti' + Fore.YELLOW + ' in the console.\nDo you still want to '\
                           + 'launch the pipeline now?' + Fore.RESET
        cprint(warning_message)
        while True:
            cprint('Your answer [yes/no]:')
            answer = input()
            if answer.lower() in ['yes', 'no']:
                break
            else:
                cprint(Fore.RED + 'You must answer yes or no' + Fore.RESET)
        if answer.lower() == 'no':
            cprint(Fore.RED + 'Clinica will now exit...' + Fore.RESET)
            sys.exit(0)


def is_centered(nii_volume, threshold_l2=50):
    """
    Tells if a NIfTI volume is centered on the origin of the world coordinate system.

    SPM has troubles to segment files if the center of the volume is not close from the origin of the world coordinate
    system. A series of experiment have been conducted: we take a volume whose center is on the origin of the world
    coordinate system. We add an offset using coordinates of affine matrix [0, 3], [1, 3], [2, 3] (or by modifying the
    header['srow_x'][3], header['srow_y'][3], header['srow_z'][3], this is strictly equivalent).

    It has been determined that volumes were still segmented with SPM when the L2 distance between origin and center of
    the volume did not exceed 100 mm. Above this distance, either the volume is either not segmented (SPM error), or the
    produced segmentation is wrong (not the shape of a brain anymore)

    Args:
        nii_volume: path to NIfTI volume
        threshold_l2: maximum distance between origin of the world coordinate system and the center of the volume to
                    be considered centered. The threshold were SPM segmentation stops working is around 100 mm
                    (it was determined empirically after several trials on a genrated dataset), so default value is 50
                    mm in order to have a security margin, even when dealing with coregistred files afterward)

    Returns:
        True or False

    """
    import numpy as np
    from os.path import basename

    center = get_world_coordinate_of_center(nii_volume)

    # Compare to the threshold and retun boolean
    # if center is a np.nan, comparison will be False, and False will be returned
    distance_from_origin = np.linalg.norm(center, ord=2)
    # if not np.isnan(distance_from_origin):
    #     print('\t' + basename(nii_volume) + ' has its center at {0:.2f} mm of the origin.'.format(distance_from_origin))
    if distance_from_origin < threshold_l2:
        return True
    else:
        # If center is a np.nan,
        return False


def get_world_coordinate_of_center(nii_volume):
    """
    Extract the world coordinates of the center of the image. Based on methods described
    here : https://brainder.org/2012/09/23/the-nifti-file-format/

    Args:
        nii_volume: path to nii volume

    Returns:

    """
    from os.path import isfile
    import nibabel as nib
    from colorama import Fore
    import numpy as np

    assert isinstance(nii_volume, str), 'input argument nii_volume must be a str'
    assert isfile(nii_volume), 'input argument must be a path to a file'

    try:
        orig_nifti = nib.load(nii_volume)
    except nib.filebasedimages.ImageFileError:
        print(Fore.RED + '[Error] ' + nii_volume
              + ' could not be read by nibabel. Is it a valid NIfTI file ?' + Fore.RESET)
        return np.nan

    head = orig_nifti.header

    if isinstance(head, nib.freesurfer.mghformat.MGHHeader):
        # If MGH volume
        center_coordinates_world = vox_to_world_space_method_3_bis(head['dims'][0:3] / 2, head)
    else:
        # Standard NIfTI volume
        center_coordinates = get_center_volume(head)

        if head['qform_code'] > 0:
            center_coordinates_world = vox_to_world_space_method_2(center_coordinates, head)
        elif head['sform_code'] > 0:
            center_coordinates_world = vox_to_world_space_method_3(center_coordinates, head)
        elif head['sform_code'] == 0:
            center_coordinates_world = vox_to_world_space_method_1(center_coordinates, head)
        else:
            center_coordinates_world = np.nan
    return center_coordinates_world


def get_center_volume(header):
    """
    Get the voxel coordinates of the center of the data, using header information
    Args:
        header: a nifti header

    Returns:
        Voxel coordinates of the center of the volume
    """
    import numpy as np

    center_x = header['dim'][1] / 2
    center_y = header['dim'][2] / 2
    center_z = header['dim'][3] / 2
    return np.array([center_x,
                     center_y,
                     center_z])


def vox_to_world_space_method_1(coordinates_vol, header):
    """
    The Method 1 is for compatibility with analyze and is not supposed to be used as the main orientation method. But it
    is used if sform_code = 0. The world coordinates are determined simply by scaling by the voxel size by their
    dimension stored in pixdim. More information here: https://brainder.org/2012/09/23/the-nifti-file-format/
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    return np.array(coordinates_vol) * np.array(header['pixdim'][1],
                                                header['pixdim'][2],
                                                header['pixdim'][3])


def vox_to_world_space_method_2(coordinates_vol, header):
    """
    The Method 2 is used when short qform_code is larger than zero. To get the coordinates, we multiply a rotation
    matrix (r_mat) by coordinates_vol, then perform hadamart with pixel dimension pixdim (like in method 1). Then we add
    an offset (qoffset_x, qoffset_y, qoffset_z)

    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    def get_r_matrix(h):
        """
        Get rotation matrix, more information here: https://brainder.org/2012/09/23/the-nifti-file-format/
        Args:
            h: header

        Returns:
            Rotation matrix
        """
        b = h['quatern_b']
        c = h['quatern_c']
        d = h['quatern_d']
        a = np.sqrt(1 - (b ** 2) - (c ** 2) - (d ** 2))
        r = np.zeros((3, 3))
        r[0, 0] = (a ** 2) + (b ** 2) - (c ** 2) - (d ** 2)
        r[0, 1] = 2 * ((b * c) - (a * d))
        r[0, 2] = 2 * ((b * d) + (a * c))
        r[1, 0] = 2 * ((b * c) + (a * d))
        r[1, 1] = (a ** 2) + (c ** 2) - (b ** 2) - (d ** 2)
        r[1, 2] = 2 * ((c * d) - (a * b))
        r[2, 0] = 2 * ((b * d) - (a * c))
        r[2, 1] = 2 * ((b * d) - (a * c))
        r[2, 2] = (a ** 2) + (d ** 2) - (b ** 2) - (c ** 2)
        return r
    i = coordinates_vol[0]
    j = coordinates_vol[1]
    k = coordinates_vol[2]
    if header['qform_code'] > 0:
        r_mat = get_r_matrix(header)
    else:
        # Should never be reached
        raise ValueError('qform_code must be greater than 0 to use this method')
    q = header['pixdim'][0]
    if q not in [-1, 1]:
        print('q was ' + str(q), ', now is 1')
        q = 1
    return np.dot(r_mat, np.array([i, j, q * k])) * np.array(header['pixdim'][1:4]) + np.array([header['qoffset_x'],
                                                                                                header['qoffset_y'],
                                                                                                header['qoffset_z']])


def vox_to_world_space_method_3(coordinates_vol, header):
    """
    This method is used when sform_code is larger than zero. It relies on a full affine matrix, stored in the header in
     the fields srow_[x,y,y], to map voxel to world coordinates.
     When a nifti file is created with raw data and affine=..., this is this method that is used to decypher the
     voxel-to-world correspondance.
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space

    """
    import numpy as np

    def get_aff_matrix(h):
        """
        Get affine transformation matrix, described here : https://brainder.org/2012/09/23/the-nifti-file-format/
        Args:
            h: header

        Returns:
            affine transformation matrix
        """
        mat = np.zeros((4, 4))
        mat[0, 0] = h['srow_x'][0]
        mat[0, 1] = h['srow_x'][1]
        mat[0, 2] = h['srow_x'][2]
        mat[0, 3] = h['srow_x'][3]
        mat[1, 0] = h['srow_y'][0]
        mat[1, 1] = h['srow_y'][1]
        mat[1, 2] = h['srow_y'][2]
        mat[1, 3] = h['srow_y'][3]
        mat[2, 0] = h['srow_z'][0]
        mat[2, 1] = h['srow_z'][1]
        mat[2, 2] = h['srow_z'][2]
        mat[2, 3] = h['srow_z'][3]
        mat[3, 3] = 1
        return mat

    if header['sform_code'] > 0:
        aff = get_aff_matrix(header)
    else:
        # Should never be reached
        raise ValueError('sform_code has a value > 0, so method 3 cannot be used')

    homogeneous_coord = np.concatenate((np.array(coordinates_vol), np.array([1])), axis=0)
    return np.dot(aff, homogeneous_coord)[0:3]


def vox_to_world_space_method_3_bis(coordinates_vol, header):
    """
    This method relies on the same technique as method 3, but for images created by FreeSurfer (MGHImage, MGHHeader)
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: nib.freesurfer.mghformat.MGHHeader object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    affine_trensformation_matrix = header.get_affine()
    homogeneous_coord = np.concatenate((np.array(coordinates_vol), np.array([1])), axis=0)
    return np.dot(affine_trensformation_matrix, homogeneous_coord)[0:3]
