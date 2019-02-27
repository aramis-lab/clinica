# coding: utf8

"""
Methods to find information in the different pipelines of Clinica
"""
__author__ = "Elina Thubeau--Sutre"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Completed"

from glob import glob
import pandas as pd
import numpy as np
from os import path
import os
import warnings


def pet_volume_pipeline(caps_dir, df, **kwargs):
    """
    This method merge the data of the PET-Volume pipeline to the merged file
    containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        atlas_selection: allows to choose the atlas to merge (default = 'all')
        pvc_restriction: gives the restriction on the inclusion or not of the file with the label 'pvc-rbv'
            1       --> only the atlases containing the label will be used
            0       --> the atlases containing the label won't be used
            None    --> all the atlases will be used

    Returns:
         final_df: a DataFrame containing the information of the bids and the pipeline
    """
    from clinica.pipelines.pet_volume.pet_volume_cli import PETVolumeCLI
    instance = PETVolumeCLI()
    pipeline_name = instance.name
    del instance

    # Recovering optional arguments
    if 'pvc_restriction' in kwargs.keys():
        pvc_restriction = kwargs['pvc_restriction']
    else:
        pvc_restriction = None

    if 'atlas_selection' in kwargs.keys():
        atlas_selection = kwargs['atlas_selection']
    else:
        atlas_selection = None

    if 'group_selection' in kwargs.keys():
        group_selection = kwargs['group_selection']
    else:
        group_selection = None

    if pvc_restriction not in [None, 0, 1]:
        raise ValueError('Error in pvc_restriction value: ' + str(pvc_restriction))

    if not path.isdir(caps_dir):
        raise IOError("The CAPS directory doesn't exist")

    pet_path = 'pet' + os.sep + 'preprocessing'
    caps_dir = caps_dir + os.sep + 'subjects'

    # Lists to be sure of the column in case of missing data
    columns = ['pipeline_name', 'group_id', 'atlas_id', 'regions_number', 'first_column_name', 'last_column_name']
    summary_df = pd.DataFrame(columns=columns)
    col_list = []
    group_list = []

    # Looking for a subject with the correct modality to initialize
    sub_paths = glob(path.join(caps_dir, '*sub-*'))
    flag_initial_subject = True
    i = 0
    while flag_initial_subject:
        sub_path = sub_paths[i]
        i += 1
        ses_paths = glob(path.join(sub_path, '*ses-*'))
        for ses_path in ses_paths:
            possible_mods = os.listdir(ses_path)
            if 'pet' in possible_mods:
                flag_initial_subject = False
                ses_init_path = ses_path
            elif i == len(sub_paths):
                raise InitException(pipeline_name)

    # Initialization must take into account the different groups !
    init_path = ses_init_path + os.sep + pet_path
    if group_selection is None:
        group_paths = glob(path.join(init_path, '*group*'))
    else:
        group_paths = []
        for group_selected in group_selection:
            new_group_path = glob(path.join(init_path, '*' + group_selected + '*'))
            group_paths = group_paths + new_group_path
            if len(new_group_path) == 0:
                warnings.warn('The group wanted does not exist for ' + pipeline_name + ': ' + group_selected, UserWarning)
    for group_path in group_paths:
        group_id = group_path.split(os.sep)[-1]
        group_list.append(group_id)
        stats_path = group_path + os.sep + 'atlas_statistics'
        # Selection of the atlas(es)
        if atlas_selection is None:
            atlas_paths = glob(path.join(stats_path, '*statistics*'))
        else:
            atlas_paths = []
            for atlas_selected in atlas_selection:
                new_atlas_path = glob(path.join(stats_path, '*' + atlas_selected + '*'))
                atlas_paths = atlas_paths + new_atlas_path
                if len(new_atlas_path) == 0:
                    warnings.warn('The atlas wanted does not exist for ' + pipeline_name + ': ' + atlas_selected + ' in group ' + group_id, UserWarning)

        if len(atlas_paths) == 0:
                raise InitException(pipeline_name)

        atlas_paths_preserved = list(atlas_paths)
        # Selection depending on the label 'pvc-rbv'
        for i, atlas_path in enumerate(atlas_paths):
            atlas_split = atlas_path.split(os.sep)[-1].split('_')
            if pvc_restriction == 1 and 'pvc-rbv' != atlas_split[6]:
                atlas_paths_preserved.remove(atlas_path)
            elif pvc_restriction == 0 and 'pvc-rbv' == atlas_split[6]:
                atlas_paths_preserved.remove(atlas_path)

        atlas_paths = atlas_paths_preserved
        atlas_paths.sort()
        for atlas_path in atlas_paths:
            row_summary_df = pd.DataFrame(index=np.arange(1), columns=columns)
            atlas_file = atlas_path.split(os.sep)[-1]
            atlas_id = '_'.join(atlas_file.split('_')[2:-1])  # We get rid of the subject and the session in the label

            init_df = pd.read_csv(atlas_path, sep='\t')
            n_regions = len(init_df['label_name'].tolist())
            label_list = [group_id + '_' + atlas_id + '_ROI-' + str(x) for x in range(n_regions)]
            for label in label_list:
                col_list.append(label)

            row_summary = [pipeline_name, group_id, atlas_id, n_regions, label_list[0], label_list[-1]]
            row_summary_df.iloc[0] = row_summary
            summary_df = pd.concat([summary_df, row_summary_df])

    # Initial empty DataFrame for the pipeline
    number_sessions = len(df)
    pipeline_df = pd.DataFrame(index=np.arange(number_sessions), columns=col_list)

    for i in range(number_sessions):
        # Reading the reference (subject + session)
        participant_id = df['participant_id'][i]
        session_id = df['session_id'][i]

        ses_path = caps_dir + os.sep + participant_id + os.sep + session_id
        mod_path = ses_path + os.sep + pet_path
        ses_df = pd.DataFrame(index=np.arange(1), columns=col_list)

        if os.path.exists(mod_path):
            # Looking for groups
            for group in group_list:
                group_path = mod_path + os.sep + group
                if os.path.exists(group_path):
                    # Looking for atlases
                    atlas_list = summary_df[summary_df.group_id == group]['atlas_id'].values
                    for atlas in atlas_list:
                        atlas_path = group_path + os.sep + 'atlas_statistics' + os.sep
                        atlas_path = atlas_path + participant_id + '_' + session_id + '_' + atlas + '_statistics.tsv'
                        if os.path.exists(atlas_path):
                            n_regions = summary_df[(summary_df.group_id == group) &
                                                   (summary_df.atlas_id == atlas)]['regions_number'].values[0]
                            atlas_df = pd.read_csv(atlas_path, sep='\t')
                            label_list = [group + '_' + atlas + '_ROI-' + str(x) for x in
                                          range(n_regions)]
                            ses_df.loc[:, label_list] = atlas_df['mean_scalar'].as_matrix()

        pipeline_df.iloc[i] = ses_df.iloc[0]

    final_df = pd.concat([df, pipeline_df], axis=1)

    return final_df, summary_df


def t1_volume_pipeline(caps_dir, df, **kwargs):
    """
    This method merge the data of the t1-volume pipeline to the merged file containing the
    BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        atlas_selection: allows to choose the atlas to merge. If None all atlases are selected.

    Returns:
        final_df: a DataFrame containing the information of the bids and the pipeline
    """
    from ...pipelines.t1_volume_new_template.t1_volume_new_template_cli import T1VolumeNewTemplateCLI
    instance = T1VolumeNewTemplateCLI()
    pipeline_name = instance.name
    del instance

    # Recovering optional arguments
    if 'atlas_selection' in kwargs.keys():
        atlas_selection = kwargs['atlas_selection']
    else:
        atlas_selection = None

    if 'group_selection' in kwargs.keys():
        group_selection = kwargs['group_selection']
    else:
        group_selection = None

    t1_spm_path = 't1' + os.sep + 'spm' + os.sep + 'dartel'
    caps_dir = caps_dir + os.sep + 'subjects'

    # Lists to be sure of the column in case of missing data
    columns = ['pipeline_name', 'group_id', 'atlas_id', 'regions_number', 'first_column_name', 'last_column_name']
    summary_df = pd.DataFrame(columns=columns)
    col_list = []
    group_list = []

    # Looking for a subject with the correct modality to initialize
    sub_paths = glob(path.join(caps_dir, '*sub-*'))
    flag_initial_subject = True
    i = 0
    while flag_initial_subject:
        sub_path = sub_paths[i]
        i += 1
        ses_paths = glob(path.join(sub_path, '*ses-*'))
        for ses_path in ses_paths:
            possible_mods = os.listdir(ses_path)
            if 't1' in possible_mods:
                mod_path = path.join(ses_path, 't1')
                possible_mods = os.listdir(mod_path)
                if 'spm' in possible_mods:
                    flag_initial_subject = False
                    ses_init_path = ses_path
            elif i == len(sub_paths):
                raise InitException(pipeline_name)

    # Initialization must take into account the different groups !
    init_path = ses_init_path + os.sep + t1_spm_path
    if group_selection is None:
        group_paths = glob(path.join(init_path, '*group*'))
    else:
        group_paths = []
        for group_selected in group_selection:
            new_group_path = glob(path.join(init_path, '*' + group_selected + '*'))
            group_paths = group_paths + new_group_path
            if len(new_group_path) == 0:
                warnings.warn('The group wanted does not exist for ' + pipeline_name + ': ' + group_selected, UserWarning)
    for group_path in group_paths:
        group_id = group_path.split(os.sep)[-1]
        group_list.append(group_id)
        stats_path = group_path + os.sep + 'atlas_statistics'
        # Selection of the atlas(es)
        if atlas_selection is None:
            atlas_paths = glob(path.join(stats_path, '*statistics*'))
        else:
            atlas_paths = []
            for atlas_selected in atlas_selection:
                new_atlas_path = glob(path.join(stats_path, '*' + atlas_selected + '*'))
                atlas_paths = atlas_paths + new_atlas_path
                if len(new_atlas_path) == 0:
                    warnings.warn('The atlas wanted does not exist for ' + pipeline_name + ': ' + atlas_selected + ' in group ' + group_id, UserWarning)

        if len(atlas_paths) == 0:
                raise InitException(pipeline_name)

        atlas_paths.sort()
        for atlas_path in atlas_paths:
            row_summary_df = pd.DataFrame(index=np.arange(1), columns=columns)
            atlas_file = atlas_path.split(os.sep)[-1]
            atlas_id = '_'.join(atlas_file.split('_')[2:-1])  # We get rid of the subject and the session in the label

            init_df = pd.read_csv(atlas_path, sep='\t')
            n_regions = len(init_df['label_name'].tolist())
            label_list = [group_id + '_' + atlas_id + '_ROI-' + str(x) for x in range(n_regions)]
            for label in label_list:
                col_list.append(label)

            row_summary = [pipeline_name, group_id, atlas_id, n_regions, label_list[0], label_list[-1]]
            row_summary_df.iloc[0] = row_summary
            summary_df = pd.concat([summary_df, row_summary_df])

    # Initial empty DataFrame for the pipeline
    number_sessions = len(df)
    pipeline_df = pd.DataFrame(index=np.arange(number_sessions), columns=col_list)

    for i in range(number_sessions):
        # Reading the reference (subject + session)
        participant_id = df['participant_id'][i]
        session_id = df['session_id'][i]

        ses_path = caps_dir + os.sep + participant_id + os.sep + session_id
        mod_path = ses_path + os.sep + t1_spm_path
        ses_df = pd.DataFrame(index=np.arange(1), columns=col_list)

        if os.path.exists(mod_path):
            # Looking for groups
            for group in group_list:
                group_path = mod_path + os.sep + group
                if os.path.exists(group_path):
                    # Looking for atlases
                    atlas_list = summary_df[summary_df.group_id == group]['atlas_id'].values
                    for atlas in atlas_list:
                        atlas_path = group_path + os.sep + 'atlas_statistics' + os.sep
                        atlas_path = atlas_path + participant_id + '_' + session_id + '_' + atlas + '_statistics.tsv'
                        if os.path.exists(atlas_path):
                            n_regions = summary_df[(summary_df.group_id == group) &
                                                   (summary_df.atlas_id == atlas)]['regions_number'].values[0]
                            atlas_df = pd.read_csv(atlas_path, sep='\t')
                            label_list = [group + '_' + atlas + '_ROI-' + str(x) for x in
                                          range(n_regions)]
                            ses_df.loc[:, label_list] = atlas_df['mean_scalar'].as_matrix()

        pipeline_df.iloc[i] = ses_df.iloc[0]

    final_df = pd.concat([df, pipeline_df], axis=1)

    return final_df, summary_df


class InitException(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return repr('The following pipeline cannot be initialized: ' + self.name)


class DatasetError(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return repr('Bad format for the sessions: ' + self.name)
