# coding: utf8

import numpy as np
import pandas as pd
import nibabel as nib
import os

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "simona.bottani@icm-institute.com"
__status__ = "Development"


def load_data(images, caps_directory, subjects, sessions, dataset):
    """

    Args:
        images:
        caps_directory:
        subjects_visits_tsv:
        dataset:

    Returns:
        np 2D array

    """


    df = pd.io.parsers.read_csv(os.path.join(caps_directory), sep='\t')

    all_vector = np.array([])

    if dataset == 'OASIS':
        df = df[df.age_bl > 61]
    #subjects_visits = pd.io.parsers.read_csv(os.path.join(subjects_visits_tsv), sep='\t')
    participant_id = subjects
    session_id = sessions

    for i in xrange(len(participant_id)):
        df_sub = df[df.participant_id == participant_id[i]]

        df_analysis = df_sub[[col for col in df_sub.columns if images in col]]

        all_vector = np.append(all_vector, df_analysis.values)


    data = np.zeros((participant_id.shape[0], df_analysis.shape[1]))
    data_temp = np.split(all_vector, participant_id.shape[0])

    for i in xrange(len(participant_id)):
        for j in xrange(df_analysis.shape[1]):
            data[i][j] = data_temp[i][j]


    return data
