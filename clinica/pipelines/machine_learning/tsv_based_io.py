import os

import numpy as np
import pandas as pd


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

    df = pd.read_csv(os.path.join(caps_directory), sep="\t")

    all_vector = np.array([])

    participant_id = subjects

    session_id = sessions

    for i in range(len(participant_id)):
        df_sub = df[df.participant_id == participant_id[i]]

        df_analysis = df_sub[[col for col in df_sub.columns if images in col]]

        all_vector = np.append(all_vector, df_analysis.values)
    data = np.zeros((len(participant_id), df_analysis.shape[1]))
    data_temp = np.split(all_vector, len(participant_id))

    for i in range(len(participant_id)):
        for j in range(df_analysis.shape[1]):
            data[i][j] = data_temp[i][j]

    return data
