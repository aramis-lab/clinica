# coding: utf8


def write(df, path, name):
    import os

    name = os.path.join(path, name) + ".tsv"
    df.to_csv(sep="\t", path_or_buf=name, index=False)


def complete_info_clinical(path_to_clinical, dfClinicalDict):
    """Determine where the clinical data that are not encoded in the Clinica standard have to be written under the BIDS format.

    Args:
       path_to_clinical: Path to the original data dictionnary 'DataDictionary_NIFD_2017.10.18.xlsx'
       dfClinicalDict: Pandas dataframe of 'clinical_info'
    """
    import pandas as pd

    def get_unclassified(dfClinicalDict):
        """Return the list of the remaining clinical data categories not encoded in the Clinica standard."""
        return list(dfClinicalDict[dfClinicalDict["BIDS_CLINICA"] == ""]["COLUMN_NAME"])

    def dicho_participants(path_to_clinical, possible_participants):
        others = []

        dfClinicalDict = pd.read_excel(path_to_clinical)
        list_patient = list(set(dfClinicalDict["LONI_ID"]))
        list_patient.sort()

        for pat in list_patient:
            df_pat = dfClinicalDict[dfClinicalDict["LONI_ID"] == pat]
            for elt in possible_participants:
                if len(set(df_pat[elt])) != 1:
                    others.append(elt)
                    possible_participants.remove(elt)

        return possible_participants, others

    possible_participants = get_unclassified(dfClinicalDict)
    participants, others = dicho_participants(path_to_clinical, possible_participants)

    return participants, others


def update_info_clinical(
    path_data_dict, path_clinicals, path_to_clinical, path_preprocessing
):
    """Build or update 'clinical_info.tsv' in the preprocessing folder.

    Args:
        path_data_dict (str): Path to the original data dictionnary 'DataDictionary_NIFD_2017.10.18.xlsx'
        path_clinicals (str): Path to the directory containing the Clinical data BIDS correspondence files
        path_to_clinical (str): [description]
        path_preprocessing (str): [description]
    """
    import os

    import numpy as np
    import pandas as pd

    dfParticipant = pd.read_csv(
        os.path.join(
            path_clinicals, "clinical_data_bids_correspondence_participant.tsv"
        ),
        sep="\t",
    )
    dfSessions = pd.read_csv(
        os.path.join(path_clinicals, "clinical_data_bids_correspondence_sessions.tsv"),
        sep="\t",
    )

    dfClinicalDict = pd.read_excel(path_data_dict)

    dfClinicalDict.insert(1, "LOCATION", "")
    dfClinicalDict.insert(1, "FILE", "")
    dfClinicalDict.insert(1, "BIDS_CLINICA", "")

    def add_row(df, name):
        serie = [np.nan] * len(df.columns.values)
        serie[0] = name
        serie = pd.Series(serie, index=list(df.columns.values))
        return df.append(serie, ignore_index=True)

    dfClinicalDict = add_row(dfClinicalDict, "Age")
    dfClinicalDict = add_row(dfClinicalDict, "Research Group")
    dfClinicalDict = add_row(dfClinicalDict, "Weight")

    def completeClinicalDict(dfCD, df, file):
        # Probably not the most elegant way of handling the issue but it works
        for index, row in dfCD.iterrows():
            for index2, row2 in df.iterrows():
                if row["COLUMN_NAME"] == row2["NIFD"]:
                    row["BIDS_CLINICA"] = row2["BIDS CLINICA"]
                    row["FILE"] = file
                    row["LOCATION"] = row2["NIFD location"]

        return dfCD

    # This part handles the clinical info not handled by the Clinica (this software) standard.
    dfClinicalDict = completeClinicalDict(dfClinicalDict, dfParticipant, "participants")
    dfClinicalDict = completeClinicalDict(dfClinicalDict, dfSessions, "sessions")

    participants, sessions = complete_info_clinical(path_to_clinical, dfClinicalDict)
    for part in participants:
        dfClinicalDict.loc[
            dfClinicalDict["COLUMN_NAME"] == part, "BIDS_CLINICA"
        ] = part.lower()
        dfClinicalDict.loc[
            dfClinicalDict["COLUMN_NAME"] == part, "FILE"
        ] = "participants"
        dfClinicalDict.loc[
            dfClinicalDict["COLUMN_NAME"] == part, "LOCATION"
        ] = "NIFD_Clinical_Data_2017_final_updated.xlsx"

    for ses in sessions:
        dfClinicalDict.loc[
            dfClinicalDict["COLUMN_NAME"] == ses, "BIDS_CLINICA"
        ] = ses.lower()
        dfClinicalDict.loc[dfClinicalDict["COLUMN_NAME"] == ses, "FILE"] = "sessions"
        dfClinicalDict.loc[
            dfClinicalDict["COLUMN_NAME"] == ses, "LOCATION"
        ] = "NIFD_Clinical_Data_2017_final_updated.xlsx"

    write(dfClinicalDict, path_preprocessing, "clinical_info")
