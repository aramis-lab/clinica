__author__ = "Simona Bottani and Sabrina Fontanella"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


# -- Methods for the clinical data --

def create_participants_df_AIBL(input_path,clinical_spec_path, clinical_data_dir, delete_non_bids_info=True):
    """
        This methods create a participants file for the AIBL dataset where information regarding the patients are reported
        
        :param input_path: path to the input directory
        :param clinical_spec_path: path to the clinical file
        :param clinical_data_dir: directory to the clinical data files
        :param delete_non_bids_info: if True delete all the rows of the subjects that are not available in the BIDS dataset
        :return: a pandas dataframe that contains the participants data and it is saved in a tsv file
    """
    import pandas as pd
    import os
    from os import path
    import re
    import numpy as np

    fields_bids = ['participant_id']
    fields_dataset = []
    prev_location = ''
    index_to_drop=[]

    location_name ='AIBL location'

    participants_specs = pd.read_excel(clinical_spec_path , sheetname='participant.tsv')
    participant_fields_db = participants_specs['AIBL']
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs['BIDS CLINICA']

    # Extract the list of the available fields for the dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])
            fields_dataset.append(participant_fields_db[i])

    # Init the dataframe that will be saved in the file participant.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

    csv_files=[]
    for i in range(0, len(participant_fields_db)):
    # If a field not empty is found
        if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            tmp = field_location[i].split('/')
            location = tmp[0]
            # If a sheet is available
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ''
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

                if file_ext == '.xlsx':
                    file_to_read = pd.read_excel(file_to_read_path, sheetname=sheet)
                elif file_ext == '.csv':
                    file_to_read = pd.read_csv(file_to_read_path)
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                if participant_fields_bids[i] == 'alternative_id_1' and\
                        (file_to_read[participant_fields_db[i]].dtype == np.float64 or file_to_read[participant_fields_db[i]].dtype == np.int64) :
                    if not pd.isnull(file_to_read.get_value(j, participant_fields_db[i])):
                        #value_to_append = str(file_to_read.get_value(j, participant_fields_db[i])).rstrip('.0')
                        value_to_append = str(file_to_read.get_value(j, participant_fields_db[i]))

                    else:
                        value_to_append = np.NaN
                else:
                    value_to_append = file_to_read.get_value(j, participant_fields_db[i])
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)


    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        value=participant_df['alternative_id_1'][i]
        participant_df['participant_id'][i] = 'sub-AIBL' + str(value)
        participant_df['date_of_birth'][i]=re.search('/([0-9].*)', str(participant_df['date_of_birth'][i])).group(1)
        if participant_df['sex'][i]==1:
            participant_df['sex'][i]='M'
        else:
            participant_df['sex'][i]='F'

    participant_df.replace('-4', np.nan)

    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info == True:
        participant_df = participant_df.drop(index_to_drop)

    participant_df.to_csv(os.path.join(input_path, 'participants.tsv'), sep='\t', index=False, encoding='utf8')

    return participant_df


def create_sessions_dict_AIBL(input_path, clinical_data_dir,clinical_spec_path):
    """
        Extract the information regarding the sessions and store them in a dictionary (session M0 only)

        :param input_path: path to the input folder
        :param clinical_spec_path: path to the clinical file
        :param clinical_data_dir: directory to the clinical data files
        :return: A dataframe saved in a tsv file which contains information for each session
    """
    import pandas as pd
    from os import path
    import numpy as np

    # Load data
    location = 'AIBL location'
    sessions = pd.read_excel(clinical_spec_path, sheetname='sessions.tsv')
    sessions_fields = sessions['AIBL']
    field_location = sessions[location]
    sessions_fields_bids = sessions['BIDS CLINICA']
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    sessions_df = pd.DataFrame(columns=fields_bids)

    files_to_read=[]
    sessions_fields_to_read=[]
    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp=field_location[i]
            location = tmp[0]
            file_to_read_path=path.join(clinical_data_dir,tmp)
            files_to_read.append(file_to_read_path)
            sessions_fields_to_read.append(sessions_fields[i])

    rid = pd.read_csv(files_to_read[0], dtype={'text': unicode}).RID
    rid = list(set(rid))
    for r in rid:
        dict = []
        for i in files_to_read:
            file_to_read = pd.read_csv(i, dtype={'text': unicode})
            #information are written following the BIDS specifications
            viscode = file_to_read.loc[(file_to_read["RID"] == r), "VISCODE"]
            viscode[viscode == 'bl'] = 'M00'
            viscode[viscode == 'm18'] = 'M18'
            viscode[viscode == 'm36'] = 'M36'
            viscode[viscode == 'm54'] = 'M54'
            for i in sessions_fields_to_read:
                if i in list(file_to_read.columns.values) and i=='MMSCORE':
                    MMSCORE = file_to_read.loc[(file_to_read["RID"] == r), i]
                    MMSCORE[MMSCORE == -4] = np.nan
                elif i in list(file_to_read.columns.values) and i=='DXCURREN':
                    DXCURREN = file_to_read.loc[(file_to_read["RID"] == r), i]
                    DXCURREN[DXCURREN == -4] = np.nan
                    DXCURREN[DXCURREN == 1] = 'CN'
                    DXCURREN[DXCURREN == 2] = 'MCI'
                    DXCURREN[DXCURREN == 3] = 'AD'
                elif i in list(file_to_read.columns.values) and i == 'EXAMDATE':
                    EXAMDATE = file_to_read.loc[(file_to_read["RID"] == r), i]
        dict = pd.DataFrame({'session_id': viscode,
                                         'MMS': MMSCORE,
                                         'diagnosis': DXCURREN,
                                        'examination_date':EXAMDATE
                                         })

        cols = dict.columns.tolist()
        dict = dict[cols[-1:] + cols[:-1]]

        bids_paths=path.join(input_path,'sub-AIBL'+str(r))
        if path.exists(bids_paths):
            dict.to_csv(path.join(input_path,'sub-AIBL'+str(r), 'sub-AIBL'+str(r)+ '_sessions.tsv'),sep='\t', index=False, encoding='utf8')

