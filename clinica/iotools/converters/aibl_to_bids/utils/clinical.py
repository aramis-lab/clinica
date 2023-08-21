def create_participants_tsv_file(
    input_path: str,
    clinical_spec_path: str,
    clinical_data_dir: str,
    delete_non_bids_info: bool = True,
):
    """Create a participants file for the AIBL dataset where information
    regarding the patients are reported.

    Parameters
    ----------
    input_path : str
        Path to the input directory.

    clinical_spec_path : str
        Path to the clinical file.

    clinical_data_dir : str
        Path to the directory to the clinical data files.

    delete_non_bids_info : bool
        If True delete all the rows of the subjects that are not
        available in the BIDS dataset.
    """
    import glob
    import os
    from os import path

    import numpy as np
    import pandas as pd

    fields_bids = ["participant_id"]
    fields_dataset = []
    prev_location = ""
    prev_sheet = ""
    index_to_drop = []

    location_name = "AIBL location"
    clinical_spec_path = clinical_spec_path + "_participant.tsv"

    if not os.path.exists(clinical_spec_path):
        raise FileNotFoundError(clinical_spec_path + " not found in clinical data.")
    participants_specs = pd.read_csv(clinical_spec_path, sep="\t")
    participant_fields_db = participants_specs["AIBL"]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs["BIDS CLINICA"]

    # Extract the list of the available fields for the dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])
            fields_dataset.append(participant_fields_db[i])

    # Init the dataframe that will be saved in the file participant.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

    for i in range(0, len(participant_fields_db)):
        # If a field not empty is found
        if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            # If a sheet is available
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

                if file_ext == ".xlsx":
                    file_to_read = pd.read_excel(
                        glob.glob(file_to_read_path)[0], sheet_name=sheet
                    )
                elif file_ext == ".csv":
                    file_to_read = pd.read_csv(glob.glob(file_to_read_path)[0])
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                if participant_fields_bids[i] == "alternative_id_1" and (
                    file_to_read[participant_fields_db[i]].dtype == np.float64
                    or file_to_read[participant_fields_db[i]].dtype == np.int64
                ):
                    if not pd.isnull(file_to_read.at[j, participant_fields_db[i]]):
                        # value_to_append = str(file_to_read.get_value(j, participant_fields_db[i])).rstrip('.0')
                        value_to_append = str(
                            file_to_read.at[j, participant_fields_db[i]]
                        )

                    else:
                        value_to_append = "n/a"
                else:
                    value_to_append = file_to_read.at[j, participant_fields_db[i]]
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    # Compute BIDS-compatible participant ID.
    participant_df["participant_id"] = "sub-AIBL" + participant_df["alternative_id_1"]

    # Keep year-of-birth only.
    participant_df["date_of_birth"] = participant_df["date_of_birth"].str.extract(
        r"/(\d{4}).*"
    )

    # Normalize sex value.
    participant_df["sex"] = participant_df["sex"].map({1: "M", 2: "F"}).fillna("n/a")

    # Normalize known NA values.
    participant_df.replace(-4, "n/a", inplace=True)

    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.drop(index_to_drop)

    participant_df.to_csv(
        os.path.join(input_path, "participants.tsv"),
        sep="\t",
        index=False,
        encoding="utf8",
    )


def create_sessions_tsv_file(
    input_path: str, clinical_data_dir: str, clinical_spec_path: str
) -> None:
    """Extract the information regarding the sessions and save them in a tsv file.

    Parameters
    ----------
    input_path : str
        Path to the input folder.

    clinical_data_dir : str
        Path to the directory to the clinical data files.

    clinical_spec_path : str
        Path to the clinical file.
    """
    import glob
    from os import path

    import pandas as pd

    # Load data
    location = "AIBL location"
    sessions = pd.read_csv(clinical_spec_path + "_sessions.tsv", sep="\t")
    sessions_fields = sessions["AIBL"]
    field_location = sessions[location]
    sessions_fields_bids = sessions["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    files_to_read = []
    sessions_fields_to_read = []
    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i]
            file_to_read_path = path.join(clinical_data_dir, tmp)
            files_to_read.append(glob.glob(file_to_read_path)[0])
            sessions_fields_to_read.append(sessions_fields[i])

    rid = pd.read_csv(files_to_read[0], dtype={"text": str}, low_memory=False).RID
    rid = list(set(rid))
    for r in rid:
        for i in files_to_read:
            file_to_read = pd.read_csv(i, dtype={"text": str})
            if len(file_to_read.columns) == 1:
                file_to_read = pd.read_csv(i, sep=";", low_memory=False)

            # information are written following the BIDS specifications
            viscode = file_to_read.loc[(file_to_read["RID"] == r), "VISCODE"]
            for j in sessions_fields_to_read:
                if j in list(file_to_read.columns.values) and j == "MMSCORE":
                    MMSCORE = file_to_read.loc[(file_to_read["RID"] == r), j]
                    MMSCORE[MMSCORE == -4] = "n/a"
                elif j in list(file_to_read.columns.values) and j == "CDGLOBAL":
                    CDGLOBAL = file_to_read.loc[(file_to_read["RID"] == r), j]
                    CDGLOBAL[CDGLOBAL == -4] = "n/a"
                elif j in list(file_to_read.columns.values) and j == "DXCURREN":
                    DXCURREN = file_to_read.loc[(file_to_read["RID"] == r), j]
                    DXCURREN[DXCURREN == -4] = "n/a"
                    DXCURREN[DXCURREN == 1] = "CN"
                    DXCURREN[DXCURREN == 2] = "MCI"
                    DXCURREN[DXCURREN == 3] = "AD"
                elif j in list(file_to_read.columns.values) and j == "EXAMDATE":
                    EXAMDATE = file_to_read.loc[(file_to_read["RID"] == r), j]
                elif j in list(file_to_read.columns.values) and j == "PTDOB":
                    PTDOB = file_to_read.loc[(file_to_read["RID"] == r), j]

        examdates = _get_examdates(
            r, EXAMDATE.to_list(), viscode.to_list(), clinical_data_dir
        )
        age = _get_ages(PTDOB.values[0], examdates)

        viscode[viscode == "bl"] = "M000"
        viscode = viscode.str.upper()

        sessions = pd.DataFrame(
            {
                "months": viscode.str[1:],
                "age": age,
                "MMS": MMSCORE,
                "cdr_global": CDGLOBAL,
                "diagnosis": DXCURREN,
                "examination_date": examdates,
            }
        )
        sessions = sessions.assign(
            session_id=lambda df: df.months.apply(lambda x: f"ses-M{int(x):03d}")
        )
        cols = sessions.columns.tolist()
        sessions = sessions[cols[-1:] + cols[:-1]]

        bids_paths = path.join(input_path, "sub-AIBL" + str(r))
        if path.exists(bids_paths):
            sessions.to_csv(
                path.join(
                    input_path,
                    "sub-AIBL" + str(r),
                    "sub-AIBL" + str(r) + "_sessions.tsv",
                ),
                sep="\t",
                index=False,
                encoding="utf8",
            )


def _get_examdates(rid, examdates, viscodes, clinical_data_dir):
    import glob
    from datetime import datetime
    from os import path

    import pandas as pd
    from dateutil.relativedelta import relativedelta

    res_examdates = []
    csv_list = [
        glob.glob(path.join(clinical_data_dir, "aibl_mri3meta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_mrimeta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_cdr_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_flutemeta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_mmse_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_pibmeta_*.csv"))[0],
    ]

    for e in range(len(examdates)):
        exam = examdates[e]

        if exam != "-4":
            res_examdates.append(exam)
            continue

        # If EXAMDATE does not exist (-4) we try to obtain it from another .csv file
        for csv_file in csv_list:
            if "aibl_flutemeta" in csv_file:
                csv_data = pd.read_csv(
                    csv_file, low_memory=False, usecols=list(range(0, 36))
                )
            else:
                csv_data = pd.read_csv(csv_file, low_memory=False)
            exam_date = csv_data[
                (csv_data.RID == rid) & (csv_data.VISCODE == viscodes[e])
            ]
            if not exam_date.empty and exam_date.iloc[0].EXAMDATE != "-4":
                exam = exam_date.iloc[0].EXAMDATE
                break

        # If EXAMDATE still does not exist (-4) we add the session months to baseline date
        if exam == "-4":
            bl_index = viscodes.index("bl")
            if bl_index > -1:
                bl_date = examdates[bl_index]
                bl_examdate = datetime.strptime(bl_date, "%m/%d/%Y")
                if viscodes[e] != "bl":
                    months = int(viscodes[e][1:])
                    examdate = bl_examdate + relativedelta(months=+months)
                    exam = examdate.strftime("%m/%d/%Y")

        if exam == "-4":
            print(f"No EXAMDATE for subject %{rid}, at session {viscodes[e]}")

        res_examdates.append(exam)

    return res_examdates


def _get_ages(pt_dob, examdates):
    """Calculate age as time passed by since DOB to EXAMDATE.

    :param pt_dob: string - Date of birth of patient ("/%Y" format)
    :param examdates: list - Exam dates ("%m/%d/%Y" format)
    :return: list - Age at each exam date
    """
    from datetime import datetime

    age = []
    dob = datetime.strptime(pt_dob, "/%Y")

    for exam in examdates:
        examdate = datetime.strptime(exam, "%m/%d/%Y")
        delta = examdate - dob
        age.append(round(delta.days / 365.25, 1))

    return age


def create_scans_tsv_file(
    input_path: str, clinical_data_dir: str, clinical_spec_path: str
) -> None:
    """Create scans.tsv files for AIBL.

    Parameters
    ----------
    input_path : str
        Path to the input folder.

    clinical_data_dir : str
        Path to the directory to the clinical data files.

    clinical_spec_path : str
        Path to the clinical file.
    """
    import glob
    from os import path

    import pandas as pd

    import clinica.iotools.bids_utils as bids

    # Load data
    location = "AIBL location"
    scans = pd.read_csv(clinical_spec_path + "_scans.tsv", sep="\t")
    scans_fields = scans["AIBL"]
    field_location = scans[location]
    scans_fields_bids = scans["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []

    # Keep only fields for which there are AIBL fields
    for i in range(0, len(scans_fields)):
        if not pd.isnull(scans_fields[i]):
            fields_bids.append(scans_fields_bids[i])
            fields_dataset.append(scans_fields[i])

    files_to_read = []
    sessions_fields_to_read = []
    for i in range(0, len(scans_fields)):
        # If the i-th field is available
        if not pd.isnull(scans_fields[i]):
            # Load the file
            tmp = field_location[i]
            file_to_read_path = path.join(clinical_data_dir, tmp)
            files_to_read.append(glob.glob(file_to_read_path)[0])
            sessions_fields_to_read.append(scans_fields[i])

    bids_ids = [
        path.basename(sub_path)
        for sub_path in glob.glob(path.join(input_path, "sub-AIBL*"))
    ]
    ses_dict = {
        bids_id: {"M000": "bl", "M018": "m18", "M036": "m36", "M054": "m54"}
        for bids_id in bids_ids
    }
    scans_dict = bids.create_scans_dict(
        clinical_data_dir,
        "AIBL",
        clinical_spec_path,
        bids_ids,
        "RID",
        "VISCODE",
        ses_dict,
    )
    bids.write_scans_tsv(input_path, bids_ids, scans_dict)
