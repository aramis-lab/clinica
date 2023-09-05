from pathlib import Path
from typing import List, Optional


def create_participants_tsv_file(
    input_path: Path,
    clinical_spec_path: str,
    clinical_data_dir: Path,
    delete_non_bids_info: bool = True,
) -> None:
    """Create a participants TSV file for the AIBL dataset where information
    regarding the patients are reported.

    Parameters
    ----------
    input_path : Path
        The path to the input directory.

    clinical_spec_path : str
        The path to the clinical file.

    clinical_data_dir : Path
        The path to the directory to the clinical data files.

    delete_non_bids_info : bool, optional
        If True delete all the rows of the subjects that are not
        available in the BIDS dataset.
        Default=True.
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

    clinical_spec_path = Path(clinical_spec_path + "_participant.tsv")

    if not clinical_spec_path.exists():
        raise FileNotFoundError(f"{clinical_spec_path} not found in clinical data.")

    participants_specs = pd.read_csv(clinical_spec_path, sep="\t")
    participant_fields_db = participants_specs["AIBL"]
    field_location = participants_specs["AIBL location"]
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
            sheet = tmp[1] if len(tmp) > 1 else ""
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
        input_path / "participants.tsv",
        sep="\t",
        index=False,
        encoding="utf8",
    )


def create_sessions_tsv_file(
    input_path: Path, clinical_data_dir: Path, clinical_spec_path: str
) -> None:
    """Extract the information regarding the sessions and save them in a tsv file.

    Parameters
    ----------
    input_path : Path
        The path to the input folder.

    clinical_data_dir : Path
        The path to the directory to the clinical data files.

    clinical_spec_path : str
        Path to the clinical file.
    """
    import glob

    import pandas as pd

    sessions = pd.read_csv(clinical_spec_path + "_sessions.tsv", sep="\t")
    sessions_fields = sessions["AIBL"]
    field_location = sessions["AIBL location"]
    sessions_fields_bids = sessions["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    files_to_read: List[str] = []
    sessions_fields_to_read: List[str] = []
    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            file_to_read_path = clinical_data_dir / field_location[i]
            files_to_read.append(glob.glob(str(file_to_read_path))[0])
            sessions_fields_to_read.append(sessions_fields[i])

    rids = pd.read_csv(files_to_read[0], dtype={"text": str}, low_memory=False).RID
    unique_rids = list(set(rids))
    for rid in unique_rids:
        for file in files_to_read:
            df = pd.read_csv(file, dtype={"text": str})
            if len(df.columns) == 1:
                df = pd.read_csv(file, sep=";", low_memory=False)

            visit_code = df.loc[(df["RID"] == rid), "VISCODE"]

            for field in sessions_fields_to_read:

                if field in list(df.columns.values) and field == "MMSCORE":
                    mm_score = df.loc[(df["RID"] == rid), field]
                    mm_score[mm_score == -4] = "n/a"

                elif field in list(df.columns.values) and field == "CDGLOBAL":
                    cd_global = df.loc[(df["RID"] == rid), field]
                    cd_global[cd_global == -4] = "n/a"

                elif field in list(df.columns.values) and field == "DXCURREN":
                    dx_curren = df.loc[(df["RID"] == rid), field]
                    dx_curren[dx_curren == -4] = "n/a"
                    dx_curren[dx_curren == 1] = "CN"
                    dx_curren[dx_curren == 2] = "MCI"
                    dx_curren[dx_curren == 3] = "AD"

                elif field in list(df.columns.values) and field == "EXAMDATE":
                    exam_date = df.loc[(df["RID"] == rid), field]

                elif field in list(df.columns.values) and field == "PTDOB":
                    patient_date_of_birth = df.loc[(df["RID"] == rid), field]

        exam_dates = _clean_exam_dates(
            rid, exam_date.to_list(), visit_code.to_list(), clinical_data_dir
        )
        age = _compute_ages_at_each_exam(patient_date_of_birth.values[0], exam_dates)

        visit_code[visit_code == "bl"] = "M000"
        visit_code = visit_code.str.upper()

        sessions = pd.DataFrame(
            {
                "months": visit_code.str[1:],
                "age": age,
                "MMS": mm_score,
                "cdr_global": cd_global,
                "diagnosis": dx_curren,
                "examination_date": exam_dates,
            }
        )
        sessions = sessions.assign(
            session_id=lambda df: df.months.apply(lambda x: f"ses-M{int(x):03d}")
        )
        cols = sessions.columns.tolist()
        sessions = sessions[cols[-1:] + cols[:-1]]

        bids_paths = input_path / f"sub-AIBL{rid}"
        if bids_paths.exists():
            sessions.to_csv(
                input_path / f"sub-AIBL{rid}" / f"sub-AIBL{rid}_sessions.tsv",
                sep="\t",
                index=False,
                encoding="utf8",
            )


def _clean_exam_dates(
    rid: str, exam_dates: List[str], visit_codes: List[str], clinical_data_dir: Path
) -> List[str]:
    """Clean the exam dates when necessary by trying to compute them from other sources."""
    from clinica.utils.stream import cprint

    exam_dates_cleaned: List[str] = []
    for visit_code, exam_date in zip(visit_codes, exam_dates):
        if exam_date == "-4":
            exam_date = _find_exam_date_in_other_csv_files(
                rid, visit_code, clinical_data_dir
            ) or _compute_exam_date_from_baseline(visit_code, exam_dates, visit_codes)
            if not exam_date:
                cprint(f"No EXAMDATE for subject %{rid}, at session {visit_code}")
                exam_date = "-4"
        exam_dates_cleaned.append(exam_date)

    return exam_dates_cleaned


def _find_exam_date_in_other_csv_files(
    rid: str, visit_code: str, clinical_data_dir: Path
) -> Optional[str]:
    """Try to find an alternative exam date by searching in other CSV files."""
    import pandas as pd

    for csv_file in _get_cvs_files(clinical_data_dir):
        if "aibl_flutemeta" in csv_file:
            csv_data = pd.read_csv(
                csv_file, low_memory=False, usecols=list(range(0, 36))
            )
        else:
            csv_data = pd.read_csv(csv_file, low_memory=False)
        exam_date = csv_data[(csv_data.RID == rid) & (csv_data.VISCODE == visit_code)]
        if not exam_date.empty and exam_date.iloc[0].EXAMDATE != "-4":
            return exam_date.iloc[0].EXAMDATE
    return None


def _get_cvs_files(clinical_data_dir: Path) -> List[str]:
    """Return a list of paths to CSV files in which an alternative exam date could be found."""
    import glob

    return [
        glob.glob(str(clinical_data_dir / pattern))[0]
        for pattern in (
            "aibl_mri3meta_*.csv",
            "aibl_mrimeta_*.csv",
            "aibl_cdr_*.csv",
            "aibl_flutemeta_*.csv",
            "aibl_mmse_*.csv",
            "aibl_pibmeta_*.csv",
        )
    ]


def _compute_exam_date_from_baseline(
    visit_code: str, exam_dates: List[str], visit_codes: List[str]
) -> Optional[str]:
    """Try to find an alternative exam date by computing the number of months from the visit code."""
    from datetime import datetime

    from dateutil.relativedelta import relativedelta

    baseline_index = visit_codes.index("bl")
    if baseline_index > -1:
        baseline_date = datetime.strptime(exam_dates[baseline_index], "%m/%d/%Y")
        if visit_code != "bl":
            try:
                months = int(visit_code[1:])
            except TypeError:
                raise ValueError(
                    f"Unexpected visit code {visit_code}. Should be in format MXXX."
                    "Ex: M000, M006, M048..."
                )
            exam_date = baseline_date + relativedelta(months=+months)
            return exam_date.strftime("%m/%d/%Y")
    return None


def _compute_ages_at_each_exam(
    patient_date_of_birth: str, exam_dates: List[str]
) -> List[float]:
    """Compute the ages of the patient at each exam date.

    Parameters
    ----------
    patient_date_of_birth : str
        Date of birth of patient ("/%Y" format)

    exam_dates : list of str
        List of exam dates ("%m/%d/%Y" format)

    Return
    ------
    list of float
        The ages of the patient at each exam date.
    """
    from datetime import datetime

    ages: List[float] = []
    date_of_birth = datetime.strptime(patient_date_of_birth, "/%Y")

    for exam_date in exam_dates:
        exam_date = datetime.strptime(exam_date, "%m/%d/%Y")
        delta = exam_date - date_of_birth
        ages.append(round(delta.days / 365.25, 1))

    return ages


def create_scans_tsv_file(
    input_path: Path, clinical_data_dir: Path, clinical_spec_path: str
) -> None:
    """Create scans.tsv files for AIBL.

    Parameters
    ----------
    input_path : Path
        The path to the input folder.

    clinical_data_dir : Path
        The path to the directory to the clinical data files.

    clinical_spec_path : str
        Path to the clinical file.
    """
    import glob
    from os import path

    import pandas as pd

    import clinica.iotools.bids_utils as bids

    scans = pd.read_csv(clinical_spec_path + "_scans.tsv", sep="\t")
    scans_fields = scans["AIBL"]
    field_location = scans["AIBL location"]
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
            file_to_read_path = clinical_data_dir / field_location[i]
            files_to_read.append(glob.glob(str(file_to_read_path))[0])
            sessions_fields_to_read.append(scans_fields[i])

    bids_ids = [
        path.basename(sub_path) for sub_path in glob.glob(str(input_path / "sub-AIBL*"))
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
