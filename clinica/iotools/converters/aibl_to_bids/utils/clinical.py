from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

__all__ = [
    "create_participants_tsv_file",
    "create_scans_tsv_file",
    "create_sessions_tsv_file",
]


def create_participants_tsv_file(
    input_path: Path,
    clinical_specifications_folder: Path,
    clinical_data_dir: Path,
    delete_non_bids_info: bool = True,
) -> None:
    """Create a participants TSV file for the AIBL dataset where information
    regarding the patients are reported.

    Parameters
    ----------
    input_path : Path
        The path to the input directory.

    clinical_specifications_folder : Path
        The path to the folder containing the clinical specification files.

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

    from clinica.iotools.bids_utils import StudyName, bids_id_factory

    fields_bids = ["participant_id"]
    fields_dataset = []
    prev_location = ""
    prev_sheet = ""
    index_to_drop = []

    specifications = _load_specifications(
        clinical_specifications_folder, "participant.tsv"
    )
    participant_fields_db = specifications[StudyName.AIBL.value]
    field_location = specifications[f"{StudyName.AIBL.value} location"]
    participant_fields_bids = specifications["BIDS CLINICA"]

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
    participant_df["participant_id"] = participant_df["alternative_id_1"].apply(
        lambda x: bids_id_factory(StudyName.AIBL).from_original_study_id(x)
    )
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


def _load_specifications(
    clinical_specifications_folder: Path, filename: str
) -> pd.DataFrame:
    specifications = clinical_specifications_folder / filename
    if not specifications.exists():
        raise FileNotFoundError(
            f"The specifications for {filename} were not found. "
            f"The should be located in {specifications}."
        )
    return pd.read_csv(specifications, sep="\t")


def _mapping_diagnosis(diagnosis: int) -> str:
    if diagnosis == 1:
        return "CN"
    elif diagnosis == 2:
        return "MCI"
    elif diagnosis == 3:
        return "AD"
    else:
        return "n/a"


def _extract_metadata_df(
    input_df: pd.DataFrame, source_id: int, bids_metadata: str, source_metadata: str
) -> pd.DataFrame:
    from clinica.iotools.converter_utils import viscode_to_session

    extract = input_df.loc[(input_df["RID"] == source_id), ["VISCODE", source_metadata]]
    extract.rename(columns={source_metadata: bids_metadata}, inplace=True)
    extract = extract.assign(
        session_id=extract.VISCODE.apply(lambda x: viscode_to_session(x))
    )
    extract.drop(labels="VISCODE", inplace=True, axis=1)
    extract.set_index("session_id", inplace=True, drop=True)

    return extract


def _compute_age_at_exam(
    birth_date: Optional[str], exam_date: Optional[str]
) -> Optional[float]:
    from datetime import datetime

    if birth_date and exam_date:
        date_of_birth = datetime.strptime(birth_date, "/%Y")
        exam_date = datetime.strptime(exam_date, "%m/%d/%Y")
        return exam_date.year - date_of_birth.year
    return None


def create_sessions_tsv_file(
    bids_dir: Path,
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
) -> None:
    """Extract the information regarding a subject sessions and save them in a tsv file.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    clinical_data_dir : Path
        The path to the directory to the clinical data files.

    clinical_specifications_folder : Path
        The path to the folder containing the clinical specification files.
    """
    # todo :rename test
    import glob

    from clinica.iotools.bids_utils import (
        StudyName,
        bids_id_factory,
        get_bids_sess_list,
        get_bids_subjs_list,
    )

    study = StudyName.AIBL.value
    specifications = _load_specifications(
        clinical_specifications_folder, "sessions.tsv"
    )[["BIDS CLINICA", f"{study} location", study]].dropna()

    for bids_id in get_bids_subjs_list(bids_dir):
        rid = int(bids_id_factory(study)(bids_id).to_original_study_id())
        test = pd.DataFrame(
            {"session_id": get_bids_sess_list(bids_dir / bids_id)}
        ).set_index("session_id", drop=False)

        for _, row in specifications.iterrows():
            df = pd.read_csv(
                glob.glob(str(clinical_data_dir / row[f"{study} location"]))[0],
                dtype={"text": str},
            )
            extract = _extract_metadata_df(
                input_df=df,
                source_id=rid,
                bids_metadata=row["BIDS CLINICA"],
                source_metadata=row[study],
            )
            test = pd.concat([test, extract], axis=1)

        test.sort_index(inplace=True)
        # -4 are considered missing values in AIBL
        test.replace([-4, "-4", np.nan], None, inplace=True)
        test["diagnosis"] = test.diagnosis.apply(lambda x: _mapping_diagnosis(x))
        test["examination_date"] = test.apply(
            lambda x: _complete_examination_dates(
                rid, x.session_id, x.examination_date, clinical_data_dir
            ),
            axis=1,
        )

        # in general age metadata is present only for baseline session
        test["age"] = test["age"].ffill()
        test["age"] = test.apply(
            lambda x: _compute_age_at_exam(x.age, x.examination_date), axis=1
        )

        # in case there is a session in clinical data that was not actually converted
        test.dropna(subset=["session_id"], inplace=True)
        test.fillna("n/a", inplace=True)

        bids_id = bids_id_factory(StudyName.AIBL).from_original_study_id(str(rid))
        test.to_csv(
            bids_dir / bids_id / f"{bids_id}_sessions.tsv",
            sep="\t",
            index=False,
            encoding="utf8",
        )


def _complete_examination_dates(
    rid: int,
    session_id: Optional[str],
    examination_date: Optional[str],
    clinical_data_dir: Path,
) -> Optional[str]:
    if examination_date:
        return examination_date
    if session_id:
        return _find_exam_date_in_other_csv_files(rid, session_id, clinical_data_dir)
    return None


def _find_exam_date_in_other_csv_files(
    rid: int, session_id: str, clinical_data_dir: Path
) -> Optional[str]:
    """Try to find an alternative exam date by searching in other CSV files."""
    from clinica.iotools.converter_utils import viscode_to_session

    for csv in _get_csv_paths(clinical_data_dir):
        csv_data = pd.read_csv(csv, low_memory=False)
        csv_data["SESSION"] = csv_data.VISCODE.apply(lambda x: viscode_to_session(x))
        exam_date = csv_data[(csv_data.RID == rid) & (csv_data.SESSION == session_id)]
        if not exam_date.empty and exam_date.iloc[0].EXAMDATE != "-4":
            return exam_date.iloc[0].EXAMDATE
    return None


def _get_csv_paths(clinical_data_dir: Path) -> Tuple[str]:
    """Return a list of paths to CSV files in which an alternative exam date could be found."""
    import glob

    pattern_list = (
        "aibl_mri3meta_*.csv",
        "aibl_mrimeta_*.csv",
        "aibl_cdr_*.csv",
        "aibl_flutemeta_*.csv",
        "aibl_mmse_*.csv",
        "aibl_pibmeta_*.csv",
    )

    paths_list = ()

    for pattern in pattern_list:
        try:
            path = glob.glob(str(clinical_data_dir / pattern))[0]
            paths_list += (path,)
        except IndexError:
            pass
    return paths_list


def create_scans_tsv_file(
    input_path: Path,
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
) -> None:
    """Create scans.tsv files for AIBL.

    Parameters
    ----------
    input_path : Path
        The path to the input folder.

    clinical_data_dir : Path
        The path to the directory to the clinical data files.

    clinical_specifications_folder : Path
        The path to the folder containing the clinical specification files.
    """
    import glob
    from os import path

    import clinica.iotools.bids_utils as bids

    specifications = _load_specifications(clinical_specifications_folder, "scans.tsv")
    scans_fields = specifications[bids.StudyName.AIBL.value]
    field_location = specifications[f"{bids.StudyName.AIBL.value} location"]
    scans_fields_bids = specifications["BIDS CLINICA"]
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
        bids.StudyName.AIBL,
        clinical_specifications_folder,
        bids_ids,
        "RID",
        "VISCODE",
        ses_dict,
    )
    bids.write_scans_tsv(input_path, bids_ids, scans_dict)
