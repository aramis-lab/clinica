from pathlib import Path
from typing import Iterator, List, Optional

import numpy as np
import pandas as pd

from clinica.iotools.bids_utils import StudyName, bids_id_factory

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


def _map_diagnosis(diagnosis: int) -> str:
    if diagnosis == 1:
        return "CN"
    elif diagnosis == 2:
        return "MCI"
    elif diagnosis == 3:
        return "AD"
    else:
        return "n/a"


def _format_metadata_for_rid(
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
) -> Optional[int]:
    from datetime import datetime

    if birth_date and exam_date:
        date_of_birth = datetime.strptime(birth_date, "/%Y")
        exam_date = datetime.strptime(exam_date, "%m/%d/%Y")
        return exam_date.year - date_of_birth.year
    return None


def _set_age_from_birth(df: pd.DataFrame) -> pd.DataFrame:
    if "date_of_birth" not in df.columns or "examination_date" not in df.columns:
        raise ValueError(
            "Columns date_of_birth or/and examination_date were not found in the sessions metadata dataframe."
            "Please check your study metadata."
        )
    if len(df["date_of_birth"].dropna().unique()) <= 1:
        df["date_of_birth"] = df["date_of_birth"].ffill()
        df["age"] = df.apply(
            lambda x: _compute_age_at_exam(x.date_of_birth, x.examination_date), axis=1
        )
    else:
        df["age"] = None
    return df.drop(labels="date_of_birth", axis=1)


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
        sessions = pd.DataFrame(
            {"session_id": get_bids_sess_list(bids_dir / bids_id)}
        ).set_index("session_id", drop=False)

        for _, row in specifications.iterrows():
            try:
                df = pd.read_csv(
                    next(clinical_data_dir.glob(row[f"{study} location"])),
                    dtype={"text": str},
                )
            except StopIteration:
                raise FileNotFoundError(
                    f"Clinical data file corresponding to pattern {row[f'{study} location']} was not found in folder "
                    f"{clinical_data_dir}"
                )

            data = _format_metadata_for_rid(
                input_df=df,
                source_id=rid,
                bids_metadata=row["BIDS CLINICA"],
                source_metadata=row[study],
            )
            sessions = pd.concat([sessions, data], axis=1)

        sessions.sort_index(inplace=True)

        # -4 are considered missing values in AIBL
        sessions.replace([-4, "-4", np.nan], None, inplace=True)
        sessions["diagnosis"] = sessions.diagnosis.apply(lambda x: _map_diagnosis(x))
        sessions["examination_date"] = sessions.apply(
            lambda x: _complete_examination_dates(
                rid, clinical_data_dir, x.session_id, x.examination_date
            ),
            axis=1,
        )
        sessions = _set_age_from_birth(sessions)

        # in case there is a session in clinical data that was not actually converted
        sessions.dropna(subset=["session_id"], inplace=True)
        sessions.fillna("n/a", inplace=True)

        bids_id = bids_id_factory(StudyName.AIBL).from_original_study_id(str(rid))
        sessions.to_csv(
            bids_dir / bids_id / f"{bids_id}_sessions.tsv",
            sep="\t",
            index=False,
            encoding="utf8",
        )


def _complete_examination_dates(
    rid: int,
    clinical_data_dir: Path,
    session_id: Optional[str] = None,
    examination_date: Optional[str] = None,
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

    for csv in _get_csv_files_for_alternative_exam_date(clinical_data_dir):
        csv_data = pd.read_csv(csv, low_memory=False)
        csv_data["SESSION"] = csv_data.VISCODE.apply(lambda x: viscode_to_session(x))
        exam_date = csv_data[(csv_data.RID == rid) & (csv_data.SESSION == session_id)]
        if not exam_date.empty and exam_date.iloc[0].EXAMDATE != "-4":
            return exam_date.iloc[0].EXAMDATE
    return None


def _get_csv_files_for_alternative_exam_date(
    clinical_data_dir: Path,
) -> Iterator[Path]:
    """Return a list of paths to CSV files in which an alternative exam date could be found."""

    for pattern in (
        "aibl_mri3meta_*.csv",
        "aibl_mrimeta_*.csv",
        "aibl_cdr_*.csv",
        "aibl_flutemeta_*.csv",
        "aibl_mmse_*.csv",
        "aibl_pibmeta_*.csv",
    ):
        try:
            yield next(clinical_data_dir.glob(pattern))
        except StopIteration:
            continue


# todo : test
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

    study = StudyName.AIBL.value

    specifications = _load_specifications(clinical_specifications_folder, "scans.tsv")
    scans_fields = specifications[study]
    field_location = specifications[f"{study} location"]
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
    scans_dict = create_scans_dict(
        clinical_data_dir,
        clinical_specifications_folder,
        bids_ids,
        "RID",
        "VISCODE",
        ses_dict,
    )
    write_scans_tsv(input_path, bids_ids, scans_dict)


# todo : test
def create_scans_dict(
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
    bids_ids: list[str],
    name_column_ids: str,
    name_column_ses: str,
    ses_dict: dict,
) -> pd.DataFrame:
    """[summary].

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the directory where the clinical data are stored.

    clinical_specifications_folder : Path
        The path to the folder containing the clinical specification files.

    bids_ids : list of str
        A list of bids ids.

    name_column_ids : str
        The name of the column where the subject id is contained.

    name_column_ses : str
        The name of the column where the viscode of the session is contained.

    ses_dict : dict
        The links the session id to the viscode of the session.

    Returns
    -------
    pd.DataFrame :
        A pandas DataFrame that contains the scans information for all sessions of all participants.
    """
    import datetime
    import os

    from clinica.utils.pet import Tracer
    from clinica.utils.stream import cprint

    scans_dict = {}
    prev_file = ""
    prev_sheet = ""

    study = StudyName.AIBL.value

    # Init the dictionary with the subject ids
    for bids_id in bids_ids:
        scans_dict[bids_id] = dict()
        for session_id in {"ses-" + key for key in ses_dict[bids_id].keys()}:
            scans_dict[bids_id][session_id] = {
                "T1/DWI/fMRI/FMAP": {},
                Tracer.PIB: {},
                Tracer.AV45: {},
                Tracer.FMM: {},
                Tracer.FDG: {},
            }

    scans_specs = pd.read_csv(clinical_specifications_folder / "scans.tsv", sep="\t")
    fields_dataset = []
    fields_location = []
    fields_bids = []
    fields_mod = []

    # Extract the fields available and the corresponding bids name, location and type
    for i in range(0, len(scans_specs[study])):
        field = scans_specs[study][i]
        if not pd.isnull(field):
            fields_dataset.append(field)
            fields_bids.append(scans_specs["BIDS CLINICA"][i])
            fields_location.append(scans_specs[f"{study} location"][i])
            fields_mod.append(scans_specs["Modalities related"][i])

    # For each field available extract the original name, extract from the file all the values and fill a data structure
    for i in range(0, len(fields_dataset)):
        # Location is composed by file/sheet
        location = fields_location[i].split("/")
        file_name = location[0]
        sheet = location[1] if len(location) > 1 else ""
        # Check if the file to read is already opened
        if file_name == prev_file and sheet == prev_sheet:
            pass
        else:
            file_ext = os.path.splitext(file_name)[1]
            files_to_read = [f for f in clinical_data_dir.glob(file_name)]
            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(files_to_read[0], sheet_name=sheet)
            elif file_ext == ".csv":
                file_path = files_to_read[0]

                # Fix for malformed flutemeta file in AIBL (see #796).
                # Some flutemeta lines contain a non-coded string value at the second-to-last position. This value
                # contains a comma which adds an extra column and shifts the remaining values to the right. In this
                # case, we just remove the erroneous content and replace it with -4 which AIBL uses as n/a value.
                on_bad_lines = lambda x: "error"  # noqa
                if "flutemeta" in file_path.name:
                    on_bad_lines = lambda bad_line: bad_line[:-3] + [-4, bad_line[-1]]  # noqa
                file_to_read = pd.read_csv(
                    file_path,
                    sep=",",
                    engine="python",
                    on_bad_lines=on_bad_lines,
                )
            prev_file = file_name
            prev_sheet = sheet

        for bids_id in bids_ids:
            original_id = bids_id.replace(
                f"sub-{study}", ""
            )  # todo : when refactoring function use BIDS id modelisation
            for session_name in {"ses-" + key for key in ses_dict[bids_id].keys()}:
                # When comparing sessions, remove the "-ses" prefix IF it exists
                row_to_extract = file_to_read[
                    (file_to_read[name_column_ids] == int(original_id))
                    & (
                        list(
                            filter(
                                None, file_to_read[name_column_ses].str.split("ses-")
                            )
                        )[0][0]
                        == ses_dict[bids_id][
                            list(filter(None, session_name.split("ses-")))[0]
                        ]
                    )
                ].index.tolist()
                if len(row_to_extract) > 0:
                    row_to_extract = row_to_extract[0]
                    # Fill the dictionary with all the information
                    value = file_to_read.iloc[row_to_extract][fields_dataset[i]]

                    # Deal with special format in AIBL
                    if value == "-4":
                        value = "n/a"
                    elif fields_bids[i] == "acq_time":
                        date_obj = datetime.datetime.strptime(value, "%m/%d/%Y")
                        value = date_obj.strftime("%Y-%m-%dT%H:%M:%S")

                    scans_dict[bids_id][session_name][fields_mod[i]][
                        fields_bids[i]
                    ] = value
                else:
                    cprint(
                        f"Scans information for {bids_id} {session_name} not found.",
                        lvl="info",
                    )
                    scans_dict[bids_id][session_name][fields_mod[i]][
                        fields_bids[i]
                    ] = "n/a"

    return scans_dict


# todo : test
def write_scans_tsv(
    bids_dir: Path, participant_ids: List[str], scans_dict: dict
) -> None:
    """Write the scans dict into TSV files.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    participant_ids : List[str]
        List of participant ids for which to write the scans TSV files.

    scans_dict : dict
        Dictionary containing scans metadata.

        .. note::
            This is the output of the function
            `clinica.iotools.bids_utils.create_scans_dict`.

    See also
    --------
    write_sessions_tsv
    """

    from clinica.iotools.bids_utils import _get_pet_tracer_from_filename

    supported_modalities = ("anat", "dwi", "func", "pet")

    for sub in participant_ids:
        for session_path in (bids_dir / sub).glob("ses-*"):
            scans_df = pd.DataFrame()
            tsv_file = (
                bids_dir
                / sub
                / session_path.name
                / f"{sub}_{session_path.name}_scans.tsv"
            )
            tsv_file.unlink(missing_ok=True)

            for mod in (bids_dir / sub / session_path.name).glob("*"):
                if mod.name in supported_modalities:
                    for file in [
                        file for file in mod.iterdir() if mod.suffix != ".json"
                    ]:
                        f_type = (
                            "T1/DWI/fMRI/FMAP"
                            if mod.name in ("anat", "dwi", "func")
                            else _get_pet_tracer_from_filename(file.name).value
                        )
                        row_to_append = pd.DataFrame(
                            scans_dict[sub][session_path.name][f_type], index=[0]
                        )
                        row_to_append.insert(
                            0, "filename", str(Path(mod.name) / Path(file.name))
                        )
                        scans_df = pd.concat([scans_df, row_to_append])
            scans_df = scans_df.set_index("filename").fillna("n/a")
            scans_df.to_csv(tsv_file, sep="\t", encoding="utf8")
