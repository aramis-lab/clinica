"""Methods used by BIDS converters."""

import json
import os
from enum import Enum
from pathlib import Path
from typing import BinaryIO, List, Optional, Union

import pandas as pd

from clinica.utils.pet import Tracer


class StudyName(str, Enum):
    """Studies supported by the converters of Clinica."""

    ADNI = "ADNI"
    AIBL = "AIBL"
    GENFI = "GENFI"
    HABS = "HABS"
    NIFD = "NIFD"
    OASIS = "OASIS"
    OASIS3 = "OASIS3"
    UKB = "UKB"


BIDS_VALIDATOR_CONFIG = {
    "ignore": [
        # Possibly dcm2nii(x) errors
        "NIFTI_UNIT",
        "INCONSISTENT_PARAMETERS",
        # fMRI-specific errors
        "SLICE_TIMING_NOT_DEFINED",
        "NIFTI_PIXDIM4",
        "BOLD_NOT_4D",
        "REPETITION_TIME_MUST_DEFINE",
        "TASK_NAME_MUST_DEFINE",
        # Won't fix errors
        "MISSING_SESSION",  # Allows subjects to have different sessions
        "INCONSISTENT_SUBJECTS",  # Allows subjects to have different modalities
        "SCANS_FILENAME_NOT_MATCH_DATASET",  # Necessary until PET is added to BIDS standard
        "CUSTOM_COLUMN_WITHOUT_DESCRIPTION",  # We won't create these JSON files as clinical description
        # is already done in TSV files of clinica.
        "NO_AUTHORS",  # Optional field in dataset_description.json
    ],
    "warn": [],
    "error": [],
    "ignoredFiles": [],
}


# -- Methods for the clinical data --
def create_participants_df(
    study_name: StudyName,
    clinical_specifications_folder: Path,
    clinical_data_dir: Path,
    bids_ids: list[str],
    delete_non_bids_info: bool = True,
) -> pd.DataFrame:
    """Create the file participants.tsv.

    Parameters
    ----------
    study_name : StudyName
        The name of the study (Ex. ADNI).

    clinical_specifications_folder : Path
        The path to the clinical file.

    clinical_data_dir : Path
        The path to the directory where the clinical data are stored.

    bids_ids : list of str
        The list of bids ids.

    delete_non_bids_info : bool, optional
        If True delete all the rows of the subjects that are not available in the BIDS dataset.
        Default=True.

    Returns
    -------
    pd.DataFrame :
        A pandas dataframe that contains the participants data.
    """
    import numpy as np

    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv
    from clinica.utils.stream import cprint

    fields_bids = ["participant_id"]
    prev_location = ""
    index_to_drop = []
    subjects_to_drop = []
    study_name = StudyName(study_name)
    location_name = f"{study_name.value} location"

    participants_specs = pd.read_csv(
        clinical_specifications_folder / "participant.tsv", sep="\t"
    )
    participant_fields_db = participants_specs[study_name.value]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs["BIDS CLINICA"]

    # Extract the list of the available BIDS fields for the dataset
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])

    # Init the dataframe that will be saved in the file participants.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

    for i in range(0, len(participant_fields_db)):
        # If a field not empty is found
        if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            # If a sheet is available
            sheet = tmp[1] if len(tmp) > 1 else ""
            # Check if the file to open for a certain field is the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = clinical_data_dir / location
                if file_ext == ".xlsx":
                    file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
                elif file_ext == ".csv":
                    file_to_read = load_clinical_csv(
                        clinical_data_dir, location.split(".")[0]
                    )
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                value_to_read = file_to_read[participant_fields_db[i]]
                if participant_fields_bids[i] == "alternative_id_1" and (
                    value_to_read.dtype == np.float64 or value_to_read.dtype == np.int64
                ):
                    if not pd.isnull(file_to_read.at[j, participant_fields_db[i]]):
                        value_to_append = str(
                            file_to_read.at[j, participant_fields_db[i]]
                        ).rstrip(".0")
                    else:
                        value_to_append = np.NaN
                else:
                    value_to_append = file_to_read.at[j, participant_fields_db[i]]
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    if study_name == StudyName.ADNI or study_name == StudyName.AIBL:
        # ADNImerge contains one row for each visits so there are duplicates
        participant_df = participant_df.drop_duplicates(
            subset=["alternative_id_1"], keep="first"
        )
    elif study_name == StudyName.OASIS:
        # OASIS provides several MRI for the same session
        participant_df = participant_df[
            ~participant_df.alternative_id_1.str.endswith("_MR2")
        ]
    participant_df.reset_index(inplace=True, drop=True)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        if study_name == StudyName.OASIS:
            value = (participant_df["alternative_id_1"][i].split("_"))[1]
        elif study_name == StudyName.OASIS3:
            value = participant_df["alternative_id_1"][i].replace("OAS3", "")
        else:
            value = remove_space_and_symbols(participant_df["alternative_id_1"][i])

        bids_id = [s for s in bids_ids if value in s]

        if len(bids_id) == 0:
            index_to_drop.append(i)
            subjects_to_drop.append(value)
        else:
            participant_df.at[i, "participant_id"] = bids_id[0]

    if len(subjects_to_drop) > 0:
        cprint(
            msg=(
                "The following subjects of dataset directory were not found in your BIDS folder :\n"
                + ", ".join(subjects_to_drop)
            ),
            lvl="info",
        )
    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.drop(index_to_drop)

    participant_df = participant_df.fillna("n/a")

    return participant_df


def create_sessions_dict_oasis(
    clinical_data_dir: Path,
    bids_dir: Path,
    study_name: StudyName,
    clinical_specifications_folder: Path,
    bids_ids: list[str],
    name_column_ids: str,
    subj_to_remove: Optional[list[str]] = None,
    participants_df: Optional[pd.DataFrame] = None,
) -> dict:
    """Extract the information regarding the sessions and store them in a dictionary (session M00 only).

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the input folder.

    bids_dir : Path
        The path to the BIDS directory.

    study_name : StudyName
        The name of the study (Ex: ADNI).

    clinical_specifications_folder : Path
        The path to the clinical file.

    bids_ids : list of str
        The list of bids ids.

    name_column_ids : str
        The name of the column where the subject ids are stored.

    subj_to_remove : list of str, optional
        The list of subject IDs to remove.

    participants_df : pd.DataFrame, optional
        A pandas dataframe that contains the participants data (required for OASIS3 only).

    Returns
    -------
    dict :
        Session dict.
    """
    import numpy as np

    from clinica.utils.stream import cprint

    subj_to_remove = subj_to_remove or []
    location = f"{study_name.value} location"
    sessions = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")
    sessions_fields = sessions[study_name.value]
    field_location = sessions[location]
    sessions_fields_bids = sessions["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""

            file_to_read_path = clinical_data_dir / location
            file_ext = os.path.splitext(location)[1]
            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
            elif file_ext == ".csv":
                file_to_read = pd.read_csv(file_to_read_path)
            else:
                raise ValueError(
                    f"Unknown file extension {file_ext}. Expecting either .xlsx or .csv."
                )

            for r in range(0, len(file_to_read.values)):
                # Extracts the subject ids columns from the dataframe
                subj_id = file_to_read.iloc[r][name_column_ids]
                if hasattr(subj_id, "dtype"):
                    if subj_id.dtype == np.int64:
                        subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)
                if study_name == StudyName.OASIS:
                    subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3] + subj_id[5:9])
                if study_name == StudyName.OASIS3:
                    subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3:])

                # Extract the corresponding BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if len(subj_bids) == 0:
                    # If the subject is not an excluded one
                    if subj_id not in subj_to_remove:
                        cprint(
                            f"{sessions_fields[i]} for {subj_id} not found in the BIDS converted.",
                            "info",
                        )
                else:
                    subj_bids = subj_bids[0]
                    subj_dir = bids_dir / subj_bids
                    session_names = get_bids_sess_list(subj_dir)
                    for s in session_names:
                        s_name = s.replace("ses-", "")
                        if study_name == StudyName.OASIS3:
                            row = file_to_read[
                                file_to_read["MR ID"].str.startswith(subj_id)
                                & file_to_read["MR ID"].str.endswith(s_name)
                            ].iloc[0]
                        else:
                            row = file_to_read.iloc[r]
                        if subj_bids not in sessions_dict:
                            sessions_dict.update({subj_bids: {}})
                        if s_name not in sessions_dict[subj_bids].keys():
                            sessions_dict[subj_bids].update({s_name: {"session_id": s}})
                        (sessions_dict[subj_bids][s_name]).update(
                            {sessions_fields_bids[i]: row[sessions_fields[i]]}
                        )
                        # Calculate the difference in months for OASIS3 only
                        if (
                            study_name == StudyName.OASIS3
                            and sessions_fields_bids[i] == "age"
                        ):
                            diff_years = (
                                float(sessions_dict[subj_bids][s_name]["age"])
                                - participants_df[
                                    participants_df["participant_id"] == subj_bids
                                ]["age_bl"]
                            )
                            (sessions_dict[subj_bids][s_name]).update(
                                {"diff_months": round(float(diff_years) * 12)}
                            )

    return sessions_dict


def create_scans_dict(
    clinical_data_dir: Path,
    study_name: StudyName,
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

    study_name : StudyName
        The name of the study (Ex ADNI).

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

    from clinica.utils.pet import Tracer
    from clinica.utils.stream import cprint

    scans_dict = {}
    prev_file = ""
    prev_sheet = ""

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
    for i in range(0, len(scans_specs[study_name.value])):
        field = scans_specs[study_name.value][i]
        if not pd.isnull(field):
            fields_dataset.append(field)
            fields_bids.append(scans_specs["BIDS CLINICA"][i])
            fields_location.append(scans_specs[f"{study_name.value} location"][i])
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
                on_bad_lines = (  # noqa: E731
                    lambda bad_line: bad_line[:-3] + [-4, bad_line[-1]]
                    if "flutemeta" in file_path and study_name == StudyName.AIBL
                    else "error"
                )

                file_to_read = pd.read_csv(
                    file_path,
                    sep=",",
                    engine="python",
                    on_bad_lines=on_bad_lines,
                )
            prev_file = file_name
            prev_sheet = sheet

        for bids_id in bids_ids:
            original_id = bids_id.replace(f"sub-{study_name.value}", "")
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

                    if study_name == StudyName.AIBL:  # Deal with special format in AIBL
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


def _write_bids_dataset_description(
    study_name: StudyName,
    bids_dir: Path,
    bids_version: Optional[str] = None,
) -> None:
    """Write `dataset_description.json` at the root of the BIDS directory."""
    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    if bids_version:
        bids_desc = BIDSDatasetDescription(name=study_name, bids_version=bids_version)
    else:
        bids_desc = BIDSDatasetDescription(name=study_name)
    with open(bids_dir / "dataset_description.json", "w") as f:
        bids_desc.write(to=f)


def _write_readme(
    study_name: StudyName,
    data_dict: dict,
    bids_dir: Path,
) -> None:
    """Write `README` at the root of the BIDS directory."""
    from clinica.iotools.bids_readme import BIDSReadme

    bids_readme = BIDSReadme(
        name=study_name,
        link=data_dict["link"],
        description=data_dict["desc"],
    )
    with open(bids_dir / "README", "w") as f:
        bids_readme.write(to=f)


def _write_bids_validator_config(bids_dir: Path) -> None:
    """Write `.bids-validator-config.json` at the root of the BIDS directory."""
    with open(bids_dir / ".bids-validator-config.json", "w") as f:
        json.dump(BIDS_VALIDATOR_CONFIG, f, skipkeys=True, indent=4)


def _write_bidsignore(bids_dir: Path) -> None:
    """Write `.bidsignore` file at the root of the BIDS directory."""
    with open(bids_dir / ".bidsignore", "w") as f:
        # pet/ is necessary until PET is added to BIDS standard
        f.write("\n".join(["swi/\n"]))
        f.write("\n".join(["conversion_info/"]))


def write_modality_agnostic_files(
    study_name: StudyName,
    readme_data: dict,
    bids_dir: Path,
    bids_version: Optional[str] = None,
) -> None:
    """
    Write the files README, dataset_description.json, .bidsignore and .bids-validator-config.json
    at the root of the BIDS directory.

    Parameters
    ----------
    study_name : StudyName
        The name of the study (Ex ADNI).

    readme_data : dict
        A dictionary containing the data specific to the dataset to write in the readme.

    bids_dir : Path
        The path to the bids directory.

    bids_version : str, optional
        The BIDS version if different from the version supported by Clinica.
    """
    _write_bids_dataset_description(study_name, bids_dir, bids_version)
    _write_readme(study_name, readme_data, bids_dir)
    _write_bids_validator_config(bids_dir)
    _write_bidsignore(bids_dir)


def write_sessions_tsv(bids_dir: Path, sessions_dict: dict) -> None:
    """Create <participant_id>_sessions.tsv files.

    Basically writes the content of the function
    `clinica.iotools.bids_utils.create_sessions_dict` in several TSV files
    following the BIDS specification.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    sessions_dict : dict
        Dictionary containing sessions metadata.

        .. note::
            This is the output of the function
            `clinica.iotools.bids_utils.create_sessions_dict`.

    See also
    --------
    create_sessions_dict
    write_scans_tsv
    """
    for subject_path in bids_dir.glob("sub-*"):
        if subject_path.name in sessions_dict:
            session_df = pd.DataFrame.from_dict(
                sessions_dict[subject_path.name], orient="index"
            )
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
        else:
            print(f"No session data available for {subject_path}")
            session_df = pd.DataFrame(columns=["session_id"])
            session_df["session_id"] = pd.Series("M000")
        session_df = session_df.set_index("session_id").fillna("n/a")
        session_df.to_csv(
            subject_path / f"{subject_path.name}_sessions.tsv",
            sep="\t",
            encoding="utf8",
        )


def _get_pet_tracer_from_filename(filename: str) -> Tracer:
    """Return the PET tracer from the provided filename.

    Parameters
    ----------
    filename : str
        The filename from which to extract the PET tracer.

    Returns
    -------
    tracer : Tracer
        The PET tracer.

    Raises
    ------
    ValueError
        If no tracer found in the filename.
    """
    import re

    tracer = None
    for entity in ("trc", "acq"):
        m = re.search(rf"({entity}-[a-zA-Z0-9]+)", filename)
        if m:
            tracer = m.group(1)[4:].upper()
    if tracer is None:
        raise ValueError(
            f"Could not extract the PET tracer from the following file name {filename}."
        )

    return Tracer(tracer)


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


def get_bids_subjs_list(bids_path: Path) -> List[str]:
    """Given a BIDS compliant dataset, return the list of all the subjects available.

    Parameters
    ----------
    bids_path : Path
        The path to the BIDS directory.

    Returns
    -------
    List[str] :
        List of subject IDs available in this BIDS dataset.

    See also
    --------
    get_bids_sess_list
    get_bids_subjs_paths
    """
    return _filter_folder_names(bids_path, "sub-*")


def _filter_folder_names(folder: Path, pattern: str) -> List[str]:
    return [d.name for d in folder.glob(pattern) if d.is_dir()]


def get_bids_sess_list(subj_path: Path) -> List[str]:
    """Given a path to a subject's folder, this function returns the
    list of sessions available.

    Parameters
    ----------
    subj_path : Path
        The path to the subject folder for which to list the sessions.

    Returns
    -------
    List[str] :
        The list of session names for this subject.

    See also
    --------
    get_bids_subjs_list
    get_bids_subjs_paths
    """
    return _filter_folder_names(subj_path, "ses-*")


def get_bids_subjs_paths(bids_path: Path) -> List[Path]:
    """Given a BIDS compliant dataset, returns the list of all paths to the subjects folders.

    Parameters
    ----------
    bids_path : str
        Path to the BIDS directory.

    Returns
    -------
    List[Path] :
        List of paths to the subjects folders.

    See also
    --------
    get_bids_subjs_list
    get_bids_sess_list
    """
    return [d for d in bids_path.glob("sub-*") if d.is_dir()]


def remove_space_and_symbols(data: Union[str, List[str]]) -> Union[str, List[str]]:
    """Remove spaces, "-", and "_" characters from a string.

    If a list of strings is provided, this function will be called
    on each item of the list.

    Parameters
    ----------
    data : str or List[str]
        String or list of strings to be cleaned.

    Returns
    -------
    str or List[str] :
        Cleaned string(s).
    """
    import re

    if isinstance(data, list):
        return [remove_space_and_symbols(d) for d in data]
    return re.sub("[-_ ]", "", data)


def json_from_dcm(dcm_dir: Path, json_path: Path) -> None:
    """Writes descriptive JSON file from DICOM header.

    Parameters
    ----------
    dcm_dir : Path
        The path to the DICOM directory.

    json_path : Path
        The path to the output JSON file.
    """
    from pydicom import dcmread
    from pydicom.tag import Tag

    from clinica.utils.stream import cprint

    fields_dict = {
        "DeviceSerialNumber": Tag(("0018", "1000")),
        "Manufacturer": Tag(("0008", "0070")),
        "ManufacturersModelName": Tag(("0008", "1090")),
        "SoftwareVersions": Tag(("0018", "1020")),
        "BodyPart": Tag(("0018", "0015")),
        "Units": Tag(("0054", "1001")),
        # Institution
        "InstitutionName": Tag(("0008", "0080")),
        "InstitutionAddress": Tag(("0008", "0081")),
        "InstitutionalDepartmentName": Tag(("0008", "1040")),
        # MRI
        "MagneticFieldStrength": Tag(("0018", "0087")),
        # PET
        "InjectedRadioactivity": Tag(("0018", "1074")),
        "MolarActivity": Tag(("0018", "1077")),
        "InjectionStart": Tag(("0018", "1042")),
        "FrameDuration": Tag(("0018", "1242")),
    }

    try:
        dcm_path = [f for f in dcm_dir.glob("*.dcm")][0]
        ds = dcmread(dcm_path)
        json_dict = {
            key: ds.get(tag).value for key, tag in fields_dict.items() if tag in ds
        }
        with open(json_path, "w") as f:
            f.write(json.dumps(json_dict, skipkeys=True, indent=4))
    except IndexError:
        cprint(msg=f"No DICOM found at {dcm_dir}", lvl="warning")


def _build_dcm2niix_command(
    input_dir: Path,
    output_dir: Path,
    output_fmt: str,
    compress: bool = False,
    bids_sidecar: bool = True,
) -> list:
    """Generate the dcm2niix command from user inputs."""
    command = ["dcm2niix", "-w", "0", "-f", output_fmt, "-o", str(output_dir)]
    command += ["-9", "-z", "y"] if compress else ["-z", "n"]
    command += ["-b", "y", "-ba", "y"] if bids_sidecar else ["-b", "n"]
    command += [str(input_dir)]

    return command


def run_dcm2niix(
    input_dir: Path,
    output_dir: Path,
    output_fmt: str,
    compress: bool = False,
    bids_sidecar: bool = True,
) -> bool:
    """Runs the `dcm2niix` command using a subprocess.

    Parameters
    ----------
    input_dir : Path
        Input folder.

    output_dir : Path
        Output folder. This will be passed to the
        dcm2niix "-o" option.

    output_fmt : str
        Output format. This will be passed to the
        dcm2niix "-f" option.

    compress : bool, optional
        Whether to compress or not.
        Default=False.

    bids_sidecar : bool, optional
        Whether to generate a BIDS sidecar or not. Default=True.

    Returns
    -------
    bool :
        True if the conversion was successful, False otherwise.
    """
    import subprocess

    from clinica.utils.stream import cprint

    cprint(f"Attempting to convert {output_fmt}.", lvl="debug")
    command = _build_dcm2niix_command(
        input_dir, output_dir, output_fmt, compress, bids_sidecar
    )
    completed_process = subprocess.run(command, capture_output=True)

    if completed_process.returncode != 0:
        output_message = (
            completed_process.stdout.decode("utf-8")
            if completed_process.stdout is not None
            else ""
        )
        cprint(
            msg=(
                "The DICOM to BIDS conversion with dcm2niix failed:\n"
                f"command: {' '.join(command)}\n"
                f"{output_message}"
            ),
            lvl="warning",
        )
        return False
    cprint(
        msg=(
            f"The DICOM to BIDS conversion with dcm2niix was successful:\n"
            f"command: {' '.join(command)}\n"
        ),
        lvl="debug",
    )
    return True


def identify_modality(filename: str) -> Optional[str]:
    """Identifies the modality of a file given its name.

    Parameters
    ----------
    filename: str
        Input filename

    Returns
    -------
    Optional[str]:
        Modality or None if parsing uns
    """
    import numpy as np

    filename = filename.lower()
    if "dwi" in filename:
        return "dwi"
    if "t1" in filename:
        return "T1"
    if "t2" in filename:
        return "T2w"
    if "fieldmap" in filename:
        return "fieldmap"
    if "fmri" in filename:
        return "rsfmri"
    else:
        return np.nan


def write_to_tsv(df: pd.DataFrame, buffer: Union[Path, BinaryIO]) -> None:
    """Save dataframe as a BIDS-compliant TSV file.

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame to write to TSV file.

    buffer : PathLike or BinaryIO
        Either a file or a stream to write the DataFrame to.
    """
    df.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")
