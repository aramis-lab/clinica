"""Methods used by BIDS converters."""

from os import PathLike
from pathlib import Path
from typing import BinaryIO, List, Optional, Union

from pandas import DataFrame

SUPPORTED_DATASETS = {
    "ADNI",
    "CLINAD",
    "PREVDEMALS",
    "INSIGHT",
    "OASIS",
    "OASIS3",
    "AIBL",
}

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
# @ToDo:test this function
def create_participants_df(
    study_name,
    clinical_spec_path,
    clinical_data_dir,
    bids_ids,
    delete_non_bids_info=True,
):
    """Create the file participants.tsv.

    Args:
        study_name: name of the study (Ex. ADNI)
        clinical_spec_path: path to the clinical file
        clinical_data_dir: path to the directory where the clinical data are stored
        bids_ids: list of bids ids
        delete_non_bids_info: if True delete all the rows of the subjects that
        are not available in the BIDS dataset

    Returns: a pandas dataframe that contains the participants data
    """
    import os
    from os import path

    import numpy as np
    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv
    from clinica.utils.stream import cprint

    fields_bids = ["participant_id"]
    prev_location = ""
    index_to_drop = []
    subjects_to_drop = []
    location_name = study_name + " location"

    # Load the data from the clincal specification file
    participants_specs = pd.read_csv(clinical_spec_path + "_participant.tsv", sep="\t")
    participant_fields_db = participants_specs[study_name]
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
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""
            # Check if the file to open for a certain field is the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

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

    if study_name == "ADNI" or study_name == "AIBL":
        # ADNImerge contains one row for each visits so there are duplicates
        participant_df = participant_df.drop_duplicates(
            subset=["alternative_id_1"], keep="first"
        )
    elif study_name == "OASIS":
        # OASIS provides several MRI for the same session
        participant_df = participant_df[
            ~participant_df.alternative_id_1.str.endswith("_MR2")
        ]
    participant_df.reset_index(inplace=True, drop=True)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        if study_name == "OASIS":
            value = (participant_df["alternative_id_1"][i].split("_"))[1]
        elif study_name == "OASIS3":
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


def create_sessions_dict_OASIS(
    clinical_data_dir,
    bids_dir,
    study_name,
    clinical_spec_path,
    bids_ids,
    name_column_ids,
    subj_to_remove=[],
    participants_df=None,
):
    """Extract the information regarding the sessions and store them in a dictionary (session M00 only).

    Args:
        clinical_data_dir: path to the input folder
        study_name: name of the study (Ex: ADNI)
        clinical_spec_path:  path to the clinical file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject ids are stored
        subj_to_remove: subjects to remove
        participants_df: a pandas dataframe that contains the participants data (required for OASIS3 only)
    """
    import os
    from os import path

    import numpy as np
    import pandas as pd

    from clinica.utils.stream import cprint

    # Load data
    location = study_name + " location"
    sessions = pd.read_csv(clinical_spec_path + "_sessions.tsv", sep="\t")
    sessions_fields = sessions[study_name]
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

            file_to_read_path = path.join(clinical_data_dir, location)
            file_ext = os.path.splitext(location)[1]

            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
            elif file_ext == ".csv":
                file_to_read = pd.read_csv(file_to_read_path)

            for r in range(0, len(file_to_read.values)):
                # Extracts the subject ids columns from the dataframe
                subj_id = file_to_read.iloc[r][name_column_ids]
                if hasattr(subj_id, "dtype"):
                    if subj_id.dtype == np.int64:
                        subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)
                if study_name == "OASIS":
                    subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3] + subj_id[5:9])
                if study_name == "OASIS3":
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

                    subj_dir = path.join(
                        bids_dir,
                        subj_bids,
                    )
                    session_names = get_bids_sess_list(subj_dir)
                    for s in session_names:
                        s_name = s.replace("ses-", "")
                        if study_name != "OASIS3":
                            row = file_to_read.iloc[r]
                        else:
                            row = file_to_read[
                                file_to_read["MR ID"].str.startswith(subj_id)
                                & file_to_read["MR ID"].str.endswith(s_name)
                            ].iloc[0]
                        if subj_bids not in sessions_dict:
                            sessions_dict.update({subj_bids: {}})
                        if s_name not in sessions_dict[subj_bids].keys():
                            sessions_dict[subj_bids].update({s_name: {"session_id": s}})
                        (sessions_dict[subj_bids][s_name]).update(
                            {sessions_fields_bids[i]: row[sessions_fields[i]]}
                        )
                        # Calculate the difference in months for OASIS3 only
                        if study_name == "OASIS3" and sessions_fields_bids[i] == "age":
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
    clinical_data_dir,
    study_name,
    clinic_specs_path,
    bids_ids,
    name_column_ids,
    name_column_ses,
    ses_dict,
):
    """[summary].

    Args:
        clinical_data_dir: path to the directory where the clinical data are stored
        study_name: name of the study (Ex ADNI)
        clinic_specs_path: path to the clinical specification file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject id is contained
        name_column_ses: name of the column where the viscode of the session is contained
        ses_dict: links the session id to the viscode of the session.

    Returns: a pandas DataFrame that contains the scans information for all sessions of all participants.
    """
    import datetime
    import glob
    from os import path

    import pandas as pd

    from clinica.utils.pet import Tracer
    from clinica.utils.stream import cprint

    scans_dict = {}
    prev_file = ""
    prev_sheet = ""

    if study_name not in SUPPORTED_DATASETS:
        raise ValueError(
            f"Dataset {study_name} is not supported. "
            f"Supported datasets are: {SUPPORTED_DATASETS}."
        )

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

    scans_specs = pd.read_csv(clinic_specs_path + "_scans.tsv", sep="\t")
    fields_dataset = []
    fields_location = []
    fields_bids = []
    fields_mod = []

    # Extract the fields available and the corresponding bids name, location and type
    for i in range(0, len(scans_specs[study_name])):
        field = scans_specs[study_name][i]
        if not pd.isnull(field):
            fields_dataset.append(field)
            fields_bids.append(scans_specs["BIDS CLINICA"][i])
            fields_location.append(scans_specs[study_name + " location"][i])
            fields_mod.append(scans_specs["Modalities related"][i])

    # For each field available extract the original name, extract from the file all the values and fill a data structure
    for i in range(0, len(fields_dataset)):
        # Location is composed by file/sheet
        location = fields_location[i].split("/")
        file_name = location[0]
        if len(location) > 1:
            sheet = location[1]
        else:
            sheet = ""
        # Check if the file to read is already opened
        if file_name == prev_file and sheet == prev_sheet:
            pass
        else:
            file_to_read_path = path.join(clinical_data_dir, file_name)
            file_ext = path.splitext(file_name)[1]
            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(
                    glob.glob(file_to_read_path)[0], sheet_name=sheet
                )
            elif file_ext == ".csv":
                file_path = glob.glob(file_to_read_path)[0]

                # Fix for malformed flutemeta file in AIBL (see #796).
                # Some flutemeta lines contain a non-coded string value at the second-to-last position. This value
                # contains a comma which adds an extra column and shifts the remaining values to the right. In this
                # case, we just remove the erroneous content and replace it with -4 which AIBL uses as n/a value.
                on_bad_lines = (
                    lambda bad_line: bad_line[:-3] + [-4, bad_line[-1]]
                    if "flutemeta" in file_path and study_name == "AIBL"
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
            original_id = bids_id.replace("sub-" + study_name, "")
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

                    if study_name == "AIBL":  # Deal with special format in AIBL
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
    study_name: str,
    bids_dir: Union[str, Path],
    bids_version: Optional[str] = None,
) -> None:
    """Write `dataset_description.json` at the root of the BIDS directory."""
    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    if bids_version:
        bids_desc = BIDSDatasetDescription(name=study_name, bids_version=bids_version)
    else:
        bids_desc = BIDSDatasetDescription(name=study_name)
    with open(Path(bids_dir) / "dataset_description.json", "w") as f:
        bids_desc.write(to=f)


def _write_readme(
    study_name: str,
    data_dict: dict,
    bids_dir: Union[str, Path],
) -> None:
    """Write `README` at the root of the BIDS directory."""
    from clinica.iotools.bids_readme import BIDSReadme

    bids_readme = BIDSReadme(
        name=study_name,
        link=data_dict["link"],
        description=data_dict["desc"],
    )

    with open(Path(bids_dir) / "README", "w") as f:
        bids_readme.write(to=f)


def _write_bids_validator_config(bids_dir: Union[str, Path]) -> None:
    """Write `.bids-validator-config.json` at the root of the BIDS directory."""
    import json

    with open(Path(bids_dir) / ".bids-validator-config.json", "w") as f:
        json.dump(BIDS_VALIDATOR_CONFIG, f, skipkeys=True, indent=4)


def _write_bidsignore(bids_dir: Union[str, Path]) -> None:
    """Write `.bidsignore` file at the root of the BIDS directory."""
    with open(Path(bids_dir) / ".bidsignore", "w") as f:
        # pet/ is necessary until PET is added to BIDS standard
        f.write("\n".join(["swi/\n"]))
        f.write("\n".join(["conversion_info/"]))


def write_modality_agnostic_files(
    study_name: str,
    readme_data: dict,
    bids_dir: Union[str, Path],
    bids_version: Optional[str] = None,
) -> None:
    """
    Write the files README, dataset_description.json, .bidsignore and .bids-validator-config.json
    at the root of the BIDS directory.

    Args:
        study_name: name of the study (Ex ADNI)
        readme_data: dictionary containing the data specific to the dataset to write in the readme
        bids_dir: path to the bids directory
        bids_version: BIDS version if different from the version supported by Clinica.
    """
    _write_bids_dataset_description(study_name, bids_dir, bids_version)
    _write_readme(study_name, readme_data, bids_dir)
    _write_bids_validator_config(bids_dir)
    _write_bidsignore(bids_dir)


def write_sessions_tsv(bids_dir: Union[str, Path], sessions_dict: dict) -> None:
    """Create <participant_id>_sessions.tsv files.

    Basically writes the content of the function
    `clinica.iotools.bids_utils.create_sessions_dict` in several TSV files
    following the BIDS specification.

    Parameters
    ----------
    bids_dir : str
        Path to the BIDS directory.

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
    from pathlib import Path

    import pandas as pd

    bids_dir = Path(bids_dir)

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


def _get_pet_tracer_from_filename(filename: str) -> str:
    """Return the PET tracer from the provided filename.

    Parameters
    ----------
    filename : str
        The filename from which to extract the PET tracer.

    Returns
    -------
    tracer : str
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

    return tracer


def write_scans_tsv(
    bids_dir: Union[str, Path], participant_ids: List[str], scans_dict: dict
) -> None:
    """Write the scans dict into TSV files.

    Parameters
    ----------
    bids_dir : str
        Path to the BIDS directory.

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
    from pathlib import Path

    import pandas as pd

    bids_dir = Path(bids_dir)
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
                            else _get_pet_tracer_from_filename(file.name)
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


def get_bids_subjs_list(bids_path: str) -> List[str]:
    """Given a BIDS compliant dataset, return the list of all the subjects available.

    Parameters
    ----------
    bids_path : str
        Path to the BIDS directory.

    Returns
    -------
    List[str] :
        List of subject IDs available in this BIDS dataset.

    See also
    --------
    get_bids_sess_list
    get_bids_subjs_paths
    """
    from pathlib import Path

    return [str(d.name) for d in Path(bids_path).glob("sub-*") if d.is_dir()]


def get_bids_sess_list(subj_path: str) -> List[str]:
    """Given a path to a subject's folder, this function returns the
    list of sessions available.

    Parameters
    ----------
    subj_path : str
        Path to the subject folder for which to list the sessions.

    Returns
    -------
    List[str] :
        The list of session names for this subject.

    See also
    --------
    get_bids_subjs_list
    get_bids_subjs_paths
    """
    from pathlib import Path

    return [str(d.name) for d in Path(subj_path).glob("ses-*") if d.is_dir()]


def get_bids_subjs_paths(bids_path: str) -> List[str]:
    """Given a BIDS compliant dataset, returns the list of all paths to the subjects folders.

    Parameters
    ----------
    bids_path : str
        Path to the BIDS directory.

    Returns
    -------
    List[str] :
        List of paths to the subjects folders.

    See also
    --------
    get_bids_subjs_list
    get_bids_sess_list
    """
    from pathlib import Path

    return [str(d) for d in Path(bids_path).glob("sub-*") if d.is_dir()]


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


def json_from_dcm(dcm_dir: str, json_path: str) -> None:
    """Writes descriptive JSON file from DICOM header.

    Parameters
    ----------
    dcm_dir : str
        Path to the DICOM directory.

    json_path : str
        Path to the output JSON file.
    """
    import json
    from glob import glob
    from os import path

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
        dcm_path = glob(path.join(dcm_dir, "*.dcm"))[0]
        ds = dcmread(dcm_path)
        json_dict = {
            key: ds.get(tag).value for key, tag in fields_dict.items() if tag in ds
        }
        json = json.dumps(json_dict, skipkeys=True, indent=4)
        with open(json_path, "w") as f:
            f.write(json)
    except IndexError:
        cprint(msg=f"No DICOM found at {dcm_dir}", lvl="warning")


def _build_dcm2niix_command(
    input_dir: str,
    output_dir: str,
    output_fmt: str,
    compress: bool = False,
    bids_sidecar: bool = True,
) -> list:
    """Generate the dcm2niix command from user inputs."""
    command = ["dcm2niix", "-w", "0", "-f", output_fmt, "-o", output_dir]
    command += ["-9", "-z", "y"] if compress else ["-z", "n"]
    command += ["-b", "y", "-ba", "y"] if bids_sidecar else ["-b", "n"]
    command += [input_dir]

    return command


def run_dcm2niix(
    input_dir: str,
    output_dir: str,
    output_fmt: str,
    compress: bool = False,
    bids_sidecar: bool = True,
) -> bool:
    """Runs the `dcm2niix` command using a subprocess.

    Parameters
    ----------
    input_dir : str
        Input folder.

    output_dir : str
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

    cprint(f"Attempting to convert {output_fmt}.", lvl="info")
    command = _build_dcm2niix_command(
        input_dir, output_dir, output_fmt, compress, bids_sidecar
    )
    completed_process = subprocess.run(command, capture_output=True)

    if completed_process.returncode != 0:
        cprint(
            msg=(
                "DICOM to BIDS conversion with dcm2niix failed:\n"
                f"command: {command}\n"
                f"{completed_process.stdout.decode('utf-8')}"
            ),
            lvl="warning",
        )
        return False
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


def parse_description(filepath: PathLike, start_line: int, end_line: int) -> str:
    """Parse the description of the dataset from the readme in the documentation.

    Parameters
    ----------
    filepath : PathLike
        Path to the readme file from which to extract the description.

    start_line : int
        Line number where the description starts.

    end_line : int
        Line number where the description ends.

    Returns
    -------
    str :
        The description extracted from the readme file.
    """
    with open(filepath, "r") as fp:
        txt = fp.readlines()
    return "\n".join(txt[start_line:end_line])


def parse_url(filepath: PathLike) -> List[str]:
    """Parse the URLs from the readme in the documentation.

    Parameters
    ----------
    filepath : PathLike
        Path to the readme file from which to extract the URLs.

    Returns
    -------
    List[str] :
        The list of found URLs.
    """
    import re

    with open(filepath, "r") as fp:
        txt = fp.readlines()
    txt = "".join(txt)
    url_extract_pattern = "https?:\\/\\/(?:www\\.)?[-a-zA-Z0-9@:%._\\+~#=]{1,256}\\.[a-zA-Z0-9()]{1,6}\\b(?:[-a-zA-Z0-9()@:%_\\+.~#?&\\/=]*)"
    return re.findall(url_extract_pattern, txt)


def write_to_tsv(df: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    """Save dataframe as a BIDS-compliant TSV file.

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame to write to TSV file.

    buffer : PathLike or BinaryIO
        Either a file or a stream to write the DataFrame to.
    """
    df.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")
