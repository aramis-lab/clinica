"""Methods used by BIDS converters."""

import json
import os
import re
from pathlib import Path
from typing import BinaryIO, List, Optional, Union

import pandas as pd

from clinica.utils.filemanip import UserProvidedPath

from .study_models import StudyName

__all__ = [
    "create_participants_df",
    "get_subjects_from_bids_dataset",
    "get_sessions_from_bids_dataset",
    "get_paths_to_subjects_in_bids_dataset",
    "json_from_dcm",
    "run_dcm2niix",
    "write_to_tsv",
    "identify_modality",
    "write_modality_agnostic_files",
    "validate_input_path",
    "viscode_to_session",
    "load_clinical_csv",
]


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

    from clinica.utils.stream import cprint

    from .study_models import bids_id_factory

    fields_bids = ["participant_id"]
    prev_location = ""
    prev_sheet = 0
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
            sheet = tmp[1] if len(tmp) > 1 else 0
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
                    # Condition to handle ADNI modification of file APOERES.csv
                    # See issue https://github.com/aramis-lab/clinica/issues/1294
                    if study_name == StudyName.ADNI and location == "APOERES.csv":
                        if (
                            participant_fields_db[i] not in file_to_read.columns
                            and "GENOTYPE" in file_to_read.columns
                        ):
                            # Split the 'GENOTYPE' column into 'APGEN1' and 'APGEN2'
                            genotype = file_to_read["GENOTYPE"].str.split(
                                "/", expand=True
                            )
                            file_to_read = file_to_read.assign(
                                APGEN1=genotype[0], APGEN2=genotype[1]
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
        value = bids_id_factory(study_name).from_original_study_id(
            participant_df["alternative_id_1"][i]
        )
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


def _write_bids_dataset_description(
    study_name: StudyName,
    bids_dir: Path,
    bids_version: Optional[str] = None,
) -> None:
    """Write `dataset_description.json` at the root of the BIDS directory."""
    from .bids_dataset_description import BIDSDatasetDescription

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
    from .bids_readme import BIDSReadme

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


def get_subjects_from_bids_dataset(bids_path: Path) -> List[str]:
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


def get_sessions_from_bids_dataset(subj_path: Path) -> List[str]:
    """Given a path to a subject's folder, this function returns the list of sessions available.

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
    get_subjects_from_bids_dataset
    get_paths_to_subjects_in_bids_dataset
    """
    return _filter_folder_names(subj_path, "ses-*")


def get_paths_to_subjects_in_bids_dataset(bids_path: Path) -> List[Path]:
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
    get_subjects_from_bids_dataset
    get_sessions_from_bids_dataset
    """
    return [d for d in bids_path.glob("sub-*") if d.is_dir()]


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


# todo : place in genfi utils ?
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


# todo : use more ?
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


def validate_input_path(input_path: UserProvidedPath, check_exist: bool = True) -> Path:
    """Take a user provided path as input and returns a resolved, existing if check_exist is True, Path."""
    input_path = Path(input_path).resolve()
    if check_exist and not input_path.exists():
        raise FileNotFoundError(f"The provided input path {input_path} does not exist.")
    return input_path


def viscode_to_session(viscode: str, baseline_identifiers: Optional[set] = None) -> str:
    """Replace the session label from 'baseline_identifiers' with 'ses-M000'.

    If `viscode` is not a baseline identifier, parse the `viscode` and return
    a session identifier of the form `ses-MXXX`.

    The `viscode` is expected to be formatted in the following way:
    "m/M{session_number}", where session_number can be casted to an integer.
    Otherwise, a `ValueError` will be raised.

    Note: This function doesn't perform any check on the number of digits
    of the session identifier. It will return at least three digits encoded
    session numbers, but doesn't raise if more digits are needed (see
    examples section below).

    Parameters
    ----------
    viscode: str
        The name of the session.

    baseline_identifiers : set of str, optional
        Possible identifiers for baseline session.
        If the `viscode` is among these identifiers, `ses-M000` is returned.
        Default={"bl", "m0"}.

    Returns
    -------
    str:
        "ses-M000" if the session is the baseline session.
        Otherwise returns the original session name capitalized.

    Raises
    ------
    ValueError :
        If the `viscode` isn't formatted as expected.

    Examples
    --------
    >>> viscode_to_session("m1")
    'ses-M001'
    >>> viscode_to_session("M123")
    'ses-M123'
    >>> viscode_to_session("m1234")
    'ses-M1234'
    """
    baseline_identifiers = baseline_identifiers or {"bl", "m0"}
    if viscode in baseline_identifiers:
        return "ses-M000"
    session_pattern = "^[m,M][0-9]*"
    if re.match(session_pattern, viscode):
        return "ses-" + f"M{(int(viscode[1:])):03d}"
    raise ValueError(
        f"The viscode {viscode} is not correctly formatted."
        "Expected a session identifier of the form 'MXXX', "
        f"or a baseline identifier among {baseline_identifiers}."
    )


def load_clinical_csv(clinical_dir: Path, filename: str) -> pd.DataFrame:
    """Load the clinical csv from ADNI. This function is able to find the csv in the
    different known format available, the old format with just the name, and the new
    format with the name and the date of download.

    Parameters
    ----------
    clinical_dir: str
        Directory containing the csv.

    filename: str
        name of the file without the suffix.

    Returns
    -------
    pd.DataFrame:
        Dataframe corresponding to the filename.
    """
    pattern = filename + r"(_\d{1,2}[A-Za-z]{3}\d{4})?.csv"
    files_matching_pattern = [
        f for f in clinical_dir.rglob("[!.]*.csv") if re.search(pattern, f.name)
    ]
    if len(files_matching_pattern) != 1:
        raise IOError(
            f"Expecting to find exactly one file in folder {clinical_dir} "
            f"matching pattern {pattern}. {len(files_matching_pattern)} "
            f"files were found instead : \n{'- '.join(str(files_matching_pattern))}"
        )
    try:
        return pd.read_csv(files_matching_pattern[0], sep=",", low_memory=False)
    except Exception:
        raise ValueError(
            f"File {str(files_matching_pattern[0])} was found but could not "
            "be loaded as a DataFrame. Please check your data."
        )
