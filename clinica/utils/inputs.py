"""This module contains utilities to grab or download files for Clinica."""

import hashlib
import os
from collections import namedtuple
from functools import partial
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

RemoteFileStructure = namedtuple("RemoteFileStructure", ["filename", "url", "checksum"])


def insensitive_glob(pattern_glob: str, recursive: Optional[bool] = False) -> List[str]:
    """This function is the glob.glob() function that is insensitive to the case.

    Parameters
    ----------
    pattern_glob : str
        Sensitive-to-the-case pattern.

    recursive : bool, optional
        Recursive parameter for `glob.glob()`.
        Default=False.

    Returns
    -------
    List[str] :
        Insensitive-to-the-case pattern.
    """
    from glob import glob

    def either(c: str) -> str:
        return "[%s%s]" % (c.lower(), c.upper()) if c.isalpha() else c

    return glob("".join(map(either, pattern_glob)), recursive=recursive)


def determine_caps_or_bids(input_dir: os.PathLike) -> bool:
    """Determine if the `input_dir` is a CAPS or a BIDS folder.

    Parameters
    ----------
    input_dir : os.PathLike
        The input folder.

    Returns
    -------
    bool :
        True if `input_dir` is a BIDS folder, False if `input_dir`
        is a CAPS folder or could not be determined.
    """
    input_dir = Path(input_dir)
    subjects_dir = input_dir / "subjects"
    groups_dir = input_dir / "groups"

    dir_to_look = subjects_dir if subjects_dir.is_dir() else input_dir
    subjects_sub_folders = _list_subjects_sub_folders(dir_to_look, groups_dir)
    if subjects_dir.is_dir():
        return False
    return len(subjects_sub_folders) > 0


def _list_subjects_sub_folders(
    root_dir: os.PathLike, groups_dir: os.PathLike
) -> List[os.PathLike]:
    from clinica.utils.stream import cprint

    root_dir = Path(root_dir)
    groups_dir = Path(groups_dir)
    warning_msg = (
        f"Could not determine if {groups_dir.parent} is a CAPS or BIDS directory. "
        "Clinica will assume this is a CAPS directory."
    )
    folder_content = [f for f in root_dir.iterdir()]
    subjects_sub_folders = [
        sub
        for sub in folder_content
        if (sub.name.startswith("sub-") and (root_dir / sub).is_dir())
    ]
    if len(subjects_sub_folders) == 0 and not groups_dir.is_dir():
        cprint(msg=warning_msg, lvl="warning")
    return subjects_sub_folders


def _common_checks(directory: os.PathLike, folder_type: str) -> None:
    """Utility function which performs checks common to BIDS and CAPS folder structures.

    Parameters
    ----------
    directory : PathLike
        Directory to check.

    folder_type : {"BIDS", "CAPS"}
        The type of directory.
    """
    from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError

    if not isinstance(directory, (os.PathLike, str)):
        raise ValueError(
            f"Argument you provided to check_{folder_type.lower()}_folder() is not a string."
        )

    error = ClinicaBIDSError if folder_type == "BIDS" else ClinicaCAPSError

    if not os.path.isdir(directory):
        raise error(
            f"The {folder_type} directory you gave is not a folder.\n"
            "Error explanations:\n"
            f"\t- Clinica expected the following path to be a folder: {directory}\n"
            "\t- If you gave relative path, did you run Clinica on the good folder?"
        )


def check_bids_folder(bids_directory: os.PathLike) -> None:
    """Check if provided `bids_directory` is a BIDS folder.

    Parameters
    ----------
    bids_directory : PathLike
        The input folder to check.

    Raises
    ------
    ValueError :
        If `bids_directory` is not a string.

    ClinicaBIDSError :
        If the provided path does not exist, or is not a directory.
        If the provided path is a CAPS folder (BIDS and CAPS could
        be swapped by user). We simply check that there is not a folder
        called 'subjects' in the provided path (that exists in CAPS hierarchy).
        If the provided folder is empty.
        If the provided folder does not contain at least one directory whose
        name starts with 'sub-'.
    """
    from clinica.utils.exceptions import ClinicaBIDSError

    bids_directory = Path(bids_directory)
    _common_checks(bids_directory, "BIDS")

    if (bids_directory / "subjects").is_dir():
        raise ClinicaBIDSError(
            f"The BIDS directory ({bids_directory}) you provided seems to "
            "be a CAPS directory due to the presence of a 'subjects' folder."
        )

    if len([f for f in bids_directory.iterdir()]) == 0:
        raise ClinicaBIDSError(
            f"The BIDS directory you provided is empty. ({bids_directory})."
        )

    subj = [f for f in bids_directory.iterdir() if f.name.startswith("sub-")]
    if len(subj) == 0:
        raise ClinicaBIDSError(
            "Your BIDS directory does not contains a single folder whose name "
            "starts with 'sub-'. Check that your folder follow BIDS standard."
        )


def check_caps_folder(caps_directory: os.PathLike) -> None:
    """Check if provided `caps_directory`is a CAPS folder.

    Parameters
    ----------
    caps_directory : os.PathLike
        The input folder to check.

    Raises
    ------
    ValueError :
        If `caps_directory` is not a string.

    ClinicaCAPSError :
        If the provided path does not exist, or is not a directory.
        If the provided path is a BIDS folder (BIDS and CAPS could be
        swapped by user). We simply check that there is not a folder
        whose name starts with 'sub-' in the provided path (that exists
        in BIDS hierarchy).

    Notes
    -----
    Keep in mind that a CAPS folder can be empty.
    """
    from clinica.utils.exceptions import ClinicaCAPSError

    caps_directory = Path(caps_directory)
    _common_checks(caps_directory, "CAPS")

    sub_folders = [f for f in caps_directory.iterdir() if f.name.startswith("sub-")]
    if len(sub_folders) > 0:
        error_string = (
            "Your CAPS directory contains at least one folder whose name "
            "starts with 'sub-'. Check that you did not swap BIDS and CAPS folders.\n"
            "Folder(s) found that match(es) BIDS architecture:\n"
        )
        for directory in sub_folders:
            error_string += f"\t{directory}\n"
        error_string += (
            "A CAPS directory has a folder 'subjects' at its root, in which "
            "are stored the output of the pipeline for each subject."
        )
        raise ClinicaCAPSError(error_string)


def find_sub_ses_pattern_path(
    input_directory: os.PathLike,
    subject: str,
    session: str,
    error_encountered: list,
    results: list,
    is_bids: bool,
    pattern: str,
) -> None:
    """Appends the output path corresponding to subject, session and pattern in results.

    If an error is encountered, its corresponding message is added to the list `error_encountered`.

    Parameters
    ----------
    input_directory : str
        Path to the root of the input directory (BIDS or CAPS).

        .. warning::
            This function does not perform any check on `input_directory`.
            It is assumed that it has been previously checked by either
            `check_bids_directory` or `check_caps_directory`, and that
            the flag `is_bids` has been set accordingly.

    subject : str
        Name given to the folder of a participant (ex: sub-ADNI002S0295).

    session : str
        Name given to the folder of a session (ex: ses-M00).

    error_encountered : List
        List to which errors encountered in this function are added.

    results : List
        List to which the output path corresponding to subject, session
        and pattern is added.

    is_bids : bool
        True if `input_dir` is a BIDS folder, False if `input_dir` is a
        CAPS folder.

    pattern : str
        Define the pattern of the final file.
    """
    from clinica.utils.stream import cprint

    input_directory = Path(input_directory)
    if is_bids:
        origin_pattern = input_directory / subject / session
    else:
        origin_pattern = input_directory / "subjects" / subject / session

    current_pattern = origin_pattern / "**" / pattern
    current_glob_found = insensitive_glob(str(current_pattern), recursive=True)
    if len(current_glob_found) > 1:
        # If we have more than one file at this point, there are two possibilities:
        #   - there is a problem somewhere which made us catch too many files
        #           --> In this case, we raise an error.
        #   - we have captured multiple runs for the same subject and session
        #           --> In this case, we need to select one of these runs to proceed.
        #               Ideally, this should be done via QC but for now, we simply
        #               select the latest run and warn the user about it.
        if _are_multiple_runs(current_glob_found):
            selected = _select_run(current_glob_found)
            list_of_found_files_for_reporting = ""
            for filename in current_glob_found:
                list_of_found_files_for_reporting += f"- {filename}\n"
            cprint(
                f"More than one run were found for subject {subject} and session {session} : "
                f"\n\n{list_of_found_files_for_reporting}\n"
                f"Clinica will proceed with the latest run available, that is \n\n-{selected}.",
                lvl="warning",
            )
            results.append(selected)
        else:
            error_str = f"\t*  ({subject} | {session}): More than 1 file found:\n"
            for found_file in current_glob_found:
                error_str += f"\t\t{found_file}\n"
            error_encountered.append(error_str)
    elif len(current_glob_found) == 0:
        error_encountered.append(f"\t* ({subject} | {session}): No file found\n")
    # Otherwise the file found is added to the result
    else:
        results.append(current_glob_found[0])


def _are_multiple_runs(files: List[str]) -> bool:
    """Returns whether the files in the provided list only differ through their run number.

    The provided files must have exactly the same parent paths, extensions, and BIDS entities
    excepted for the 'run' entity which must be different.

    Parameters
    ----------
    files : List of str
        The files to analyze.

    Returns
    -------
    bool :
        True if the provided files only differ through their run number, False otherwise.
    """
    from pathlib import Path

    files = [Path(_) for _ in files]
    # Exit quickly if less than one file or if at least one file does not have the entity run
    if len(files) < 2 or any(["_run-" not in f.name for f in files]):
        return False
    try:
        _check_common_parent_path(files)
        _check_common_extension(files)
        common_suffix = _check_common_suffix(files)
    except ValueError:
        return False
    found_entities = _get_entities(files, common_suffix)
    for entity_name, entity_values in found_entities.items():
        if entity_name != "run":
            # All entities except run numbers should be the same
            if len(entity_values) != 1:
                return False
        else:
            # Run numbers should differ otherwise this is a BIDS violation at this point
            if len(entity_values) != len(files):
                return False
    return True


def _get_entities(files: List[Path], common_suffix: str) -> dict:
    """Compute a dictionary where the keys are entity names and the values
    are sets of all the corresponding entity values found while iterating over
    the provided files.

    Parameters
    ----------
    files : List of Path
        List of paths to get entities of.

    common_suffix : str
        The suffix common to all the files. This suffix will be stripped
        from the file names in order to only analyze BIDS entities.

    Returns
    -------
    dict :
        The entities dictionary.
    """
    from clinica.utils.filemanip import get_filename_no_ext

    found_entities = dict()
    for f in files:
        entities = get_filename_no_ext(f.name).rstrip(common_suffix).split("_")
        for entity in entities:
            entity_name, entity_value = entity.split("-")
            if entity_name in found_entities:
                found_entities[entity_name].add(entity_value)
            else:
                found_entities[entity_name] = {entity_value}

    return found_entities


def _check_common_properties_of_files(
    files: List[Path],
    property_name: str,
    property_extractor: Callable,
) -> str:
    """Verify that all provided files share the same property and return its value.

    Parameters
    ----------
    files : List of Paths
        List of file paths for which to verify common property.

    property_name : str
        The name of the property to verify.

    property_extractor : Callable
        The function which is responsible for the property extraction.
        It must implement the interface `property_extractor(filename: Path) -> str`

    Returns
    -------
    str :
        The value of the common property.

    Raises
    ------
    ValueError :
        If the provided files do not have the same given property.
    """
    extracted_properties = {property_extractor(f) for f in files}
    if len(extracted_properties) != 1:
        raise ValueError(
            f"The provided files do not share the same {property_name}."
            f"The following {property_name}s were found: {extracted_properties}"
        )
    return extracted_properties.pop()


def _get_parent_path(filename: Path) -> str:
    return str(filename.parent)


def _get_extension(filename: Path) -> str:
    return "".join(filename.suffixes)


def _get_suffix(filename: Path) -> str:
    from clinica.utils.filemanip import get_filename_no_ext

    return f"_{get_filename_no_ext(filename.name).split('_')[-1]}"


_check_common_parent_path = partial(
    _check_common_properties_of_files,
    property_name="parent path",
    property_extractor=_get_parent_path,
)


_check_common_extension = partial(
    _check_common_properties_of_files,
    property_name="extension",
    property_extractor=_get_extension,
)


_check_common_suffix = partial(
    _check_common_properties_of_files,
    property_name="suffix",
    property_extractor=_get_suffix,
)


def _select_run(files: List[str]) -> str:
    import numpy as np

    runs = [int(_get_run_number(f)) for f in files]
    return files[np.argmax(runs)]


def _get_run_number(filename: str) -> str:
    import re

    matches = re.match(r".*_run-(\d+).*", filename)
    if matches:
        return matches[1]
    raise ValueError(f"Filename {filename} should contain one and only one run entity.")


def _check_information(information: Dict) -> None:
    if not isinstance(information, (dict, list)):
        raise TypeError(
            "A dict or list of dicts must be provided for the argument 'information'"
        )

    if isinstance(information, list):
        for item in information:
            if not all(elem in item for elem in ["pattern", "description"]):
                raise ValueError(
                    "'information' must contain the keys 'pattern' and 'description'"
                )

            if not all(
                elem in ["pattern", "description", "needed_pipeline"]
                for elem in item.keys()
            ):
                raise ValueError(
                    "'information' can only contain the keys 'pattern', 'description' and 'needed_pipeline'"
                )

            if item["pattern"][0] == "/":
                raise ValueError(
                    "pattern argument cannot start with char: / (does not work in os.path.join function). "
                    "If you want to indicate the exact name of the file, use the format "
                    "directory_name/filename.extension or filename.extension in the pattern argument."
                )
    else:
        if not all(elem in information for elem in ["pattern", "description"]):
            raise ValueError(
                "'information' must contain the keys 'pattern' and 'description'"
            )

        if not all(
            elem in ["pattern", "description", "needed_pipeline"]
            for elem in information.keys()
        ):
            raise ValueError(
                "'information' can only contain the keys 'pattern', 'description' and 'needed_pipeline'"
            )

        if information["pattern"][0] == "/":
            raise ValueError(
                "pattern argument cannot start with char: / (does not work in os.path.join function). "
                "If you want to indicate the exact name of the file, use the format "
                "directory_name/filename.extension or filename.extension in the pattern argument."
            )


def _format_errors(errors: List, information: Dict) -> str:
    error_message = (
        f"Clinica encountered {len(errors)} "
        f"problem(s) while getting {information['description']}:\n"
    )
    if "needed_pipeline" in information and information["needed_pipeline"]:
        error_message += (
            "Please note that the following clinica pipeline(s) must "
            f"have run to obtain these files: {information['needed_pipeline']}\n"
        )
    error_message += "\n".join(errors)

    return error_message


def clinica_file_reader(
    subjects: List[str],
    sessions: List[str],
    input_directory: os.PathLike,
    information: Dict,
    raise_exception: Optional[bool] = True,
    n_procs: Optional[int] = 1,
):
    """Read files in BIDS or CAPS directory based on participant ID(s).

    This function grabs files relative to a subject and session list according to a glob pattern (using *)

    Parameters
    ----------
    subjects : List[str]
        List of subjects.

    sessions : List[str]
        List of sessions. Must be same size as `subjects` and must correspond.

    input_directory : PathLike
        Path to the BIDS or CAPS directory to read from.

    information : Dict
        Dictionary containing all the relevant information to look for the files.
        The possible keys are:

            - `pattern`: Required. Define the pattern of the final file.
            - `description`: Required. String to describe what the file is.
            - `needed_pipeline` : Optional. String describing the pipeline(s)
              needed to obtain the related file.

    raise_exception : bool, optional
        If True, an exception is raised if errors happen. If not, we return the file
        list as it is. Default=True.

    n_procs : int, optional
        Number of cores used to fetch files in parallel.
        If set to 1, subjects and sessions will be processed sequentially.
        Default=1.

    Returns
    -------
    results : List[str]
        List of files respecting the subject/session order provided in input.

    error_message : str
        Error message which contains all errors encountered while reading the files.

    Raises
    ------
    TypeError
        If `information` is not a dictionary.

    ValueError
        If `information` is not formatted correctly. See function `_check_information`
        for more details.
        If the length of `subjects` is different from the length of `sessions`.

    ClinicaCAPSError or ClinicaBIDSError
        If multiples files are found for 1 subject/session, or if no file is found.

        .. note::
            If `raise_exception` is False, no exception is raised.

    Notes
    -----
    This function is case-insensitive, meaning that the pattern argument can, for example,
    contain upper case letters that do not exist in the existing file path.

    You should always use `clinica_file_reader` in the following manner:

    .. code-block:: python

         try:
            file_list = clinica_file_reader(...)
         except ClinicaException as e:
            # Deal with the error

    Examples
    --------
    The paths are shortened for readability.

    You have the full name of a file.

    File `orig_nu.mgz` from FreeSurfer of subject `sub-ADNI011S4105`, session `ses-M00`
    located in mri folder of FreeSurfer output :

    >>> clinica_file_reader(
            ['sub-ADNI011S4105'],
            ['ses-M00'],
            caps_directory,
            {
                'pattern': 'freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz',
                'description': 'freesurfer file orig_nu.mgz',
                'needed_pipeline': 't1-freesurfer'
            }
        )
    ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/mri/orig_nu.mgz']

    You have a partial name of the file.

    File `sub-ADNI011S4105_ses-M00_trc-18FFDG_pet.nii.gz` in BIDS directory.
    Here, filename depends on subject and session name :

    >>> clinica_file_reader(
            ['sub-ADNI011S4105'],
            ['ses-M00'],
            bids_directory,
            {
                'pattern': '*18FFDG_pet.nii*',
                'description': 'FDG PET data'
            }
        )
    ['/bids/sub-ADNI011S4105/ses-M00/pet/sub-ADNI011S4105_ses-M00_trc-18FFDG_pet.nii.gz']

    Tricky example.

    Get the file `rh.white` from FreeSurfer :

    This will fail :

    >>> clinica_file_reader(
            ['sub-ADNI011S4105'],
            ['ses-M00'],
            caps,
            {
                'pattern': 'rh.white',
                'description': 'right hemisphere of outer cortical surface.',
                'needed_pipeline': 't1-freesurfer'
            }
        )
    * More than 1 file found::
            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/fsaverage/surf/rh.white
            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/rh.EC_average/surf/rh.white
            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white

    Correct usage (e.g. in pet-surface): pattern string must be 'sub-*_ses-*/surf/rh.white',
    or even more precise: 't1/freesurfer_cross_sectional/sub-*_ses-*/surf/rh.white'
    It then gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white']
    """
    from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError

    _check_information(information)
    pattern = information["pattern"]

    input_directory = Path(input_directory)
    is_bids = determine_caps_or_bids(input_directory)
    if is_bids:
        check_bids_folder(input_directory)
    else:
        check_caps_folder(input_directory)

    if len(subjects) != len(sessions):
        raise ValueError("Subjects and sessions must have the same length.")

    if len(subjects) == 0:
        return [], ""

    file_reader = _read_files_parallel if n_procs > 1 else _read_files_sequential
    results, errors_encountered = file_reader(
        input_directory,
        subjects,
        sessions,
        is_bids,
        pattern,
        n_procs=n_procs,
    )
    error_message = _format_errors(errors_encountered, information)

    if len(errors_encountered) > 0 and raise_exception:
        if is_bids:
            raise ClinicaBIDSError(error_message)
        else:
            raise ClinicaCAPSError(error_message)

    return results, error_message


def _read_files_parallel(
    input_directory: os.PathLike,
    subjects: List[str],
    sessions: List[str],
    is_bids: bool,
    pattern: str,
    n_procs: int,
) -> Tuple[List[str], List[str]]:
    from multiprocessing import Manager

    from joblib import Parallel, delayed

    manager = Manager()
    shared_results = manager.list()
    shared_errors_encountered = manager.list()
    Parallel(n_jobs=n_procs)(
        delayed(find_sub_ses_pattern_path)(
            input_directory,
            sub,
            ses,
            shared_errors_encountered,
            shared_results,
            is_bids,
            pattern,
        )
        for sub, ses in zip(subjects, sessions)
    )
    results = list(shared_results)
    errors_encountered = list(shared_errors_encountered)
    return results, errors_encountered


def _read_files_sequential(
    input_directory: os.PathLike,
    subjects: List[str],
    sessions: List[str],
    is_bids: bool,
    pattern: str,
    **kwargs,
) -> Tuple[List[str], List[str]]:
    errors_encountered, results = [], []
    for sub, ses in zip(subjects, sessions):
        find_sub_ses_pattern_path(
            input_directory, sub, ses, errors_encountered, results, is_bids, pattern
        )
    return results, errors_encountered


def clinica_list_of_files_reader(
    participant_ids: List[str],
    session_ids: List[str],
    bids_or_caps_directory: os.PathLike,
    list_information: List[Dict],
    raise_exception: Optional[bool] = True,
) -> List[List[str]]:
    """Read list of BIDS or CAPS files.

    This function iterates calls of `clinica_file_reader` to extract input files based
    on information given by `list_information`.

    Parameters
    ----------
    participant_ids : List[str]
        List of participant IDs.
        Example: ['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02']

    session_ids : List[str]
        List of sessions ID associated to `participant_ids`
        Example: ['ses-M00', 'ses-M18', 'ses-M00']

    bids_or_caps_directory : PathLike
        Path to the BIDS of CAPS directory to read from.

    list_information : List[Dict]
        List of information dictionaries described in `clinica_file_reader`.

    raise_exception : bool, optional
        Raise Exception or not. Defaults to True.

    Returns
    -------
    list_found_files : List[List[str]]
        List of lists of found files following order of `list_information`
    """
    from .exceptions import ClinicaBIDSError, ClinicaException

    all_errors = []
    list_found_files = []
    for info_file in list_information:
        try:
            list_found_files.append(
                clinica_file_reader(
                    participant_ids,
                    session_ids,
                    bids_or_caps_directory,
                    info_file,
                    True,
                )[0]
            )
        except ClinicaException as e:
            list_found_files.append([])
            all_errors.append(e)

    if len(all_errors) > 0 and raise_exception:
        error_message = "Clinica faced error(s) while trying to read files in your BIDS or CAPS directory.\n"
        for msg in all_errors:
            error_message += str(msg)
        raise ClinicaBIDSError(error_message)

    return list_found_files


def clinica_group_reader(
    caps_directory: os.PathLike,
    information: Dict,
    raise_exception: Optional[bool] = True,
) -> str:
    """Read files from CAPS directory based on group ID(s).

    This function grabs files relative to a group, according to a glob pattern (using *).
    Only one file can be returned, as order is arbitrary in `glob.glob()`.

    Parameters
    ----------
    caps_directory : PathLike
        Path to the input CAPS directory.

    information : Dict
        Dictionary containing all the relevant information to look for the files.
        The possible keys are:

            - `pattern`: Required. Define the pattern of the final file.
            - `description`: Required. String to describe what the file is.
            - `needed_pipeline` : Optional. String describing the pipeline(s)
              needed to obtain the related file.

    raise_exception : bool, optional
        If True, an exception is raised if errors happen.
        If not, we return the file list as it is.
        Default=True.

    Returns
    -------
    str :
        Path to the found file.

    Raises
    ------
    ClinicaCAPSError :
        If no file is found, or more than 1 files are found.
    """
    _check_information(information)
    pattern = information["pattern"]
    caps_directory = Path(caps_directory)
    check_caps_folder(caps_directory)

    current_pattern = caps_directory / "**" / pattern
    found_files = insensitive_glob(str(current_pattern), recursive=True)

    # Since we are returning found_files[0], force raising even if raise_exception is False
    # Otherwise we'll get an uninformative IndexError...
    if (len(found_files) == 0) or (len(found_files) > 1 and raise_exception is True):
        _format_and_raise_group_reader_errors(caps_directory, found_files, information)

    return found_files[0]


def _format_and_raise_group_reader_errors(
    caps_directory: os.PathLike,
    found_files: List,
    information: Dict,
) -> None:
    from clinica.utils.exceptions import ClinicaCAPSError

    error_string = (
        f"Clinica encountered a problem while getting {information['description']}. "
    )
    if len(found_files) == 0:
        error_string += "No file was found"
    else:
        error_string += f"{len(found_files)} files were found:"
        for found_file in found_files:
            error_string += f"\n\t{found_file}"
        error_string += (
            f"\n\tCAPS directory: {caps_directory}\n"
            "Please note that the following clinica pipeline(s) must have run to obtain these files: "
            f"{information['needed_pipeline']}\n"
        )
    raise ClinicaCAPSError(error_string)


def _sha256(path):
    """Calculate the sha256 hash of the file at path."""
    sha256hash = hashlib.sha256()
    chunk_size = 8192
    with open(path, "rb") as f:
        while True:
            buffer = f.read(chunk_size)
            if not buffer:
                break
            sha256hash.update(buffer)
    return sha256hash.hexdigest()


def fetch_file(remote: RemoteFileStructure, dirname: Optional[str]) -> str:
    """Download a specific file and save it into the resources folder of the package.

    Parameters
    ----------
    remote : RemoteFileStructure
        Structure containing url, filename and checksum.

    dirname : str
        Absolute path where the file will be downloaded.

    Returns
    -------
    file_path : str
        Absolute file path.
    """
    import os.path
    import shutil
    import ssl
    from urllib.error import URLError
    from urllib.request import Request, urlopen

    from clinica.utils.stream import cprint

    if not os.path.exists(dirname):
        cprint(msg="Path to the file does not exist", lvl="warning")
        cprint(msg="Stop Clinica and handle this error", lvl="warning")

    file_path = os.path.join(dirname, remote.filename)
    # Download the file from `url` and save it locally under `file_name`:
    gcontext = ssl.SSLContext()
    req = Request(remote.url + remote.filename)
    try:
        response = urlopen(req, context=gcontext)
    except URLError as e:
        if hasattr(e, "reason"):
            cprint(msg=f"We failed to reach a server. Reason: {e.reason}", lvl="error")
        elif hasattr(e, "code"):
            cprint(
                msg=f"The server could not fulfill the request. Error code: {e.code}",
                lvl="error",
            )
    else:
        try:
            with open(file_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)
        except OSError as err:
            cprint(msg="OS error: {0}".format(err), lvl="error")

    checksum = _sha256(file_path)
    if remote.checksum != checksum:
        raise IOError(
            f"{file_path} has an SHA256 checksum ({checksum}) from expected "
            f"({remote.checksum}), file may be corrupted."
        )
    return file_path


def get_file_from_server(
    remote_file: RemoteFileStructure,
    cache_path: Optional[str] = None,
) -> str:
    """Download file from server.

    Parameters
    ----------
    remote_file : RemoteFileStructure
        Remote file structure defined in `clinica.utils.inputs` specifying
        the file to retrieve.

    cache_path : str, optional
        Path to cache. Default="~/.cache/clinica/data".

    Returns
    -------
    local_file : str
        Path to the downloaded file.
    """
    import os
    from pathlib import Path

    from clinica.utils.stream import cprint

    home = str(Path.home())
    if cache_path:
        cache_clinica = os.path.join(home, ".cache", cache_path)
    else:
        cache_clinica = os.path.join(home, ".cache", "clinica", "data")

    os.makedirs(cache_clinica, exist_ok=True)

    local_file = os.path.join(cache_clinica, remote_file.filename)

    if not (os.path.exists(local_file)):
        try:
            local_file = fetch_file(remote_file, cache_clinica)
        except IOError as err:
            cprint(
                msg=f"Unable to download {remote_file.filename} from {remote_file.url}: {err}",
                lvl="error",
            )

    return local_file
