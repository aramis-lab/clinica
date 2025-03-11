"""This module contains utilities to grab or download files for Clinica."""

import hashlib
import os
from collections import namedtuple
from functools import partial
from pathlib import Path
from typing import Callable, Iterable, Optional, Sequence, Union

__all__ = [
    "RemoteFileStructure",
    "InvalidSubjectSession",
    "insensitive_glob",
    "find_images_path",
    "clinica_file_filter",
    "clinica_file_reader",
    "format_clinica_file_reader_errors",
    "clinica_list_of_files_reader",
    "clinica_group_reader",
    "fetch_file",
    "compute_sha256_hash",
    "get_file_from_server",
]

RemoteFileStructure = namedtuple("RemoteFileStructure", ["filename", "url", "checksum"])
InvalidSubjectSession = namedtuple("InvalidSubjectSession", ["subject", "session"])


def insensitive_glob(pattern_glob: str, recursive: Optional[bool] = False) -> list[str]:
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


def find_images_path(
    input_directory: os.PathLike,
    subject: str,
    session: str,
    errors: list[InvalidSubjectSession],
    valid_paths: list[str],
    is_bids: bool,
    pattern: str,
) -> None:
    """Appends the resulting path corresponding to subject, session and pattern in valid_paths.
    If an error is encountered, its (subject,session) couple is added to the list `errors`.

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

    errors : List
        List to which errors encountered in this function are added.

    valid_paths : List
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
            valid_paths.append(selected)
        else:
            errors.append(InvalidSubjectSession(subject, session))
    elif len(current_glob_found) == 0:
        errors.append(InvalidSubjectSession(subject, session))
    # Otherwise the file found is added to the result
    else:
        valid_paths.append(current_glob_found[0])


def _are_multiple_runs(files: list[str]) -> bool:
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


def _get_entities(files: list[Path], common_suffix: str) -> dict:
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
    files: list[Path],
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


def _select_run(files: list[str]) -> str:
    import numpy as np

    runs = [int(_get_run_number(f)) for f in files]
    return files[np.argmax(runs)]


def _get_run_number(filename: str) -> str:
    import re

    matches = re.match(r".*_run-(\d+).*", filename)
    if matches:
        return matches[1]
    raise ValueError(f"Filename {filename} should contain one and only one run entity.")


def _check_information(information: dict) -> None:
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

            if isinstance(item["pattern"], str) and item["pattern"][0] == "/":
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

        if isinstance(information["pattern"], str) and information["pattern"][0] == "/":
            raise ValueError(
                "pattern argument cannot start with char: / (does not work in os.path.join function). "
                "If you want to indicate the exact name of the file, use the format "
                "directory_name/filename.extension or filename.extension in the pattern argument."
            )


def clinica_file_filter(
    subjects: list[str],
    sessions: list[str],
    input_directory: Path,
    information: dict,
    n_procs: int = 1,
) -> tuple[list[str], list[str], list[str]]:
    from clinica.utils.stream import cprint

    files, errors = clinica_file_reader(
        subjects, sessions, input_directory, information, n_procs
    )
    cprint(format_clinica_file_reader_errors(errors, information), "warning")
    filtered_subjects, filtered_sessions = _remove_sub_ses_from_list(
        subjects, sessions, errors
    )
    return files, filtered_subjects, filtered_sessions


def format_clinica_file_reader_errors(
    errors: Sequence[InvalidSubjectSession], information: dict
) -> str:
    message = (
        f"Clinica encountered {len(errors)} "
        f"problem(s) while getting {information['description']}:\n"
    )
    if "needed_pipeline" in information and information["needed_pipeline"]:
        message += (
            "Please note that the following clinica pipeline(s) must "
            f"have run to obtain these files: {information['needed_pipeline']}\n"
        )
    if errors:
        message += "".join(f"\t* ({err.subject} | {err.session})\n" for err in errors)
        message += (
            "Clinica could not identify which file to use (missing or too many) for these sessions. "
            "They will not be processed."
        )
    return message


def _remove_sub_ses_from_list(
    subjects: list[str],
    sessions: list[str],
    errors: Iterable[InvalidSubjectSession],
) -> tuple[list[str], list[str]]:
    subjects = subjects.copy()
    sessions = sessions.copy()
    for invalid in errors:
        sub_indexes = [
            i for i, subject in enumerate(subjects) if subject == invalid.subject
        ]
        session_indexes = [
            i for i, session in enumerate(sessions) if session == invalid.session
        ]
        to_remove = list(set(sub_indexes) & set(session_indexes))
        to_remove.sort(reverse=True)
        for index in to_remove:
            subjects.pop(index)
            sessions.pop(index)
    return subjects, sessions


# todo : generalize
def clinica_file_reader(
    subjects: Sequence[str],
    sessions: Sequence[str],
    input_directory: os.PathLike,
    information: dict,
    n_procs: int = 1,
) -> tuple[list[str], list[InvalidSubjectSession]]:
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

    n_procs : int, optional
        Number of cores used to fetch files in parallel.
        If set to 1, subjects and sessions will be processed sequentially.
        Default=1.

    Returns
    -------
    results : List[str]
        List of files respecting the subject/session order provided in input.
              Iterable[InvalidSubjectSession]
        List of tuples (subject, session) which were identified as invalid (too many files or none).

    Raises
    ------
    TypeError
        If `information` is not a dictionary.

    ValueError
        If `information` is not formatted correctly. See function `_check_information`
        for more details.
        If the length of `subjects` is different from the length of `sessions`.

    Notes
    -----
    This function is case-insensitive, meaning that the pattern argument can, for example,
    contain upper case letters that do not exist in the existing file path.

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
    from clinica.dataset import DatasetType, check_dataset, get_dataset_type

    input_directory = Path(input_directory)
    _check_information(information)
    pattern = information["pattern"]
    check_dataset(input_directory)
    if len(subjects) != len(sessions):
        raise ValueError("Subjects and sessions must have the same length.")
    if len(subjects) == 0:
        return [], []
    file_reader = _read_files_parallel if n_procs > 1 else _read_files_sequential
    return file_reader(
        input_directory,
        subjects,
        sessions,
        get_dataset_type(input_directory) == DatasetType.RAW,
        pattern,
        n_procs=n_procs,
    )


def _read_files_parallel(
    input_directory: os.PathLike,
    subjects: Sequence[str],
    sessions: Sequence[str],
    is_bids: bool,
    pattern: str,
    n_procs: int,
) -> tuple[list[str], list[InvalidSubjectSession]]:
    from multiprocessing import Manager

    from joblib import Parallel, delayed

    manager = Manager()
    shared_results = manager.list()
    shared_errors_encountered = manager.list()
    Parallel(n_jobs=n_procs)(
        delayed(find_images_path)(
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
    subjects: Iterable[str],
    sessions: Iterable[str],
    is_bids: bool,
    pattern: str,
    **kwargs,
) -> tuple[list[str], list[InvalidSubjectSession]]:
    errors_encountered, results = [], []
    for sub, ses in zip(subjects, sessions):
        find_images_path(
            input_directory, sub, ses, errors_encountered, results, is_bids, pattern
        )
    return results, errors_encountered


def clinica_list_of_files_reader(
    participant_ids: list[str],
    session_ids: list[str],
    bids_or_caps_directory: os.PathLike,
    list_information: Iterable[dict],
    raise_exception: Optional[bool] = True,
) -> list[list[str]]:
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
    from .exceptions import ClinicaBIDSError

    all_errors = []
    list_found_files = []
    for info_file in list_information:
        files, errors = clinica_file_reader(
            participant_ids,
            session_ids,
            bids_or_caps_directory,
            info_file,
        )
        all_errors.append(errors)
        list_found_files.append([] if errors else files)

    if any(all_errors) and raise_exception:
        error_message = "Clinica faced error(s) while trying to read files in your BIDS or CAPS directory.\n"
        for error, info in zip(all_errors, list_information):
            error_message += format_clinica_file_reader_errors(error, info)
        raise ClinicaBIDSError(error_message)

    return list_found_files


def clinica_group_reader(
    caps_directory: os.PathLike,
    information: dict,
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
    from clinica.dataset import check_caps_dataset

    _check_information(information)
    pattern = information["pattern"]
    caps_directory = Path(caps_directory)
    check_caps_dataset(caps_directory)

    current_pattern = caps_directory / "**" / pattern
    found_files = insensitive_glob(str(current_pattern), recursive=True)

    # Since we are returning found_files[0], force raising even if raise_exception is False
    # Otherwise we'll get an uninformative IndexError...
    if (len(found_files) == 0) or (len(found_files) > 1 and raise_exception is True):
        _format_and_raise_group_reader_errors(caps_directory, found_files, information)

    return found_files[0]


def _format_and_raise_group_reader_errors(
    caps_directory: os.PathLike,
    found_files: list,
    information: dict,
) -> None:
    # todo : TEST
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


def compute_sha256_hash(file_path: Path) -> str:
    """Calculate the sha256 hash of the file at path."""

    sha256hash = hashlib.sha256()
    chunk_size = 8192
    with open(file_path, "rb") as f:
        while True:
            buffer = f.read(chunk_size)
            if not buffer:
                break
            sha256hash.update(buffer)
    return sha256hash.hexdigest()


def fetch_file(
    remote: RemoteFileStructure, output_folder: Union[str, os.PathLike]
) -> Path:
    """Download a specific file and save it into the resources folder of the package.

    Parameters
    ----------
    remote : RemoteFileStructure
        Structure containing url, filename and checksum.

    output_folder : str or PathLike
        Absolute path where the file will be downloaded.
        The name of the file will be the same as the downloaded file.

    Returns
    -------
    file_path : Path
        The path to the downloaded file.
    """
    import shutil
    import ssl
    from urllib.error import URLError
    from urllib.request import Request, urlopen

    from clinica.utils.stream import cprint

    output_folder = Path(output_folder)
    if not output_folder.exists():
        cprint(
            msg=f"The path {output_folder} to store the downloaded file does not exist",
            lvl="warning",
        )
        cprint(msg="Stop Clinica and handle this error", lvl="warning")
    file_path = output_folder / remote.filename
    # Download the file from `url` and save it locally under `file_name`:
    gcontext = ssl.SSLContext(protocol=ssl.PROTOCOL_TLS_CLIENT)
    gcontext.load_default_certs()
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

    if (checksum := compute_sha256_hash(file_path)) != remote.checksum:
        raise IOError(
            f"{file_path} has an SHA256 checksum ({checksum}) from expected "
            f"({remote.checksum}), file may be corrupted."
        )
    return file_path


def get_file_from_server(
    remote_file: RemoteFileStructure,
    cache_path: Optional[str] = None,
) -> Path:
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
    local_file : Path
        The path to the downloaded file.
    """
    from clinica.utils.stream import cprint

    if cache_path:
        cache_clinica = Path.home() / ".cache" / cache_path
    else:
        cache_clinica = Path.home() / ".cache" / "clinica" / "data"

    cache_clinica.mkdir(exist_ok=True, parents=True)
    local_file = cache_clinica / remote_file.filename

    if not local_file.exists():
        try:
            local_file = fetch_file(remote_file, cache_clinica)
        except IOError as err:
            cprint(
                msg=f"Unable to download {remote_file.filename} from {remote_file.url}: {err}",
                lvl="error",
            )

    return local_file
