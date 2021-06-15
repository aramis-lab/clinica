# coding: utf8

"""This module contains utilities to grab or download files for Clinica."""

import hashlib
from collections import namedtuple

RemoteFileStructure = namedtuple("RemoteFileStructure", ["filename", "url", "checksum"])


def insensitive_glob(pattern_glob, recursive=False):
    """This function is the glob.glob() function that is insensitive to the case.

    Args:
        pattern_glob: sensitive-to-the-case pattern
        recursive: recursive parameter for glob.glob()

    Returns:
         insensitive-to-the-case pattern
    """
    from glob import glob

    def either(c):
        return "[%s%s]" % (c.lower(), c.upper()) if c.isalpha() else c

    return glob("".join(map(either, pattern_glob)), recursive=recursive)


def determine_caps_or_bids(input_dir):
    """Determine if the input is a CAPS or a BIDS folder.

    Args
        input_dir: input folder

    Returns:
        True if input_dir is a bids, False if input_dir is a CAPS

    Raise:
        RuntimeError if function could not determine if BIDS or CAPS or whatever else
    """
    from os import listdir
    from os.path import isdir, join
    from colorama import Fore

    from clinica.utils.stream import cprint

    if isdir(join(input_dir, "subjects")):
        if len(
            [
                sub
                for sub in listdir(join(input_dir, "subjects"))
                if (sub.startswith("sub-") and isdir(join(input_dir, "subjects", sub)))
            ]
        ) > 0 or isdir(join(input_dir, "groups")):
            return False
        else:
            cprint(f"{Fore.YELLOW}[Warning] Could not determine if {input_dir} is a CAPS or BIDS directory. Clinica will asume this is a CAPS directory.{Fore.RESET}")
            return False

    else:
        if (
            len(
                [
                    sub
                    for sub in listdir(input_dir)
                    if (sub.startswith("sub-") and isdir(join(input_dir, sub)))
                ]
            )
            > 0
        ):
            return True
        else:
            if isdir(join(input_dir, "groups")):
                return False
            else:
                cprint(f"{Fore.YELLOW}[Warning] Could not determine if {input_dir} is a CAPS or BIDS directory. Clinica will asume this is a CAPS directory.{Fore.RESET}")
                return False


def check_bids_folder(bids_directory):
    """Check BIDS folder.

    This function checks the following items:
        - bids_directory is a string
        - the provided path exists and is a directory
        - provided path is not a CAPS folder (BIDS and CAPS could be swapped by user). We simply check that there is
          not a folder called 'subjects' in the provided path (that exists in CAPS hierarchy)
        - provided folder is not empty
        - provided folder must contains at least one directory whose name starts with 'sub-'
    """
    from os import listdir
    from os.path import isdir, join

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaBIDSError

    assert isinstance(
        bids_directory, str
    ), "Argument you provided to check_bids_folder() is not a string."

    if not isdir(bids_directory):
        raise ClinicaBIDSError(
            f"{Fore.RED}\n[Error] The BIDS directory you gave is not a folder.\n{Fore.RESET}"
            f"{Fore.YELLOW}\nError explanations:\n{Fore.RESET}"
            f" - Clinica expected the following path to be a folder: {Fore.BLUE}{bids_directory}{Fore.RESET}\n"
            f" - If you gave relative path, did you run Clinica on the good folder?"
        )

    if isdir(join(bids_directory, "subjects")):
        raise ClinicaBIDSError(
            f"{Fore.RED}\n[Error] The BIDS directory ({bids_directory}) you provided seems to "
            f"be a CAPS directory due to the presence of a 'subjects' folder.{Fore.RESET}"
        )

    if len(listdir(bids_directory)) == 0:
        raise ClinicaBIDSError(
            f"{Fore.RED}\n[Error] The BIDS directory you provided  is empty. ({bids_directory}).{Fore.RESET}"
        )

    if len([item for item in listdir(bids_directory) if item.startswith("sub-")]) == 0:
        raise ClinicaBIDSError(
            f"{Fore.RED}\n[Error] Your BIDS directory does not contains a single folder whose name "
            f"starts with 'sub-'. Check that your folder follow BIDS standard.{Fore.RESET}"
        )


def check_caps_folder(caps_directory):
    """Check CAPS folder.

    This function checks the following items:
        - caps_directory is a string
        - the provided path exists and is a directory
        - provided path is not a BIDS folder (BIDS and CAPS could be swapped by user). We simply check that there is
          not a folder whose name starts with 'sub-' in the provided path (that exists in BIDS hierarchy)
    Keep in mind that CAPS folder can be empty.
    """
    import os
    from os import listdir

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaCAPSError

    assert isinstance(
        caps_directory, str
    ), "Argument you provided to check_caps_folder() is not a string."

    if not os.path.isdir(caps_directory):
        raise ClinicaCAPSError(
            f"{Fore.RED}\n[Error] The CAPS directory you gave is not a folder.\n{Fore.RESET}"
            f"{Fore.YELLOW}\nError explanations:\n{Fore.RESET}"
            f" - Clinica expected the following path to be a folder: {Fore.BLUE}{caps_directory}{Fore.RESET}\n"
            f" - If you gave relative path, did you run Clinica on the good folder?"
        )

    sub_folders = [item for item in listdir(caps_directory) if item.startswith("sub-")]
    if len(sub_folders) > 0:
        error_string = (
            "\n[Error] Your CAPS directory contains at least one folder whose name "
            "starts with 'sub-'. Check that you did not swap BIDS and CAPS folders.\n"
            " Folder(s) found that match(es) BIDS architecture:\n"
        )
        for directory in sub_folders:
            error_string += f"\t{directory}\n"
        error_string += (
            "A CAPS directory has a folder 'subjects' at its root, in which "
            "are stored the output of the pipeline for each subject."
        )
        raise ClinicaCAPSError(error_string)


def clinica_file_reader(
    subjects, sessions, input_directory, information, raise_exception=True
):
    """Read files in BIDS or CAPS directory based on participant ID(s).

    This function grabs files relative to a subject and session list according to a glob pattern (using *)
    Args:
        subjects: list of subjects
        sessions: list of sessions (must be same size as subjects, and must correspond )
        input_directory: location of the bids or caps directory
        information: dictionary containing all the relevant information to look for the files. Dict must contains the
                     following keys : pattern, description. The optional key is: needed_pipeline
                             pattern: define the pattern of the final file
                             description: string to describe what the file is
                             needed_pipeline (optional): string describing the pipeline(s) needed to obtain the related
                                                        file
        raise_exception: if True (normal behavior), an exception is raised if errors happen. If not, we return the file
                        list as it is

    Returns:
         list of files respecting the subject/session order provided in input,
         You should always use clinica_file_reader in the following manner:
         try:
            file_list = clinica_file_reader(...)
         except ClinicaException as e:
            # Deal with the error

    Raises:
        ClinicaCAPSError or ClinicaBIDSError if multiples files are found for 1 subject/session, or no file is found
        If raise_exception is False, no exception is raised

        Examples: (path are shortened for readability)
            - You have the full name of a file:
                File orig_nu.mgz from FreeSurfer of subject sub-ADNI011S4105 session ses-M00 located in mri folder of
                FreeSurfer output :
                    clinica_file_reader(['sub-ADNI011S4105'],
                                        ['ses-M00'],
                                        caps_directory,
                                        {'pattern': 'freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz',
                                         'description': 'freesurfer file orig_nu.mgz',
                                         'needed_pipeline': 't1-freesurfer'})
                    gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/mri/orig_nu.mgz']

            - You have a partial name of the file:
                File sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz in BIDS directory. Here, filename depends on
                subject and session name :
                     clinica_file_reader(['sub-ADNI011S4105'],
                                         ['ses-M00'],
                                         bids_directory,
                                         {'pattern': '*fdg_pet.nii*',
                                          'description': 'FDG PET data'})
                     gives: ['/bids/sub-ADNI011S4105/ses-M00/pet/sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz']

            - Tricky example:
                Get the file rh.white from FreeSurfer:
                If you try:
                    clinica_file_reader(['sub-ADNI011S4105'],
                                        ['ses-M00'],
                                        caps,
                                        {'pattern': 'rh.white',
                                         'description': 'right hemisphere of outter cortical surface.',
                                         'needed_pipeline': 't1-freesurfer'})
                        the following error will arise:
                        * More than 1 file found::
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/fsaverage/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/rh.EC_average/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white
                Correct usage (e.g. in pet-surface): pattern string must be 'sub-*_ses-*/surf/rh.white' or even more precise:
                        't1/freesurfer_cross_sectional/sub-*_ses-*/surf/rh.white'
                    It then gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white']

        Note:
            This function is case insensitive, meaning that the pattern argument can, for example, contain maj letter
            that do not exists in the existing file path.

    """
    from os.path import join

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError

    assert isinstance(
        information, dict
    ), "A dict must be provided for the argument 'dict'"
    assert all(
        elem in information.keys() for elem in ["pattern", "description"]
    ), "'information' must contain the keys 'pattern' and 'description'"
    assert all(
        elem in ["pattern", "description", "needed_pipeline"]
        for elem in information.keys()
    ), "'information' can only contain the keys 'pattern', 'description' and 'needed_pipeline'"

    pattern = information["pattern"]
    is_bids = determine_caps_or_bids(input_directory)

    if is_bids:
        check_bids_folder(input_directory)
    else:
        check_caps_folder(input_directory)

    # Some check on the formatting on the data
    assert pattern[0] != "/", (
        "pattern argument cannot start with char: / (does not work in os.path.join function). "
        "If you want to indicate the exact name of the file, use the format"
        " directory_name/filename.extension or filename.extension in the pattern argument"
    )
    assert len(subjects) == len(
        sessions
    ), "Subjects and sessions must have the same length"
    if len(subjects) == 0:
        return []

    # results is the list containing the results
    results = []
    # error is the list of the errors that happen during the whole process
    error_encountered = []
    for sub, ses in zip(subjects, sessions):
        if is_bids:
            origin_pattern = join(input_directory, sub, ses)
        else:
            origin_pattern = join(input_directory, "subjects", sub, ses)

        current_pattern = join(origin_pattern, "**/", pattern)
        current_glob_found = insensitive_glob(current_pattern, recursive=True)

        # Error handling if more than 1 file are found, or when no file is found
        if len(current_glob_found) > 1:
            error_str = f"\t* {Fore.BLUE}  ({sub} | {ses}) {Fore.RESET}: More than 1 file found:\n"
            for found_file in current_glob_found:
                error_str += f"\t\t{found_file}\n"
            error_encountered.append(error_str)
        elif len(current_glob_found) == 0:
            error_encountered.append(
                f"\t* {Fore.BLUE} ({sub} | {ses}) {Fore.RESET}: No file found\n"
            )
        # Otherwise the file found is added to the result
        else:
            results.append(current_glob_found[0])

    # We do not raise an error, so that the developper can gather all the problems before Clinica crashes
    if len(error_encountered) > 0 and raise_exception is True:
        error_message = (
            f"{Fore.RED}\n[Error] Clinica encountered {len(error_encountered)} "
            f"problem(s) while getting {information['description']}:\n{Fore.RESET}"
        )
        if "needed_pipeline" in information.keys():
            if information["needed_pipeline"]:
                error_message += (
                    f"{Fore.YELLOW}Please note that the following clinica pipeline(s) must "
                    f"have run to obtain these files: {information['needed_pipeline']}{Fore.RESET}\n"
                )
        for msg in error_encountered:
            error_message += msg
        if is_bids:
            raise ClinicaBIDSError(error_message)
        else:
            raise ClinicaCAPSError(error_message)
    return results


def clinica_group_reader(caps_directory, information, raise_exception=True):
    """Read files CAPS directory based on group ID(s).

    This function grabs files relative to a group, according to a glob pattern (using *). Only one file can be returned,
    as order is arbitrary in glob.glob().
    Args:
        caps_directory: input caps directory
        information: dictionary containing all the relevant information to look for the files. Dict must contains the
                     following keys : pattern, description, needed_pipeline
                             pattern: define the pattern of the final file
                             description: string to describe what the file is
                             needed_pipeline (optional): string describing the pipeline needed to obtain the file beforehand
        raise_exception: if True (normal behavior), an exception is raised if errors happen. If not, we return the file
                        list as it is

    Returns:
          string of the found file

    Raises:
        ClinicaCAPSError if no file is found, or more than 1 files are found
    """
    from os.path import join

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaCAPSError

    assert isinstance(
        information, dict
    ), "A dict must be provided for the argument 'dict'"
    assert all(
        elem in information.keys()
        for elem in ["pattern", "description", "needed_pipeline"]
    ), "'information' must contain the keys 'pattern', 'description', 'needed_pipeline'"

    pattern = information["pattern"]
    # Some check on the formatting on the data
    assert pattern[0] != "/", (
        "pattern argument cannot start with char: / (does not work in os.path.join function). "
        "If you want to indicate the exact name of the file, use the format"
        " directory_name/filename.extension or filename.extension in the pattern argument."
    )

    check_caps_folder(caps_directory)

    current_pattern = join(caps_directory, "**/", pattern)
    current_glob_found = insensitive_glob(current_pattern, recursive=True)

    if len(current_glob_found) != 1 and raise_exception is True:
        error_string = f"{Fore.RED}\n[Error] Clinica encountered a problem while getting {information['description']}. "
        if len(current_glob_found) == 0:
            error_string += "No file was found"
        else:
            error_string += f"{len(current_glob_found)} files were found:"
            for found_files in current_glob_found:
                error_string += f"\n\t{found_files}"
            error_string += (
                f"{Fore.RESET}\n\tCAPS directory: {caps_directory}\n{Fore.YELLOW}"
                f"Please note that the following clinica pipeline(s) must have run to obtain these files: "
                f"{information['needed_pipeline']} {Fore.RESET}\n"
            )
        raise ClinicaCAPSError(error_string)
    return current_glob_found[0]


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


def fetch_file(remote, dirname=None):
    """Download a specific file and save it into the resources folder of the package.

    Args:
        remote: structure containing url, filename and checksum
        dirname: absolute path where the file will be downloaded

    Returns:
        file_path: absolute file path
    """
    import os.path
    import shutil
    import ssl
    from urllib.error import URLError
    from urllib.request import Request, urlopen

    from clinica.utils.stream import cprint

    if not os.path.exists(dirname):
        cprint("Path to the file does not exist")
        cprint("Stop Clinica and handle this error")

    file_path = os.path.join(dirname, remote.filename)
    # Download the file from `url` and save it locally under `file_name`:
    gcontext = ssl.SSLContext()
    req = Request(remote.url + remote.filename)
    try:
        response = urlopen(req, context=gcontext)
    except URLError as e:
        if hasattr(e, "reason"):
            cprint("We failed to reach a server.")
            cprint([f"Reason: {e.reason}"])
        elif hasattr(e, "code"):
            cprint("The server could not fulfill the request.")
            cprint([f"Error code: {e.code}"])
    else:
        try:
            with open(file_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)
        except OSError as err:
            cprint("OS error: {0}".format(err))

    checksum = _sha256(file_path)
    if remote.checksum != checksum:
        raise IOError(
            f"{file_path} has an SHA256 checksum ({checksum}) from expected "
            f"({remote.checksum}), file may be corrupted."
        )
    return file_path


def get_file_from_server(remote_file, cache_path=None):
    """Download file from server.

    Args:
        remote_file (str): RemoteFileStructure defined in clinica.utils.inputs
        cache_path (str): (default: ~/.cache/clinica/data)

    Returns:
        Path to downloaded file.
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
                f"Unable to download {remote_file.filename} from {remote_file.url}: {err}"
            )

    return local_file
