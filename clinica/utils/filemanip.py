import os
import re
from pathlib import Path, PurePath
from typing import Callable, List, Optional, Tuple, Union

import pandas as pd
from nipype.utils.filemanip import split_filename

from clinica.utils.exceptions import ClinicaException


def _zip_unzip_nii(
    in_file: Union[PurePath, List[PurePath]],
    func: Callable,
    same_dir: bool = False,
) -> Union[PurePath, List[PurePath]]:
    import gzip
    import operator
    import shutil
    from os import getcwd
    from os.path import abspath, join

    from traits.trait_base import _Undefined

    zipping = func.__name__ == "zip_nii"

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if isinstance(in_file, list):
        return [func(f, same_dir) for f in in_file]

    op = operator.eq if zipping else operator.ne
    if op(in_file.suffix, ".gz"):
        return in_file

    if not in_file.exists():
        raise FileNotFoundError(f"File {in_file} does not exist.")

    orig_dir, base, ext = split_filename(str(in_file))
    new_ext = ext + ".gz" if zipping else ext[:-3]
    out_file = abspath(join(orig_dir if same_dir else getcwd(), base + new_ext))

    outer = open if zipping else gzip.open
    inner = gzip.open if zipping else open
    with outer(in_file, "rb") as f_in:
        with inner(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return Path(out_file)


def zip_nii(
    in_file: Union[PurePath, List[PurePath]], same_dir: bool = False
) -> Union[PurePath, List[PurePath]]:
    """Compress the provided file(s).

    Parameters
    ----------
    in_file: PurePath or List[PurePath]
        File to be zipped, or list of files to be zipped.

    same_dir: bool, optional
        If True, the zip file is written in the same directory as the
        provided file. Otherwise, it is written in the current working
        directory. Default=False.

    Returns
    -------
    out_file: PurePath or List[PurePath]
        Path to the resulting zip file, or list of paths to the
        resulting zip files if a list was provided.

    Notes
    -----
    If the provided file is already zipped (that is, its extension
    is '.gz'), then nothing is done and the path to the provided file
    is returned.

    Raises
    ------
    FileNotFoundError
        If the provided file does not exist.
    """
    return _zip_unzip_nii(in_file, zip_nii, same_dir)


def unzip_nii(
    in_file: Union[PurePath, List[PurePath]],
    same_dir: bool = False,
) -> Union[PurePath, List[PurePath]]:
    """Decompress the provided file(s).

    Parameters
    ----------
    in_file: PurePath or List[PurePath]
        File to be zipped, or list of files to be zipped.

    same_dir: bool, optional
        If True, the unzipped file is written in the same directory as the
        provided zip file. Otherwise, it is written in the current working
        directory. Default=False.

    Returns
    -------
    out_file: PurePath or List[PurePath]
        Path to the resulting zip file, or list of paths to the
        resulting zip files if a list was provided.

    Notes
    -----
    If the provided file is not zipped (that is, its extension
    is different from '.gz'), then nothing is done and the path
    to the provided file is returned.

    Raises
    ------
    FileNotFoundError
        If the provided file does not exist.
    """
    return _zip_unzip_nii(in_file, unzip_nii, same_dir)


def save_participants_sessions(
    participant_ids: List[str],
    session_ids: List[str],
    out_folder: PurePath,
    out_file: Optional[str] = "participants.tsv",
) -> None:
    """Save the participants and sessions to TSV.

    Parameters
    ----------
    participant_ids: List[str]
        List of participant IDs.

    session_ids: List[str]
        List of session IDs.

    out_folder: PurePath
        Folder where the resulting TSV file should be written.

    out_file: str, optional
        Name of the TSV file. Default="participants.tsv".

    Raises
    ------
    ValueError
        If provided `participant_ids` and `session_ids` do not have
        the same length.
    """
    from clinica.utils.stream import cprint

    if len(participant_ids) != len(session_ids):
        raise ValueError(
            "The number of participant IDs is not equal to the number of session IDs."
        )

    out_folder = Path(out_folder)
    os.makedirs(out_folder, exist_ok=True)
    tsv_file = out_folder / out_file

    data = pd.DataFrame(
        {
            "participant_id": participant_ids,
            "session_id": session_ids,
        }
    )
    try:
        data.to_csv(tsv_file, sep="\t", index=False, encoding="utf-8")
    except Exception as e:
        cprint(msg=f"Impossible to save {out_file} with pandas", lvl="error")
        raise e


def _raise_non_bids_or_caps_compliant_filename(filename: PurePath) -> None:
    raise ValueError(
        f"Input filename {filename} is not in a BIDS or CAPS compliant format."
        " It does not contain the subject and session information."
    )


def get_subject_id(bids_or_caps_file: PurePath) -> str:
    """Extract the subject ID from a BIDS or CAPS file path.

    The subject ID is defined as

        sub-<participant_id>_ses-<session_label>

    In other words, it is the concatenation of the subject and session labels,
    separated by "_".

    Parameters
    ----------
    bids_or_caps_file: PurePath
        Path to a file from a BIDS or CAPS folder.

    Returns
    -------
    subject_id: str
        The subject ID corresponding to the given file.

    Examples
    --------
    >>> get_subject_id(Path("sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz"))
    'sub-01_ses-M000'
    >>> get_subject_id(Path("sub-01_ses-M000_trc-18FAV45_pet.nii.gz"))
    Traceback (most recent call last):
    ValueError: Input filename sub-01_ses-M000_trc-18FAV45_pet.nii.gz is not in a BIDS or CAPS compliant format. It does not contain the subject and session information.

    See also
    --------
    extract_image_ids
    """
    m = re.search(r"(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)", str(bids_or_caps_file))
    if not m:
        _raise_non_bids_or_caps_compliant_filename(bids_or_caps_file)
    subject_id = m.group(1) + "_" + m.group(2)

    return subject_id


def get_filename_no_ext(filename: PurePath) -> str:
    """Get the filename without the extension.

    Parameters
    ----------
    filename: PurePath
        The full filename from which to extract the extension out.

    Returns
    -------
    filename_no_ext: str
        The filename with extension removed.

    Examples
    --------
    >>> get_filename_no_ext(Path("foo.nii.gz"))
    'foo'
    >>> get_filename_no_ext(Path("sub-01/ses-M000/sub-01_ses-M000.tar.gz"))
    'sub-01_ses-M000'
    """
    _, filename_no_ext, _ = split_filename(str(filename))

    return filename_no_ext


def extract_image_ids(bids_or_caps_files: List[PurePath]) -> List[str]:
    """Extract the image IDs from a list of BIDS or CAPS files.

    .. warning::
        The image ID is the same as the subject ID from function
        `get_subject_id()` but is extracted from the filename rather
        than from the path. (See examples section).

    Parameters
    ----------
    bids_or_caps_files: List[PurePath]
        List of file paths from which to extract the image IDs.

    Returns
    -------
    id_bids_or_caps_files: List[str]
        List of extracted image IDs.

    Examples
    --------
    >>> extract_image_ids([Path("sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz")])
    ['sub-01_ses-M000']

    Beware, the image IDs is extracted from the filename and not from the folder names
    as it is the case for `get_subject_id()`:

    >>> extract_image_ids([Path("foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json")])
    ['sub-01_ses-M000']

    See also
    --------
    get_subject_id
    """
    id_bids_or_caps_files = []
    for f in bids_or_caps_files:
        m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)", str(f))
        if not m:
            _raise_non_bids_or_caps_compliant_filename(f)
        id_bids_or_caps_files.append(m.group())

    return id_bids_or_caps_files


def extract_subjects_sessions_from_filename(
    bids_or_caps_files: List[PurePath],
) -> Tuple[List[str], List[str]]:
    """Extract the subject and session labels from a list of BIDS or CAPS files.

    Parameters
    ----------
    bids_or_caps_files: List[PurePath]
        List of file paths for which to extract the subject and session labels.

    Returns
    -------
    subject_ids: List[str]
        List of subject labels in the same order as `bids_or_caps_files`.

    session_ids: List[str]
        List of session labels in the same order as `bids_or_caps_files`.

    Examples
    --------
    >>> extract_subjects_sessions_from_filename([Path("sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz"), Path("foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json", "sub-01_ses-M000.tar.gz")])
    (['sub-01', 'sub-01', 'sub-01'], ['ses-M000', 'ses-M000', 'ses-M000'])

    See also
    --------
    extract_image_ids
    """
    id_bids_or_caps_files = extract_image_ids(bids_or_caps_files)
    split = [image_id.split("_") for image_id in id_bids_or_caps_files]
    subject_ids = [p_id[0] for p_id in split]
    session_ids = [s_id[1] for s_id in split]
    return subject_ids, session_ids


def extract_crash_files_from_log_file(filename: PurePath) -> List[str]:
    """Extract crash files (*.pklz) from `filename`.

    Parameters
    ----------
    filename: PurePath
        Path to log file for which the crash files should be extracted.

    Returns
    -------
    crash_files: List[str]
        List of crash files.
    """
    filename = Path(filename)
    if not filename.is_file():
        raise ValueError(
            f"extract_crash_files_from_log_file: filename parameter is not a file ({filename})"
        )
    crash_files = []
    with open(filename, "r") as log_file:
        for line in log_file:
            if re.match("(.*)crashfile:(.*)", line):
                crash_files.append(line.replace("\t crashfile:", "").replace("\n", ""))

    return crash_files


def read_participant_tsv(tsv_file: PurePath) -> Tuple[List[str], List[str]]:
    """Extract participant IDs and session IDs from TSV file.

    Parameters
    ----------
    tsv_file: PurePath
        Participant TSV file from which to extract the participant and session IDs.

    Returns
    -------
    participants: List[str]
        List of participant IDs.

    sessions: List[str]
        List of session IDs.

    Raises
    ------
    ClinicaException
        If `tsv_file` is not a file.
        If `participant_id` or `session_id` column is missing from TSV file.

    Examples
    --------
    >>> df = pd.DataFrame({"participant_id": ["sub-01", "sub-01", "sub-02"], "session_id": ["ses-M000", "ses-M006", "ses-M000"]})
    >>> df.to_csv("participants.tsv", sep="\t")
    >>> read_participant_tsv(Path("participant.tsv"))
    (["sub-01", "sub-01", "sub-02"], ["ses-M000", "ses-M006", "ses-M000"])
    """
    tsv_file = Path(tsv_file)
    if not tsv_file.is_file():
        raise ClinicaException(
            "The TSV file you gave is not a file.\nError explanations:\n"
            f"\t- Clinica expected the following path to be a file: {tsv_file}\n"
            "\t- If you gave relative path, did you run Clinica on the good folder?"
        )

    df = pd.read_csv(tsv_file, sep="\t")
    for column in ("participant_id", "session_id"):
        if column not in list(df.columns.values):
            raise ClinicaException(
                f"The TSV file does not contain {column} column (path: {tsv_file})"
            )

    return (
        [sub.strip(" ") for sub in list(df.participant_id)],
        [ses.strip(" ") for ses in list(df.session_id)],
    )


def extract_metadata_from_json(json_file: PurePath, list_keys: List[str]) -> List[str]:
    """Extract fields from JSON file.

    Parameters
    ----------
    json_file: PurePath
        Path to JSON file from which to extract metadata.

    list_keys: List[str]
        List of keys to extract.

    Returns
    -------
    list_values: List[str]
        List of values extracted from the JSON file and corresponding to the keys.
    """
    import json

    json_file = Path(json_file)
    if not json_file.exists():
        raise FileNotFoundError(
            f"Clinica could not open the following JSON file: {json_file}"
        )
    with open(json_file, "r") as file:
        data = json.load(file)
    list_values = []
    missing_keys = []
    for key in list_keys:
        if key in data:
            list_values.append(data[key])
        else:
            missing_keys.append(key)
    if len(missing_keys) > 0:
        error_message = f"Clinica could not find the following keys in the following JSON file: {json_file}:\n- "
        error_message += "\n- ".join(missing_keys)
        raise ClinicaException(error_message)

    return list_values
