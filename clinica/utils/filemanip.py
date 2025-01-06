from os import PathLike
from pathlib import Path
from typing import Callable, List, Optional, Union

from .bids import Visit

__all__ = [
    "UserProvidedPath",
    "delete_directories",
    "delete_directories_task",
    "extract_crash_files_from_log_file",
    "extract_image_ids",
    "extract_visits",
    "extract_metadata_from_json",
    "extract_subjects_sessions_from_filename",
    "get_filename_no_ext",
    "get_parent",
    "get_subject_id",
    "load_volume",
    "save_participants_sessions",
    "unzip_nii",
    "zip_nii",
]


UserProvidedPath = Union[str, PathLike]


def zip_nii(in_file: str, same_dir: bool = False) -> str:
    """Compress the provided file(s).

    .. note::
        If a list of file paths is provided, then this function
        will be called on each item and a list of paths to the
        resulting zip files will be returned.
        The type hint is not reflecting this fact because of the
        self-contained requirement from Nipype which is preventing
        the usage of types outside of Python's builtins.

    Parameters
    ----------
    in_file: str
        Path to file to be zipped as a string.

    same_dir: bool, optional
        If True, the zip file is written in the same directory as the
        provided file. Otherwise, it is written in the current working
        directory. Default=False.

    Returns
    -------
    out_file: str
        Path to the resulting zip file.

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
    from clinica.utils.filemanip import _zip_unzip_nii  # noqa

    return _zip_unzip_nii(in_file, same_dir, compress=True)


def unzip_nii(
    in_file: str,
    same_dir: bool = False,
) -> str:
    """Decompress the provided file(s).

    .. note::
        If a list of file paths is provided, then this function
        will be called on each item and a list of paths to the
        resulting unzipped files will be returned.
        The type hint is not reflecting this fact because of the
        self-contained requirement from Nipype which is preventing
        the usage of types outside of Python's builtins.

    Parameters
    ----------
    in_file: str
        Path to file to be zipped.

    same_dir: bool, optional
        If True, the unzipped file is written in the same directory as the
        provided zip file. Otherwise, it is written in the current working
        directory. Default=False.

    Returns
    -------
    out_file: str
        Path to the resulting unzipped file.

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
    from clinica.utils.filemanip import _zip_unzip_nii  # noqa

    return _zip_unzip_nii(in_file, same_dir, compress=False)


def _zip_unzip_nii(in_file: str, same_dir: bool, compress: bool):
    import gzip
    import operator
    import shutil
    from os import getcwd
    from os.path import abspath, join
    from pathlib import Path

    from nipype.utils.filemanip import split_filename

    try:
        # Assuming in_file is path-like.
        in_file = Path(in_file)
    except TypeError:
        try:
            # Assuming in_file is a sequence type.
            return [_zip_unzip_nii(f, same_dir, compress) for f in in_file]
        except TypeError:
            # All other cases.
            return None

    op = operator.eq if compress else operator.ne
    if op(in_file.suffix, ".gz"):
        return str(in_file)

    if not in_file.exists():
        raise FileNotFoundError(f"File {in_file} does not exist.")

    orig_dir, base, ext = split_filename(str(in_file))
    new_ext = ext + ".gz" if compress else ext[:-3]
    out_file = abspath(join(orig_dir if same_dir else getcwd(), base + new_ext))

    outer = open if compress else gzip.open
    inner = gzip.open if compress else open
    with outer(in_file, "rb") as f_in:
        with inner(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return out_file


def load_volume(image_path: str):
    """Load a 3D nifti image from its path.

    If the image is 4D with a dummy fourth dimension,
    then "squeeze" the image into a proper 3D image.

    Parameters
    ----------
    image_path : str
        Path to the image to load.

    Returns
    -------
    img : Nifti1Image
        The loaded 3D image.

    Raises
    ------
    ValueError
        If the loaded image isn't 3D.
    """
    import copy

    import nibabel as nib

    img = nib.load(image_path)
    dim = len(img.shape)
    if dim != 3:
        if dim == 4 and img.shape[3] == 1:
            data = img.get_fdata()
            klass = img.__class__
            header = copy.deepcopy(img.header)
            return klass(data[:, :, :, 0], img.affine, header=header)
        raise ValueError(f"The image is not 3D but {dim}D.")
    return img


def save_participants_sessions(
    participant_ids: list[str],
    session_ids: list[str],
    out_folder: str,
    out_file: Optional[str] = "participants.tsv",
) -> None:
    """Save the participants and sessions to TSV.

    Parameters
    ----------
    participant_ids: List[str]
        List of participant IDs.

    session_ids: List[str]
        List of session IDs.

    out_folder: str
        Folder where the resulting TSV file should be written.

    out_file: str, optional
        Name of the TSV file. Default="participants.tsv".

    Raises
    ------
    ValueError
        If provided `participant_ids` and `session_ids` do not have
        the same length.
    """
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    if len(participant_ids) != len(session_ids):
        raise ValueError(
            "The number of participant IDs is not equal to the number of session IDs."
        )

    out_folder = Path(out_folder)
    out_folder.mkdir(parents=True, exist_ok=True)
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


def _check_bids_or_caps_compliance(filename: str, sep: str):
    import re

    m = re.search(sep.join([r"(sub-[a-zA-Z0-9]+)", r"(ses-[a-zA-Z0-9]+)"]), filename)
    if not m:
        raise ValueError(
            f"Input filename {filename} is not in a BIDS or CAPS compliant format."
            " It does not contain the subject and session information."
        )

    return m


def get_subject_id(bids_or_caps_file: Union[str, Path]) -> str:
    """Extract the subject ID from a BIDS or CAPS file path.

    The subject ID is defined as

        sub-<participant_id>_ses-<session_label>

    In other words, it is the concatenation of the subject and session labels,
    separated by "_".

    Parameters
    ----------
    bids_or_caps_file: str or Path
        The path to a file from a BIDS or CAPS folder.

    Returns
    -------
    subject_id: str
        The subject ID corresponding to the given file.

    Examples
    --------
    >>> get_subject_id("sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz")
    'sub-01_ses-M000'
    >>> get_subject_id("sub-01_ses-M000_trc-18FAV45_pet.nii.gz")
    Traceback (most recent call last):
    ValueError: Input filename sub-01_ses-M000_trc-18FAV45_pet.nii.gz is not in a BIDS or CAPS compliant format. It does not contain the subject and session information.  # noqa

    See also
    --------
    extract_image_ids
    """
    match = _check_bids_or_caps_compliance(str(bids_or_caps_file), sep="/")
    subject_id = match.group(1) + "_" + match.group(2)

    return subject_id


def get_filename_no_ext(filename: Union[str, Path]) -> str:
    """Get the filename without the extension.

    Parameters
    ----------
    filename: str or Path
        The full filename from which to extract the extension out.

    Returns
    -------
    filename_no_ext: str
        The filename with extension removed.

    Examples
    --------
    >>> get_filename_no_ext("foo.nii.gz")
    'foo'
    >>> get_filename_no_ext("sub-01/ses-M000/sub-01_ses-M000.tar.gz")
    'sub-01_ses-M000'
    """
    if not isinstance(filename, Path):
        filename = Path(filename)
    stem = filename.stem
    while "." in stem:
        stem = Path(stem).stem

    return stem


def extract_image_ids(bids_or_caps_files: list[str]) -> list[str]:
    """Extract the image IDs from a list of BIDS or CAPS files.

    .. warning::
        The image ID is the same as the subject ID from function
        `get_subject_id()` but is extracted from the filename rather
        than from the path. (See examples section).

    Parameters
    ----------
    bids_or_caps_files: List[str]
        List of file paths from which to extract the image IDs.

    Returns
    -------
    id_bids_or_caps_files: List[str]
        List of extracted image IDs.

    Examples
    --------
    >>> extract_image_ids(["sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz"])
    ['sub-01_ses-M000']

    Beware, the image IDs is extracted from the filename and not from the folder names
    as it is the case for `get_subject_id()`:

    >>> extract_image_ids(["foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json"])
    ['sub-01_ses-M000']

    See also
    --------
    get_subject_id
    """
    id_bids_or_caps_files = []
    for f in bids_or_caps_files:
        match = _check_bids_or_caps_compliance(f, sep="_")
        id_bids_or_caps_files.append(match.group())

    return id_bids_or_caps_files


def extract_visits(bids_or_caps_files: list[str]) -> list[Visit]:
    return [
        Visit(*image_id.split("_"))
        for image_id in extract_image_ids(bids_or_caps_files)
    ]


def extract_subjects_sessions_from_filename(
    bids_or_caps_files: list[str],
) -> tuple[list[str], list[str]]:
    """Extract the subject and session labels from a list of BIDS or CAPS files.

    Parameters
    ----------
    bids_or_caps_files: List[str]
        List of file paths for which to extract the subject and session labels.

    Returns
    -------
    subject_ids: List[str]
        List of subject labels in the same order as `bids_or_caps_files`.

    session_ids: List[str]
        List of session labels in the same order as `bids_or_caps_files`.

    Examples
    --------
    >>> extract_subjects_sessions_from_filename([
    ...     "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
    ...     "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json", "sub-01_ses-M000.tar.gz",
    ...])
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


def extract_crash_files_from_log_file(filename: str) -> list[str]:
    """Extract crash files (*.pklz) from `filename`.

    Parameters
    ----------
    filename: str
        Path to log file for which the crash files should be extracted.

    Returns
    -------
    crash_files: List[str]
        List of crash files.
    """
    import re
    from pathlib import Path

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


def extract_metadata_from_json(
    json_file: Union[str, Path],
    list_keys: list[str],
    handle_missing_keys: Optional[Callable] = None,
) -> list[str]:
    """Extract fields from JSON file.

    Parameters
    ----------
    json_file: str or Path
        Path to a json file containing metadata needed for the pipeline.

    list_keys: list of str
        List of fields the users wants to obtain from the json.

    handle_missing_keys: Optional[Callable]
        Function to use to handle the fields of the json that may be missing.

    Returns
    -------
    list of str:
        Contains the values for the requested fields.
    """
    import json

    from clinica.utils.exceptions import ClinicaException

    try:
        with open(json_file, "r") as file:
            data = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Clinica could not open the following JSON file: {json_file}"
        )
    missing_keys = set(list_keys).difference(set(data.keys()))
    if len(missing_keys) > 0 and handle_missing_keys is None:
        raise ClinicaException(
            f"Clinica could not find the following keys in the following JSON file: {missing_keys}."
        )
    if handle_missing_keys:
        patch = handle_missing_keys(data, missing_keys)
        data.update(patch)
    return [data[k] for k in list_keys]


def get_parent(path: str, n: int = 1) -> Path:
    """Get the path to the nth parent.

    Parameters
    ----------
    path: str
        path to a file.
    n: int
        depth we want to go up in the parents.

    Returns
    -------
    Path
        Path to a parent directory.

    Examples
    --------
    >>> get_parent('/path/to/a/file', 2)
    /path/to
    """
    if n <= 0:
        return Path(path)
    return get_parent(Path(path).parent, n - 1)


def _get_folder_size(folder: Union[str, Path]) -> int:
    """Compute the size in bytes recursively of the given folder.

    Parameters
    ----------
    folder : str or Path
        The path to the folder for which to compute size.

    Returns
    -------
    int :
        The size of the folder in bytes.

    Examples
    --------
    >>> _get_folder_size("./test/instantiation/")
    52571
    """
    import os
    from functools import partial

    prepend = partial(os.path.join, folder)
    return sum(
        [
            (os.path.getsize(f) if os.path.isfile(f) else _get_folder_size(f))
            for f in map(prepend, os.listdir(folder))
        ]
    )


def _get_folder_size_human(folder: Union[str, Path]) -> str:
    """Computes the size of the given folder in human-readable form.

    Parameters
    ----------
    folder : str
        Path to the folder for which to compute size.

    Returns
    -------
    str :
        The size in human-readable form.

    Examples
    --------
    >>> _get_folder_size_human("./test/instantiation/")
    '51.3388671875 KB'
    """
    return _humanize_bytes(_get_folder_size(folder))


def _humanize_bytes(size: int) -> str:
    """Convert a number of bytes in a human-readable form.

    Parameters
    ----------
    size : int
        The number of bytes to convert.

    Returns
    -------
    str :
        The number converted in the best unit for easy reading.
    """
    units = ("B", "KB", "MB", "GB", "TB")

    for unit in units[:-1]:
        if size < 1024.0:
            return f"{size} {unit}"
        size /= 1024.0

    return f"{size} {units[-1]}"


def delete_directories(directories: list[Union[str, Path]]) -> None:
    """This function deletes the directories of the given list".

    Parameters
    ----------
    directories : list of str
        Names of the directories we want to delete.
    """
    import shutil

    total_size: int = 0
    for directory in directories:
        total_size += _get_folder_size(str(directory))
        size = _get_folder_size_human(str(directory))
        shutil.rmtree(directory)
        _print_and_warn(f"Folder {directory} deleted. Freeing {size} of disk space...")
    _print_and_warn(f"Was able to remove {_humanize_bytes(total_size)} of data.")


def _print_and_warn(msg: str, lvl: str = "info") -> None:
    """Print the given message with the given level and warns with the same message."""
    import warnings

    from clinica.utils.stream import cprint

    cprint(msg=msg, lvl=lvl)
    warnings.warn(msg)


def delete_directories_task(directories: list) -> None:
    """Task for Nipype."""
    from clinica.utils.filemanip import delete_directories  # noqa

    return delete_directories(directories)
