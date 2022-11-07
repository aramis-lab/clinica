import re
from typing import List


def zip_nii(in_file: str, same_dir: bool = False):
    import gzip
    import shutil
    from os import getcwd
    from os.path import abspath, join

    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [zip_nii(f, same_dir) for f in in_file]

    orig_dir, base, ext = split_filename(str(in_file))

    # Already compressed
    if ext[-3:].lower() == ".gz":
        return in_file
    # Not compressed

    out_file = abspath(join(orig_dir if same_dir else getcwd(), base + ext + ".gz"))

    with open(in_file, "rb") as f_in, gzip.open(out_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


def unzip_nii(in_file: str):
    from nipype.algorithms.misc import Gunzip
    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [unzip_nii(f) for f in in_file]

    _, base, ext = split_filename(in_file)

    # Not compressed
    if ext[-3:].lower() != ".gz":
        return in_file
    # Compressed
    gunzip = Gunzip(in_file=in_file)
    gunzip.run()
    return gunzip.aggregate_outputs().out_file


def save_participants_sessions(participant_ids, session_ids, out_folder, out_file=None):
    """Save <participant_ids> <session_ids> in <out_folder>/<out_file> TSV file."""
    import os

    import pandas

    from clinica.utils.stream import cprint

    assert len(participant_ids) == len(session_ids)

    os.makedirs(out_folder, exist_ok=True)

    if out_file:
        tsv_file = os.path.join(out_folder, out_file)
    else:
        tsv_file = os.path.join(out_folder, "participants.tsv")

    try:
        data = pandas.DataFrame(
            {
                "participant_id": participant_ids,
                "session_id": session_ids,
            }
        )
        data.to_csv(tsv_file, sep="\t", index=False, encoding="utf-8")
    except Exception as e:
        cprint(msg=f"Impossible to save {out_file} with pandas", lvl="error")
        raise e


def _raise_non_bids_or_caps_compliant_filename(filename: str) -> None:
    raise ValueError(
        f"Input filename {filename} is not in a BIDS or CAPS compliant format."
        " It does not contain the subject and session information."
    )


def get_subject_id(bids_or_caps_file: str) -> str:
    """Extract the subject ID from a BIDS or CAPS file path.

    The subject ID is defined as

        sub-<participant_id>_ses-<session_label>

    In other words, it is the concatenation of the subject and session labels,
    separated by "_".

    Parameters
    ----------
    bids_or_caps_file: str
        Path to a file from a BIDS or CAPS folder.

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
    ValueError: Input filename sub-01_ses-M000_trc-18FAV45_pet.nii.gz is not in a BIDS or CAPS compliant format. It does not contain the subject and session information.

    See also
    --------
    extract_image_ids
    """
    m = re.search(r"(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)", bids_or_caps_file)
    if not m:
        _raise_non_bids_or_caps_compliant_filename(bids_or_caps_file)
    subject_id = m.group(1) + "_" + m.group(2)

    return subject_id


def get_filename_no_ext(filename: str) -> str:
    """Get the filename without the extension.

    Parameters
    ----------
    filename: str
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
    from nipype.utils.filemanip import split_filename

    _, filename_no_ext, _ = split_filename(filename)

    return filename_no_ext


def extract_image_ids(bids_or_caps_files: List[str]) -> List[str]:
    """Extract the image IDs from a list of BIDS or CAPS files.

    .. warning::
        The image ID is the same as the subject ID from function
        `get_subject_id()` but is extracted from the filename rather
        than from the path. (See examples section).

    Parameters
    ----------
    bids_or_caps_files: List[str]
        List of file names from which to extract the image IDs.

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
        m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)", f)
        if not m:
            _raise_non_bids_or_caps_compliant_filename(f)
        id_bids_or_caps_files.append(m.group())

    return id_bids_or_caps_files


def extract_subjects_sessions_from_filename(bids_or_caps_files):
    """Extract subjects/sessions (e.g. ['sub-CLNC01', 'sub-CLNC01']/['ses-M00', 'ses-M18'] from `bids_or_caps_files`."""
    id_bids_or_caps_files = extract_image_ids(bids_or_caps_files)
    split = [image_id.split("_") for image_id in id_bids_or_caps_files]
    subject_ids = [p_id[0] for p_id in split]
    session_ids = [s_id[1] for s_id in split]
    return subject_ids, session_ids


def extract_crash_files_from_log_file(filename):
    """Extract crash files (*.pklz) from `filename`."""
    import os
    import re

    assert os.path.isfile(
        filename
    ), f"extract_crash_files_from_log_file: filename parameter is not a file ({filename})"

    log_file = open(filename, "r")
    crash_files = []
    for line in log_file:
        if re.match("(.*)crashfile:(.*)", line):
            crash_files.append(line.replace("\t crashfile:", "").replace("\n", ""))

    return crash_files


def read_participant_tsv(tsv_file):
    """Extract participant IDs and session IDs from TSV file.

    Raise:
        ClinicaException if tsv_file is not a file
        ClinicaException if participant_id or session_id column is missing from TSV file
    """
    import os

    import pandas as pd

    from clinica.utils.exceptions import ClinicaException

    if not os.path.isfile(tsv_file):
        raise ClinicaException(
            "The TSV file you gave is not a file.\n"
            "Error explanations:\n"
            f"\t- Clinica expected the following path to be a file: {tsv_file}\n"
            "\t- If you gave relative path, did you run Clinica on the good folder?"
        )
    ss_df = pd.read_csv(tsv_file, sep="\t")
    if "participant_id" not in list(ss_df.columns.values):
        raise ClinicaException(
            f"The TSV file does not contain participant_id column (path: {tsv_file})"
        )
    if "session_id" not in list(ss_df.columns.values):
        raise ClinicaException(
            f"The TSV file does not contain session_id column (path: {tsv_file})"
        )
    participants = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    # Remove potential whitespace in participant_id or session_id
    return [sub.strip(" ") for sub in participants], [
        ses.strip(" ") for ses in sessions
    ]


def extract_metadata_from_json(json_file, list_keys):
    """Extract fields from JSON file."""
    import datetime
    import json

    from clinica.utils.exceptions import ClinicaException

    list_values = []
    try:
        with open(json_file, "r") as file:
            data = json.load(file)
            for key in list_keys:
                list_values.append(data[key])
    except EnvironmentError:
        raise EnvironmentError(
            f"[Error] Clinica could not open the following JSON file: {json_file}"
        )
    except KeyError:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_message = f"[{now}] Error: Clinica could not find the e key in the following JSON file: {json_file}"
        raise ClinicaException(error_message)
    finally:
        file.close()

    return list_values
