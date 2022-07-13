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


def get_subject_id(bids_or_caps_file: str) -> str:
    """Extract "sub-<participant_id>_ses-<session_label>" from BIDS or CAPS file."""
    import re

    m = re.search(r"(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)", bids_or_caps_file)

    if not m:
        raise ValueError(
            f"Input filename {bids_or_caps_file} is not in a BIDS or CAPS compliant format."
            " It does not contain the subject and session information."
        )

    subject_id = m.group(1) + "_" + m.group(2)

    return subject_id


def get_filename_no_ext(filename):
    """Get filename without extension [".nii.gz", ".tar.gz", ".niml.dset"]."""
    from nipype.utils.filemanip import split_filename

    _, filename_no_ext, _ = split_filename(filename)

    return filename_no_ext


def extract_image_ids(bids_or_caps_files):
    """Extract image IDs (e.g. ['sub-CLNC01_ses-M00', 'sub-CLNC01_ses-M18']  from `bids_or_caps_files`."""
    import re

    id_bids_or_caps_files = [
        re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)", file).group()
        for file in bids_or_caps_files
    ]
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
    except KeyError as e:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_message = f"[{now}] Error: Clinica could not find the e key in the following JSON file: {json_file}"
        raise ClinicaException(error_message)
    finally:
        file.close()

    return list_values
