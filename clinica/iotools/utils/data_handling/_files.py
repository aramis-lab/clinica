from os import PathLike
from pathlib import Path
from typing import List, Optional

import pandas as pd


def create_subs_sess_list(
    input_dir: str,
    output_dir: str,
    file_name: Optional[str] = None,
    is_bids_dir: bool = True,
    use_session_tsv: bool = False,
):
    """Create the file subject_session_list.tsv that contains the list of the
    visits for each subject for a BIDS or CAPS compliant dataset.

    Args:
        input_dir (str): Path to the BIDS or CAPS directory.
        output_dir (str): Path to the output directory
        file_name: name of the output file
        is_bids_dir (boolean): Specify if input_dir is a BIDS directory or
            not (i.e. a CAPS directory)
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files
    """
    import os
    from glob import glob
    from os import path

    os.makedirs(output_dir, exist_ok=True)

    if not file_name:
        file_name = "subjects_sessions_list.tsv"
    subjs_sess_tsv = open(path.join(output_dir, file_name), "w")
    subjs_sess_tsv.write("participant_id" + "\t" + "session_id" + "\n")

    if is_bids_dir:
        path_to_search = input_dir
    else:
        path_to_search = path.join(input_dir, "subjects")
    subjects_paths = glob(path.join(path_to_search, "*sub-*"))

    # Sort the subjects list
    subjects_paths.sort()

    if len(subjects_paths) == 0:
        raise IOError("Dataset empty or not BIDS/CAPS compliant.")

    for sub_path in subjects_paths:
        subj_id = sub_path.split(os.sep)[-1]

        if use_session_tsv:
            session_df = pd.read_csv(
                path.join(sub_path, subj_id + "_sessions.tsv"), sep="\t"
            )
            session_df.dropna(how="all", inplace=True)
            session_list = list(session_df["session_id"].to_numpy())
            for session in session_list:
                subjs_sess_tsv.write(subj_id + "\t" + session + "\n")

        else:
            sess_list = glob(path.join(sub_path, "*ses-*"))

            for ses_path in sess_list:
                session_name = ses_path.split(os.sep)[-1]
                subjs_sess_tsv.write(subj_id + "\t" + session_name + "\n")

    subjs_sess_tsv.close()


def write_list_of_files(file_list: List[PathLike], output_file: PathLike) -> Path:
    """Save `file_list` list of files into `output_file` text file.

    Parameters
    ----------
    file_list : list of PathLike objects
        List of path to files to write.

    output_file : PathLike
        Path to the output txt file.

    Returns
    -------
    output_file : PathLike
        The path to the output file.

    Raises
    ------
    TypeError :
        If something else than a list was provided for file_list.
    IOError :
        If output_file already exists.
    """
    output_file = Path(output_file)
    if not isinstance(file_list, list):
        raise TypeError(
            f"`file_list` argument must be a list of paths. Instead {type(file_list)} was provided."
        )
    if output_file.is_file():
        raise IOError(f"Output file {output_file} already exists.")

    with open(output_file, "w") as fp:
        for created_file in file_list:
            fp.write(f"{created_file}\n")

    return output_file
