from os import PathLike
from pathlib import Path
from typing import List, Optional

import pandas as pd


def create_subs_sess_list(
    input_dir: PathLike,
    output_dir: PathLike,
    file_name: Optional[str] = None,
    is_bids_dir: bool = True,
    use_session_tsv: bool = False,
):
    """Create the file subject_session_list.tsv that contains the list of the
    visits for each subject for a BIDS or CAPS compliant dataset.

    Parameters
    ----------
    input_dir : PathLike
        Path to the BIDS or CAPS directory.

    output_dir : PathLike
        Path to the output directory.

    file_name : str, optional
        The name of the output file.

    is_bids_dir : bool, optional
        Specify if input_dir is a BIDS directory or not (i.e. a CAPS directory).
        Default=True.

    use_session_tsv : bool
        Specify if the list uses the sessions listed in the sessions.tsv files.
        Default=False.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    path_to_search = input_dir if is_bids_dir else input_dir / "subjects"
    txt = _create_subs_sess_list_as_text(path_to_search, use_session_tsv)
    file_name = file_name or "subjects_sessions_list.tsv"
    (output_dir / file_name).write_text(txt)


def _create_subs_sess_list_as_text(path_to_search: Path, use_session_tsv: bool) -> str:
    txt = "participant_id\tsession_id\n"
    subjects_paths = [f for f in path_to_search.glob("*sub-*")]
    subjects_paths.sort()
    if len(subjects_paths) == 0:
        raise IOError("Dataset empty or not BIDS/CAPS compliant.")
    for subject_path in subjects_paths:
        if use_session_tsv:
            txt += _create_session_list_for_subject_from_tsv(subject_path)
        else:
            sessions = [f for f in subject_path.glob("*ses-*")]
            sessions.sort()
            for ses_path in sessions:
                txt += f"{subject_path.name}\t{ses_path.name}\n"
    return txt


def _create_session_list_for_subject_from_tsv(subject_path: Path) -> str:
    if not (subject_path / f"{subject_path.name}_sessions.tsv").exists():
        raise ValueError(
            f"In dataset located at {subject_path.parent}, there is no session "
            f"TSV file for subject {subject_path.name}. Consider setting the "
            "argument `use_session_tsv` to False in order to rely on folder "
            "names parsing."
        )
    session_df = pd.read_csv(
        subject_path / f"{subject_path.name}_sessions.tsv", sep="\t"
    )
    session_df.dropna(how="all", inplace=True)
    session_list = list(session_df["session_id"].to_numpy())
    return (
        "\n".join([f"{subject_path.name}\t{session}" for session in session_list])
        + "\n"
    )


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
