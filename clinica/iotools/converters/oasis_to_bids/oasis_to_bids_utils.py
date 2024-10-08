import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from clinica.iotools.bids_utils import StudyName, get_bids_sess_list
from clinica.utils.stream import cprint

__all__ = ["create_sessions_dict", "write_sessions_tsv"]


def create_sessions_dict(
    clinical_data_dir: Path,
    bids_dir: Path,
    clinical_specifications_folder: Path,
    bids_ids: list[str],
) -> dict:
    """Extract the information regarding the sessions and store them in a dictionary (session M000 only).

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the input folder.

    bids_dir : Path
        The path to the BIDS directory.

    clinical_specifications_folder : Path
        The path to the clinical file.

    bids_ids : list of str
        The list of bids ids.

    Returns
    -------
    dict :
        Session dict.
    """

    location = f"{StudyName.OASIS.value} location"
    sessions = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")
    sessions_fields = sessions[StudyName.OASIS.value]
    field_location = sessions[location]
    sessions_fields_bids = sessions["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""

            file_to_read_path = clinical_data_dir / location
            file_ext = os.path.splitext(location)[1]
            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
            elif file_ext == ".csv":
                file_to_read = pd.read_csv(file_to_read_path)
            else:
                raise ValueError(
                    f"Unknown file extension {file_ext}. Expecting either .xlsx or .csv."
                )

            for r in range(0, len(file_to_read.values)):
                # Extracts the subject ids columns from the dataframe
                subj_id = file_to_read.iloc[r]["ID"]
                if hasattr(subj_id, "dtype"):
                    if subj_id.dtype == np.int64:
                        subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3] + subj_id[5:9])

                # Extract the corresponding BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if subj_bids:
                    subj_bids = subj_bids[0]
                    subj_dir = bids_dir / subj_bids
                    session_names = get_bids_sess_list(subj_dir)
                    for s in session_names:
                        s_name = s.replace("ses-", "")
                        row = file_to_read.iloc[r]
                        if subj_bids not in sessions_dict:
                            sessions_dict.update({subj_bids: {}})
                        if s_name not in sessions_dict[subj_bids].keys():
                            sessions_dict[subj_bids].update({s_name: {"session_id": s}})
                        (sessions_dict[subj_bids][s_name]).update(
                            {sessions_fields_bids[i]: row[sessions_fields[i]]}
                        )

    return sessions_dict


def write_sessions_tsv(bids_dir: Path, sessions_dict: dict) -> None:
    """Create <participant_id>_sessions.tsv files.

    Basically writes the content of the function
    `clinica.iotools.bids_utils.create_sessions_dict` in several TSV files
    following the BIDS specification.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    sessions_dict : dict
        Dictionary containing sessions metadata.

        .. note::
            This is the output of the function
            `clinica.iotools.bids_utils.create_sessions_dict`.

    See also
    --------
    create_sessions_dict
    write_scans_tsv
    """
    for subject_path in bids_dir.glob("sub-*"):
        if subject_path.name in sessions_dict:
            session_df = pd.DataFrame.from_dict(
                sessions_dict[subject_path.name], orient="index"
            )
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
        else:
            print(f"No session data available for {subject_path}")
            session_df = pd.DataFrame(columns=["session_id"])
            session_df["session_id"] = pd.Series("M000")
        session_df = session_df.set_index("session_id").fillna("n/a")
        session_df.to_csv(
            subject_path / f"{subject_path.name}_sessions.tsv",
            sep="\t",
            encoding="utf8",
        )
