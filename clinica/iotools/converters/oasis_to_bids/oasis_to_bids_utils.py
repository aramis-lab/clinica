from pathlib import Path
from typing import Iterable

import pandas as pd

from clinica.iotools.bids_utils import StudyName, bids_id_factory

__all__ = ["create_sessions_dict", "write_sessions_tsv", "write_scans"]


def create_sessions_dict(
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
    bids_ids: Iterable[str],
) -> dict:
    """Extract the information regarding the sessions and store them in a dictionary (session M000 only).

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the input folder.

    clinical_specifications_folder : Path
        The path to the clinical file folder.

    bids_ids : list of str
        The list of bids ids.

    Returns
    -------
    dict :
        Session dict.
    """

    study = StudyName.OASIS.value
    location = f"{study} location"
    spec = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")[
        [study, location, "BIDS CLINICA"]
    ].dropna()
    sessions_dict = {}

    for loc in spec[location].unique():
        file = pd.read_excel(clinical_data_dir / loc)
        file["BIDS ID"] = file.ID.apply(
            lambda x: bids_id_factory(StudyName.OASIS).from_original_study_id(x)
        )
        file.set_index("BIDS ID", drop=True, inplace=True)
        sessions_df = pd.DataFrame()
        for _, row in spec[spec[location] == loc].iterrows():
            sessions_df[row["BIDS CLINICA"]] = file[row[[study]]]

        sessions_df = sessions_df.loc[bids_ids]
        sessions_df["diagnosis"] = sessions_df["diagnosis"].apply(
            lambda x: "AD" if x > 0 else "CN"
        )
        sessions_df["session_id"] = "ses-M000"

        for bids_id, row in sessions_df.iterrows():
            sessions_dict.update(
                {bids_id: {"M000": {label: value for label, value in row.items()}}}
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


def write_scans(
    bids_dir: Path,
) -> None:
    """
    Write the scans.tsv file at the root of baseline sessions (ses-M000).

    Parameters
    ----------
    bids_dir : Path to the BIDS output

    """

    for subject_path in bids_dir.rglob("sub-*"):
        if subject_path.is_dir():
            to_write = pd.DataFrame(
                {
                    "filename": [
                        f"{path.parent.name}/{path.name}"
                        for path in subject_path.rglob("*ses-M000*.nii.gz")
                    ]
                }
            )

            to_write.to_csv(
                subject_path / "ses-M000" / f"{subject_path.name}_ses-M000_scans.tsv",
                sep="\t",
                index=False,
            )
