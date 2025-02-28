from pathlib import Path
from typing import Iterable, Union

import pandas as pd

from clinica.converters.study_models import StudyName, bids_id_factory

__all__ = ["create_sessions_df", "write_sessions_tsv", "write_scans_tsv"]


def _convert_cdr_to_diagnosis(cdr: Union[int, str]) -> str:
    if cdr == 0:
        return "CN"
    elif isinstance(cdr, int) and cdr > 0:
        return "AD"
    else:
        return "n/a"


def create_sessions_df(
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
    bids_ids: Iterable[str],
) -> pd.DataFrame:
    """Extract the information regarding sessions M000 and store them in a dataframe.

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the input folder.

    clinical_specifications_folder : Path
        The path to the clinical file folder.

    bids_ids : list of str
        The list of bids ids which are in the BIDS directory.

    Returns
    -------
    pd.Dataframe :
        Session df.
    """

    study = StudyName.OASIS.value
    location = f"{study} location"
    spec = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")[
        [study, location, "BIDS CLINICA"]
    ].dropna()

    sessions_df = pd.DataFrame()
    if len(spec[location].unique()) == 1:
        loc = spec[location].unique()[0]
    else:
        raise ValueError(
            f"OASIS1 metadata is supposed to be contained in only 1 file, {len(spec[location].unique())} were detected : {spec[location].unique()}"
        )

    file = pd.read_excel(clinical_data_dir / loc)
    file["BIDS ID"] = file.ID.apply(
        lambda x: bids_id_factory(StudyName.OASIS).from_original_study_id(x)
    )
    file.set_index("BIDS ID", drop=True, inplace=True)

    for _, row in spec[spec[location] == loc].iterrows():
        sessions_df[row["BIDS CLINICA"]] = file[row[[study]]]

    missing_subjects = set(bids_ids) - set(sessions_df.index)
    for ms in missing_subjects:
        sessions_df.loc[ms] = ["n/a" for _ in sessions_df.columns]

    sessions_df = sessions_df.loc[bids_ids]

    sessions_df["diagnosis"] = sessions_df["diagnosis"].apply(
        lambda x: _convert_cdr_to_diagnosis(x)
    )

    sessions_df.insert(loc=0, column="session_id", value="ses-M000")

    return sessions_df


def write_sessions_tsv(bids_dir: Path, sessions_df: pd.DataFrame) -> None:
    """Writes the content of the function `clinica.iotools.bids_utils.create_sessions_df`
    in several TSV files following the BIDS specification.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    sessions_df : DataFrame
        Contains sessions metadata.

        .. note::
            This is the output of the function
            `clinica.iotools.bids_utils.create_sessions_df`.

    See also
    --------
    create_sessions_df
    """
    for subject, data in sessions_df.iterrows():
        session_path = bids_dir / subject
        data.to_frame().T.to_csv(
            session_path / f"{subject}_sessions.tsv",
            sep="\t",
            encoding="utf8",
            index=False,
        )


def write_scans_tsv(bids_dir: Path) -> None:
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
