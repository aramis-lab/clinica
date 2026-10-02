from pathlib import Path
from typing import Iterable

import pandas as pd
from clinicaio import BIDSDataset

from clinica.converters.study_models import OASIS2BIDSSubjectID
from clinica.utils.stream import cprint

__all__ = [
    "read_clinical_data",
    "read_imaging_data",
    "intersect_data",
    "split_clinical_data",
    "populate_bids_with_info",
    "convert_analyze_to_nifti",
]

"""

OASIS-2 dataset structure expected on disk
------------------------------------------
<imaging_dir>/
    OAS2_XXXX_MR1/
        RAW/
            mpr-1.nifti.img
            mpr-1.nifti.hdr
            mpr-2.nifti.img
            mpr-2.nifti.hdr
            ...
    OAS2_XXXX_MR2/
        RAW/
            ...

<clinical_dir>/
    oasis_longitudinal_demographics.csv   # or .xlsx
    (downloaded from https://sites.wustl.edu/oasisbrains/)

Clinical file columns
------------------------------------------
Subject ID | MRI ID | Group | Visit | MR Delay | M/F | Hand | Age | EDUC | SES |
MMSE | CDR | eTIV | nWBV | ASF
"""

# ---------------------------------------------------------------------------
# Required columns in the OASIS-2 demographics CSV
# ---------------------------------------------------------------------------
_REQUIRED_CLINICAL_COLUMNS = {
    "Subject ID",
    "MRI ID",
    "Group",
    "M/F",
    "Hand",
    "Age",
    "EDUC",
    "SES",
    "MMSE",
    "CDR",
    "eTIV",
    "nWBV",
    "ASF",
}


# ---------------------------------------------------------------------------
# Reading raw data
# ---------------------------------------------------------------------------


def read_clinical_data(clinical_data_directory: Path) -> pd.DataFrame:
    """Read the OASIS-2 longitudinal demographics CSV file.

    The OASIS-2 website distributes the demographics file as an XLSX.
    Both .xlsx and .csv are supported; the directory must contain exactly
    one file of either format.

    Parameters
    ----------
    clinical_data_directory:
        Directory that contains the ``oasis_longitudinal_demographics.csv``
        (or any single ``*.csv`` file).

    Returns
    -------
    pd.DataFrame
        One row per MRI session with all clinical variables.
    """
    xlsx_files = list(clinical_data_directory.glob("*.xlsx"))
    if not xlsx_files:
        raise FileNotFoundError(
            f"No clinical data file found in {clinical_data_directory}.\n"
            "Please place the OASIS-2 demographics file (.csv or .xlsx) "
            "in that directory."
        )
    if len(xlsx_files) > 1:
        raise ValueError(
            f"Multiple clinical data files found in {clinical_data_directory}: "
            f"{[f.name for f in xlsx_files]}. "
            "Only one file is expected."
        )

    df = pd.read_excel(xlsx_files[0])
    df.columns = df.columns.str.strip()

    missing = _REQUIRED_CLINICAL_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Clinical file is missing expected columns: {sorted(missing)}.\n"
            "Make sure you are providing the OASIS-2 longitudinal demographics file."
        )
    return df


def read_imaging_data(imaging_data_directory: Path) -> pd.DataFrame:
    """Discover all raw T1w ANALYZE images in the OASIS-2 directory tree.

    Parameters
    ----------
    imaging_data_directory:
        Root of the OASIS-2 imaging directory (contains ``OAS2_XXXX_MRY/``
        session folders).

    Returns
    -------
    pd.DataFrame
        One row per ``.img`` file with columns:
        ``source_path``, ``Subject ID``, ``session_label``, ``run_number``,
        ``participant_id``, ``session_id``.
    """
    from clinica.utils.stream import cprint

    records = []
    for img_path in _find_imaging_data(imaging_data_directory):
        # img_path is relative to imaging_data_directory
        # e.g.  OAS2_0001_MR1/RAW/mpr-1.nifti.img
        mri_id = img_path.parts[0]
        tokens = mri_id.split("_")  # ['OAS2', '0001', 'MR1']
        if len(tokens) < 3:
            cprint(
                f"Invalid image located at {img_path}. Will be skipped for conversion.",
                lvl="warning",
            )
            continue

        records.append(
            {
                "source_path": img_path,
                "participant_id": OASIS2BIDSSubjectID.from_original_study_id(
                    f"{tokens[0]}_{tokens[1]}"
                ),
                "session_id": "ses-" + tokens[2],  # ses-MR1 | MR2 | ...
                "run_number": _identify_run(img_path.name),
                "MRI ID": mri_id,
            }
        )

    if not records:
        raise FileNotFoundError(
            f"No ANALYZE (.img) T1w acquisitions found under {imaging_data_directory}.\n"
            "Expected files matching the pattern: "
            "OAS2_XXXX_MRY/RAW/mpr-N.nifti.img"
        )

    return (
        pd.DataFrame(records)
        .drop_duplicates()
        .sort_values(by=["source_path"])
        .reset_index(drop=True)
    )


def _find_imaging_data(path_to_source_data: Path) -> Iterable[Path]:
    """Yield relative paths to every raw MPRAGE ANALYZE acquisition.

    Files follow the naming convention: mpr-N.nifti.img
    e.g. mpr-1.nifti.img, mpr-2.nifti.img, ...
    """
    for image in path_to_source_data.rglob("OAS2_*_MR*/RAW/mpr-*.img"):
        yield image.relative_to(path_to_source_data)


def _identify_run(image_file_name: str) -> str:
    """Return a BIDS run number from the mpr-N suffix in the filename.

    Examples
    --------
    ``mpr-1.nifti.img``  ->  ``01``
    ``mpr-4.nifti.img``  ->  ``04``
    """
    import re

    match = re.search(r"mpr-(\d+)", image_file_name)
    return f"{int(match.group(1)):02d}" if match else "01"


# ---------------------------------------------------------------------------
# Merging imaging and clinical data
# ---------------------------------------------------------------------------


def intersect_data(
    df_imaging: pd.DataFrame,
    df_clinical: pd.DataFrame,
) -> pd.DataFrame:
    """Join the imaging inventory with clinical metadata.

    Parameters
    ----------
    df_imaging:
        Output of :func:`read_imaging_data`.
    df_clinical:
        Output of :func:`read_clinical_data`.

    Returns
    -------
    df_merged:
        One row per image file, enriched with all clinical variables and
        the target BIDS ``filename``.
    """
    df_clinical = df_clinical.copy()
    # Inner join: keep only subjects/sessions (info in MRI ID OASX_YYYY_MRZ) present in both data sources
    df_merged = df_imaging.merge(
        df_clinical,
        on=["MRI ID"],
        how="inner",
    )

    df_merged = df_merged.assign(
        filename=lambda df: df.apply(
            lambda row: (
                f"anat/{row.participant_id}_{row.session_id}_run-{row.run_number}_T1w.nii.gz"
            ),
            axis=1,
        )
    )
    return df_merged


# ---------------------------------------------------------------------------
# Building BIDS metadata tables
# ---------------------------------------------------------------------------


def split_clinical_data(
    df_merged: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the three BIDS metadata tables from the merged data.

    Returns
    -------
    participants : pd.DataFrame  (index = participant_id)
    sessions     : pd.DataFrame  (index = [participant_id, session_id])
    scans        : pd.DataFrame  (index = BIDS filename)
    """
    return (
        _build_participants_df(df_merged),
        _build_sessions_df(df_merged),
        _build_scans_df(df_merged),
    )


def _build_participants_df(df_merged: pd.DataFrame) -> pd.DataFrame:
    """One row per subject: sex, handedness, education, socioeconomic status."""
    return (
        df_merged[["participant_id", "M/F", "Hand", "EDUC", "SES"]]
        .rename(
            columns={
                "M/F": "sex",
                "Hand": "handedness",
                "EDUC": "education",
                "SES": "socioeconomic_status",
            }
        )
        .drop_duplicates()
        .set_index("participant_id", verify_integrity=True)
    )


def _build_sessions_df(df_merged: pd.DataFrame) -> pd.DataFrame:
    """One row per (subject, session): clinical variables measured at that visit."""
    # MR Delay = days since first visit; Visit = session label in source CSV.
    # Both are optional: not all exports include them.
    optional_cols = {"MR Delay": "days_since_first_visit", "Visit": "visit"}
    base_cols = {
        "Age": "age",
        "MMSE": "mmse",
        "CDR": "cdr",
        "eTIV": "etiv",
        "nWBV": "nwbv",
        "ASF": "asf",
        "Group": "group",
    }
    rename_map = {
        k: v
        for k, v in {**base_cols, **optional_cols}.items()
        if k in df_merged.columns
    }
    cols = ["participant_id", "session_id"] + list(rename_map.keys())
    return (
        df_merged[cols]
        .rename(columns=rename_map)
        .drop_duplicates(subset=["participant_id", "session_id"])
        .set_index(["participant_id", "session_id"], verify_integrity=True)
    )


def _build_scans_df(df_merged: pd.DataFrame) -> pd.DataFrame:
    """One row per image file: maps BIDS filename -> source ANALYZE path."""
    return df_merged[
        ["participant_id", "session_id", "source_path", "filename", "run_number"]
    ].set_index(["participant_id", "session_id"])


# ---------------------------------------------------------------------------
# ANALYZE → NIfTI conversion
# ---------------------------------------------------------------------------


def convert_analyze_to_nifti(source_path: Path, target_path: Path) -> None:
    """Convert an ANALYZE 7.5 .img/.hdr pair to a compressed NIfTI file.

    Some OASIS-2 ANALYZE files are stored as 4D arrays with a dummy
    trailing dimension of size 1, e.g. shape (256, 256, 128, 1).
    BIDS requires exactly 3 dimensions for T1w files, so any such
    trailing dimensions are squeezed out before saving.

    Parameters
    ----------
    source_path:
        Absolute path to the ``.img`` file (the ``.hdr`` must be alongside it).
    target_path:
        Destination path for the output ``.nii.gz`` file.
        Parent directories are created automatically.
    """
    import nibabel as nib
    import numpy as np

    img = nib.load(str(source_path))
    data = np.squeeze(np.asarray(img.dataobj))

    if data.ndim != 3:
        raise ValueError(
            f"Expected a 3D volume after squeezing but got shape {data.shape} "
            f"for file {source_path}."
        )

    nii = nib.Nifti1Image(data, img.affine, img.header)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nii, str(target_path))


# ---------------------------------------------------------------------------
# WRITE BIDS
# ---------------------------------------------------------------------------


def populate_bids_with_info(
    bids_dataset: BIDSDataset,
    participants: pd.DataFrame,
    sessions: pd.DataFrame,
    scans: pd.DataFrame,
) -> None:
    """
        Populates the BIDSDataset structure with information from participants, sessions and scans.
    After being applied each image is described and thus accessible without having to loop over each
    participant and session.

    Parameters
    ----------
    bids_dataset : BIDSDataset
    participants : pd.DataFrame
        Contains information that should go into participants.tsv. Indexed on participants.
    sessions : pd.DataFrame
        Contains information that should go into sessions.tsv. Indexed on participants and sessions.
    scans : pd.DataFrame
        Contains information that should go into scans.tsv. Indexed on participants and sessions.
    """
    # todo : chance that this would be useful for several datasets
    from clinicaio import DataType, FileExtension, ImageScanInfo

    for participant in participants.index:
        cprint(f"Converting OASIS2 subject {participant} to BIDS", lvl="debug")

        subject = bids_dataset.add_subject(participant)

        for session in sessions.loc[participant].index:
            ses = subject.add_session(session)

            for _, image in scans.loc[participant, session].iterrows():
                ses.write_image(
                    data_type=DataType.ANAT,
                    nifti_extension=FileExtension.NII_GZ,
                    entities={"run": image.run_number},
                    suffix="T1w",
                    scan_info=ImageScanInfo(
                        {"source_path": image.source_path, "run": image.run_number}
                    ),
                )

        subject.populate_sessions_info_from_df(sessions.reset_index(drop=False))

    bids_dataset.populate_subjects_info_from_df(participants.reset_index(drop=False))
