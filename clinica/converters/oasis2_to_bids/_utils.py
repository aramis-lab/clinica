from pathlib import Path
from typing import Iterable

import pandas as pd

__all__ = [
    "read_clinical_data",
    "read_imaging_data",
    "intersect_data",
    "dataset_to_bids",
    "write_bids",
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
    csv_files = list(clinical_data_directory.glob("*.csv"))
    xlsx_files = list(clinical_data_directory.glob("*.xlsx"))

    all_files = csv_files + xlsx_files
    if not all_files:
        raise FileNotFoundError(
            f"No clinical data file found in {clinical_data_directory}.\n"
            "Please place the OASIS-2 demographics file (.csv or .xlsx) "
            "in that directory."
        )
    if len(all_files) > 1:
        raise ValueError(
            f"Multiple clinical data files found in {clinical_data_directory}: "
            f"{[f.name for f in all_files]}. "
            "Only one file is expected."
        )

    clinical_file = all_files[0]
    if clinical_file.suffix == ".xlsx":
        df = pd.read_excel(clinical_file)
    else:
        df = pd.read_csv(clinical_file)

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
    records = []
    for img_path in _find_imaging_data(imaging_data_directory):
        # img_path is relative to imaging_data_directory
        # e.g.  OAS2_0001_MR1/RAW/OAS2_0001_MR1_mpr-1_anon.img
        session_folder = img_path.parts[0]  # OAS2_0001_MR1
        tokens = session_folder.split("_")  # ['OAS2', '0001', 'MR1']
        if len(tokens) < 3:
            continue
        subject_id = f"{tokens[0]}_{tokens[1]}"  # OAS2_0001
        session_label = tokens[2]  # MR1 | MR2 | ...
        run_number = _identify_run(img_path)

        records.append(
            {
                "source_path": img_path,
                "Subject ID": subject_id,
                "session_label": session_label,
                "run_number": run_number,
            }
        )

    if not records:
        raise FileNotFoundError(
            f"No ANALYZE (.img) T1w acquisitions found under {imaging_data_directory}.\n"
            "Expected files matching the pattern: "
            "OAS2_XXXX_MRY/RAW/mpr-N.nifti.img"
        )

    df = pd.DataFrame(records)
    df["participant_id"] = df["Subject ID"].apply(
        lambda x: "sub-" + x.replace("_", "")  # OAS2_0001 → sub-OAS20001
    )
    df["session_id"] = df["session_label"].apply(lambda x: f"ses-{x}")
    return df.drop_duplicates().sort_values(by=["source_path"])


def _find_imaging_data(path_to_source_data: Path) -> Iterable[Path]:
    """Yield relative paths to every raw MPRAGE ANALYZE acquisition.

    Files follow the naming convention: mpr-N.nifti.img
    e.g. mpr-1.nifti.img, mpr-2.nifti.img, ...
    """
    for image in path_to_source_data.rglob("RAW/mpr-*.img"):
        yield image.relative_to(path_to_source_data)


def _identify_run(source_path: Path) -> str:
    """Return a BIDS run label from the mpr-N suffix in the filename.

    Examples
    --------
    ``mpr-1.nifti.img``  ->  ``run-01``
    ``mpr-4.nifti.img``  ->  ``run-04``
    """
    import re

    match = re.search(r"mpr-(\d+)", source_path.name)
    return f"run-{int(match.group(1)):02d}" if match else "run-01"


# ---------------------------------------------------------------------------
# Merging imaging and clinical data
# ---------------------------------------------------------------------------


def intersect_data(
    df_imaging: pd.DataFrame,
    df_clinical: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    df_subjects:
        One row per subject (baseline visit), used to build
        ``participants.tsv``.
    """
    df_clinical = df_clinical.copy()
    # Derive session label (MR1, MR2, ...) from the MRI ID column
    # MRI ID format: OAS2_XXXX_MRY
    df_clinical["session_label"] = df_clinical["MRI ID"].str.split("_").str[2]

    # Inner join: keep only subjects/sessions present in both data sources
    df_merged = df_imaging.merge(
        df_clinical,
        on=["Subject ID", "session_label"],
        how="inner",
    )

    # Build the target BIDS filename for each acquisition
    # Pattern: sub-OAS2XXXX/ses-MRY/anat/sub-OAS2XXXX_ses-MRY_run-0N_T1w.nii.gz
    df_merged = df_merged.assign(
        filename=lambda df: df.apply(
            lambda row: (
                f"{row.participant_id}/{row.session_id}/anat/"
                f"{row.participant_id}_{row.session_id}_{row.run_number}_T1w.nii.gz"
            ),
            axis=1,
        )
    )

    # Per-subject baseline record (earliest session) -> participants.tsv
    df_subjects = (
        df_clinical.sort_values("session_label")
        .drop_duplicates(subset=["Subject ID"], keep="first")[
            ["Subject ID", "M/F", "Hand", "EDUC", "SES"]
        ]
        .copy()
    )
    df_subjects["participant_id"] = df_subjects["Subject ID"].apply(
        lambda x: "sub-" + x.replace("_", "")
    )

    return df_merged, df_subjects


# ---------------------------------------------------------------------------
# Building BIDS metadata tables
# ---------------------------------------------------------------------------


def dataset_to_bids(
    df_merged: pd.DataFrame,
    df_subjects: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the three BIDS metadata tables from the merged data.

    Returns
    -------
    participants : pd.DataFrame  (index = participant_id)
    sessions     : pd.DataFrame  (index = [participant_id, session_id])
    scans        : pd.DataFrame  (index = BIDS filename)
    """
    return (
        _build_participants_df(df_subjects),
        _build_sessions_df(df_merged),
        _build_scans_df(df_merged),
    )


def _build_participants_df(df_subjects: pd.DataFrame) -> pd.DataFrame:
    """One row per subject: sex, handedness, education, socioeconomic status."""
    return (
        df_subjects[["participant_id", "M/F", "Hand", "EDUC", "SES"]]
        .rename(
            columns={
                "M/F": "sex",
                "Hand": "handedness",
                "EDUC": "education",
                "SES": "socioeconomic_status",
            }
        )
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
    return (
        df_merged[["filename", "source_path"]]
        .drop_duplicates(subset=["filename"])
        .set_index("filename", verify_integrity=True)
    )


# ---------------------------------------------------------------------------
# ANALYZE → NIfTI conversion
# ---------------------------------------------------------------------------


def _convert_analyze_to_nifti(source_path: Path, target_path: Path) -> None:
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
# Writing the BIDS directory
# ---------------------------------------------------------------------------


def write_bids(
    to: Path,
    participants: pd.DataFrame,
    sessions: pd.DataFrame,
    scans: pd.DataFrame,
    dataset_directory: Path,
) -> list[str]:
    """Write the full BIDS output directory.

    Parameters
    ----------
    to:
        Root of the target BIDS directory (created if absent).
    participants:
        Output of :func:`_build_participants_df`.
    sessions:
        Output of :func:`_build_sessions_df`.
    scans:
        Output of :func:`_build_scans_df`.
    dataset_directory:
        Root of the original OASIS-2 imaging directory.

    Returns
    -------
    list[str]
        List of BIDS filenames that were written.
    """
    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters._utils import write_to_tsv
    from clinica.dataset import BIDSDatasetDescription
    from clinica.utils.stream import cprint

    fs = LocalFileSystem(auto_mkdir=True)

    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as f:
            BIDSDatasetDescription(name="OASIS-2").write(to=f)
        with fs.open(to / "participants.tsv", "w") as f:
            write_to_tsv(participants, f)

    for participant_id, sessions_group in sessions.groupby("participant_id"):
        sessions_group = sessions_group.droplevel("participant_id")
        sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
        with fs.open(sessions_filepath, "w") as sf:
            write_to_tsv(sessions_group, sf)

    written = []
    for filename, metadata in scans.iterrows():
        source_full = dataset_directory / metadata.source_path
        target_full = to / filename

        if not source_full.exists():
            cprint(
                f"Source file not found, skipping: {source_full}",
                lvl="warning",
            )
            continue

        _convert_analyze_to_nifti(source_full, target_full)
        written.append(filename)

    return written
