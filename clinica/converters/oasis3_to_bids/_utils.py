from pathlib import Path
from typing import Iterable

import pandas as pd

from clinica.utils.stream import cprint

__all__ = [
    "read_clinical_data",
    "read_imaging_data",
    "intersect_data",
    "dataset_to_bids",
    "write_bids",
]

# Hardcode relevant .csv filenames from the standardized OASIS3_data_files directory.
# Each key maps to the stem(s) of the expected CSV file(s) within that subdirectory.
_CLINICAL_FILES = {
    "pet": ["OASIS3_PET_json", "OASIS3_AV1451_PET_json", "OASIS3_AV1451L_json"],
    "mri": ["OASIS3_MR_json"],
    "clinical": ["OASIS3_UDSb1_physical_eval", "OASIS3_UDSb4_cdr"],
    "demo": ["OASIS3_demographics"],
}

# Columns shared across the UDS clinical sub-files used to join visits.
_CLINICAL_MERGE_KEYS = ["OASISID", "days_to_visit", "age at visit"]

# Mapping from OASIS3 modality strings to BIDS datatype / suffix / tracer label.
_MODALITY_TO_BIDS = {
    "dwi_MR": {"datatype": "dwi", "suffix": "dwi"},
    "T1w_MR": {"datatype": "anat", "suffix": "T1w"},
    "T2star_MR": {"datatype": "anat", "suffix": "T2starw"},
    "FLAIR_MR": {"datatype": "anat", "suffix": "FLAIR"},
    "bold_MR": {"datatype": "func", "suffix": "bold"},
    "pet_FDG": {"datatype": "pet", "suffix": "pet", "trc_label": "18FFDG"},
    "pet_PIB": {"datatype": "pet", "suffix": "pet", "trc_label": "11CPIB"},
    "pet_AV45": {"datatype": "pet", "suffix": "pet", "trc_label": "18FAV45"},
    "pet_AV1451": {"datatype": "pet", "suffix": "pet", "trc_label": "18FAV1451"},
}

# Clinical Score column names
_CLINICAL_SCORE_COLUMNS = [
    "session",
    "mmse",
    "cdr",
    "commun",
    "dx1",
    "homehobb",
    "judgment",
    "memory",
    "orient",
    "perscare",
    "sumbox",
]


def _build_file_map(data_directory: Path) -> dict[str, pd.DataFrame]:
    """Recursively find all CSVs under data_directory and return {stem: DataFrame}."""
    csv_files = list(data_directory.rglob("*.csv"))
    if not csv_files:
        cprint(f"No CSV files found under {data_directory}.", lvl="error")
    file_map = {}
    for f in csv_files:
        cprint(f"  Loading {f.name}", lvl="debug")
        file_map[f.stem] = pd.read_csv(f)
    return file_map


def _read_clinical_data(
    file_map: dict[str, pd.DataFrame], filenames: list[str]
) -> pd.DataFrame:
    dfs = [file_map[name].copy() for name in filenames if name in file_map]
    if not dfs:
        raise FileNotFoundError(f"None of {filenames} found in data directory.")

    df_clinical = _merge_clinical_df(dfs) if len(dfs) > 1 else dfs[0]
    # Rename new-schema columns to the names expected by downstream functions.
    df_clinical = df_clinical.rename(
        columns={
            "OASISID": "Subject",
            "MMSE": "mmse",
            "CDRTOT": "cdr",
            "CDRSUM": "sumbox",
        }
    )
    # Derive a d0000-style session_id from days_to_visit so we can identify the
    # baseline visit and compute session bins consistently with the imaging data.
    df_clinical = df_clinical.assign(
        session_id=lambda df: "d" + df["days_to_visit"].astype(str).str.zfill(4)
    )
    return df_clinical


def _read_pet_data(
    file_map: dict[str, pd.DataFrame], filenames: list[str]
) -> pd.DataFrame:
    dfs = [file_map[name] for name in filenames if name in file_map]
    if not dfs:
        raise FileNotFoundError(f"None of {filenames} found in data directory.")
    return pd.concat(dfs, ignore_index=True)


def _read_demo_data(file_map: dict[str, pd.DataFrame], filename: str) -> pd.DataFrame:
    if filename not in file_map:
        raise FileNotFoundError(f"{filename} not found in data directory.")
    return (
        file_map[filename]
        .copy()
        .rename(
            columns={
                "OASISID": "Subject",
                "AgeatEntry": "ageAtEntry",
                "GENDER": "M/F",
                "HAND": "Hand",
                "EDUC": "Education",
                "APOE": "apoe",
            }
        )
    )


def _read_mri_data(file_map: dict[str, pd.DataFrame], filename: str) -> pd.DataFrame:
    if filename not in file_map:
        raise FileNotFoundError(f"{filename} not found in data directory.")
    return file_map[filename]


def _merge_clinical_df(list_dfs: list) -> pd.DataFrame:
    merged_df = list_dfs[0]
    for df in list_dfs[1:]:
        merged_df = merged_df.merge(df, on=_CLINICAL_MERGE_KEYS, how="outer")
    return merged_df


def read_clinical_data(clinical_data_directory: Path) -> dict[str, pd.DataFrame]:
    """Read clinical data from the OASIS3_data_files FTP directory structure."""
    cprint(f"Reading clinical data from {clinical_data_directory}", lvl="info")
    file_map = _build_file_map(clinical_data_directory)
    return {
        "pet": _read_pet_data(file_map, _CLINICAL_FILES["pet"]),
        "clinical": _read_clinical_data(file_map, _CLINICAL_FILES["clinical"]),
        "mri": _read_mri_data(file_map, _CLINICAL_FILES["mri"][0]),
        "demo": _read_demo_data(file_map, _CLINICAL_FILES["demo"][0]),
    }


def _find_imaging_data(path_to_source_data: Path) -> Iterable[Path]:
    for image in path_to_source_data.rglob("*.nii.gz"):
        yield image.relative_to(path_to_source_data)


def _identify_modality(source_path: Path) -> str:
    try:
        return source_path.name.split(".")[0].split("_")[-1]
    except Exception:
        return "nan"


def _identify_runs(source_path: Path) -> str:
    import re

    try:
        return re.search(r"run-\d+", str(source_path))[0]
    except Exception:
        return "run-01"


def read_imaging_data(imaging_data_directory: Path) -> pd.DataFrame:
    from clinica.converters.study_models import StudyName, bids_id_factory

    source_path_series = pd.Series(
        _find_imaging_data(imaging_data_directory), name="source_path"
    )
    source_dir_series = source_path_series.apply(lambda x: x.parent).rename(
        "source_dir"
    )
    file_spec_series = source_path_series.apply(lambda x: x.parts[0]).rename("path")
    source_file_series = source_path_series.apply(
        lambda x: _identify_modality(x)
    ).rename("modality")
    source_run_series = source_path_series.apply(lambda x: _identify_runs(x)).rename(
        "run_number"
    )
    df_source = pd.concat(
        [
            source_path_series,
            file_spec_series,
            source_dir_series,
            file_spec_series.str.split("_", expand=True),
            source_file_series,
            source_run_series,
        ],
        axis=1,
    )
    df_source = (
        df_source.rename({0: "Subject", 1: "modality_2", 2: "Date"}, axis="columns")
        .drop_duplicates()
        .sort_values(by=["source_path"])
    )
    df_source = df_source.assign(
        participant_id=lambda df: df.Subject.apply(
            lambda x: bids_id_factory(StudyName.OASIS3).from_original_study_id(x)
        )
    )
    df_source["modality"] = df_source[["modality", "modality_2"]].apply(
        "_".join, axis=1
    )
    return df_source


def _filter_to_imaging_subjects(
    df: pd.DataFrame, imaging_subjects: pd.Series
) -> pd.DataFrame:
    """Keep only rows for subjects present in the imaging data."""
    return df.merge(imaging_subjects, how="inner", on="Subject").drop_duplicates()


def _get_baseline_clinical(
    df_clinical: pd.DataFrame, imaging_subjects: pd.Series
) -> pd.DataFrame:
    """Filter clinical data to the baseline visit (d0000) for imaging subjects only."""
    return _filter_to_imaging_subjects(
        df_clinical[df_clinical.session_id == "d0000"], imaging_subjects
    )


def _add_demo_and_age_at_scan(
    df_source: pd.DataFrame, df_demo: pd.DataFrame
) -> pd.DataFrame:
    """Merge demographics into imaging data and compute age at each scan."""
    return df_source.merge(
        df_demo[["Subject", "ageAtEntry", "apoe"]], how="inner", on="Subject"
    ).assign(
        age=lambda df: df["ageAtEntry"] + df["Date"].str[1:].astype("float") / 365.25
    )


def _add_mri_scanner_metadata(
    df_source: pd.DataFrame, df_mri: pd.DataFrame
) -> pd.DataFrame:
    """Left-join MRI scanner metadata (manufacturer, model, field strength) onto imaging sessions."""
    mri_scanner = df_mri[
        ["label", "Manufacturer", "ManufacturersModelName", "MagneticFieldStrength"]
    ].drop_duplicates(subset=["label"])
    return df_source.merge(mri_scanner, how="left", left_on="path", right_on="label")


def _merge_baseline_data(
    df_subject_small: pd.DataFrame, df_clinical_small: pd.DataFrame
) -> pd.DataFrame:
    """Combine subject-level demographics with baseline clinical scores."""
    return df_subject_small.merge(df_clinical_small, how="inner", on="Subject").assign(
        participant_id=lambda df: "sub-" + df.Subject
    )


def _days_to_session_bin(days: pd.Series) -> pd.Series:
    """Convert a Series of days-from-entry to 6-month BIDS session bin integers."""
    return round(days / (365.25 / 2)) * 6


def _compute_session_bins(df_source: pd.DataFrame) -> pd.DataFrame:
    """Convert image acquisition days to 6-month BIDS session bins (ses-MXXX)."""
    return df_source.assign(
        session=lambda df: _days_to_session_bin(df["Date"].str[1:].astype("int"))
    ).assign(ses=lambda df: df.session.apply(lambda x: f"ses-M{int(x):03d}"))


def _map_modality_to_bids(df_source: pd.DataFrame) -> pd.DataFrame:
    """Expand the modality string into BIDS datatype, suffix, and tracer label columns."""
    return df_source.join(df_source.modality.map(_MODALITY_TO_BIDS).apply(pd.Series))


def _build_bids_filenames(df_source: pd.DataFrame) -> pd.DataFrame:
    """Construct the BIDS-relative filename for each scan."""

    def _make_filename(x) -> str:
        trc = (
            f"_trc-{x.trc_label}"
            if "trc_label" in x.index and pd.notna(x.trc_label)
            else ""
        )
        return (
            f"{x.participant_id}/{x.ses}/{x.datatype}/"
            f"{x.participant_id}_{x.ses}{trc}_{x.run_number}_{x.suffix}.nii.gz"
        )

    return df_source.assign(filename=lambda df: df.apply(_make_filename, axis=1))


def _merge_clinical_scores(
    df_source: pd.DataFrame, df_clinical: pd.DataFrame
) -> pd.DataFrame:
    """Bin clinical visits into session bins and merge scores onto imaging sessions."""
    df_clinical = (
        df_clinical.merge(df_source["Subject"], how="inner", on="Subject")
        .assign(session=lambda df: _days_to_session_bin(df["days_to_visit"]))
        .drop_duplicates()
        .set_index(["Subject", "session_id"])
    )
    return df_source.merge(
        df_clinical[_CLINICAL_SCORE_COLUMNS], how="left", on="session"
    )


def intersect_data(
    df_source: pd.DataFrame, dict_df: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_clinical = dict_df["clinical"]
    df_demo = dict_df["demo"]

    df_clinical_small = _get_baseline_clinical(df_clinical, df_source["Subject"])
    df_source = _add_demo_and_age_at_scan(df_source, df_demo)
    df_subject_small = _filter_to_imaging_subjects(df_demo, df_source["Subject"])
    df_source = _add_mri_scanner_metadata(df_source, dict_df["mri"])
    df_baseline = _merge_baseline_data(df_subject_small, df_clinical_small)
    df_source = _compute_session_bins(df_source)
    df_source = _map_modality_to_bids(df_source)
    df_source = _build_bids_filenames(df_source)
    df_source = _merge_clinical_scores(df_source, df_clinical)
    return df_source, df_baseline


def dataset_to_bids(
    df_source: pd.DataFrame, df_small: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return (
        _build_participants_df(df_small),
        _build_sessions_df(df_source),
        _build_scans_df(df_source),
    )


def _build_participants_df(df_small: pd.DataFrame) -> pd.DataFrame:
    df_participants = (
        df_small[["participant_id", "ageAtEntry", "M/F", "Hand", "Education", "apoe"]]
        .rename(
            columns={
                "ageAtEntry": "age",
                "M/F": "sex",
                "Hand": "handedness",
            }
        )
        .set_index("participant_id", verify_integrity=True)
    )
    return df_participants


def _build_sessions_df(df_source: pd.DataFrame) -> pd.DataFrame:
    df_session = (
        df_source[
            [
                "participant_id",
                "ses",
                "Date",
                "age",
                "mmse",
                "cdr",
                "commun",
                "dx1",
                "homehobb",
                "judgment",
                "memory",
                "orient",
                "perscare",
                "sumbox",
            ]
        ]
        .rename(
            columns={
                "ses": "session_id",
                "Date": "source_session_id",
            }
        )
        .astype({"age": "int"})
        .drop_duplicates(subset=["participant_id", "session_id"])
        .set_index(["participant_id", "session_id"], verify_integrity=True)
    )
    return df_session


def _build_scans_df(df_source: pd.DataFrame) -> pd.DataFrame:
    df_scan = (
        df_source[["filename", "source_dir"]]
        .drop_duplicates(subset=["filename"])
        .set_index("filename", verify_integrity=True)
    )
    return df_scan


def _install_bids(sourcedata_dir: Path, bids_filename: Path) -> None:
    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)

    source_file = fs.open(fs.ls(sourcedata_dir)[0], mode="rb")
    cprint(
        f"  Copying {fs.ls(sourcedata_dir)[0]} -> {bids_filename}",
        lvl="debug",
    )
    target_file = fs.open(bids_filename, mode="wb")

    with source_file as sf, target_file as tf:
        tf.write(sf.read())

    source_basename = Path(Path(Path(fs.ls(sourcedata_dir)[0]).stem).stem)
    target_basename = Path(bids_filename.stem).stem

    # The following part adds the sidecar files related to the nifti with the same name: it can be tsv or json files.
    # It may or may not be used, since there might not be any sidecars.
    sidecar_dir = sourcedata_dir.parent / "BIDS"
    for source_sidecar in sidecar_dir.rglob(f"{source_basename}*"):
        target_sidecar = (bids_filename.parent / target_basename).with_name(
            f"{target_basename}{source_sidecar.suffix}"
        )
        source_file = fs.open(source_sidecar, mode="rb")
        target_file = fs.open(target_sidecar, mode="wb")
        with source_file as sf, target_file as tf:
            tf.write(sf.read())


def write_bids(
    to: Path,
    participants: pd.DataFrame,
    sessions: pd.DataFrame,
    scans: pd.DataFrame,
    dataset_directory: Path,
) -> list[str]:
    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters._utils import write_to_tsv
    from clinica.converters.study_models import StudyName
    from clinica.dataset import BIDSDatasetDescription

    fs = LocalFileSystem(auto_mkdir=True)

    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name=StudyName.OASIS3).write(
                to=dataset_description_file
            )
        with fs.open(to / "participants.tsv", "w") as participant_file:
            write_to_tsv(participants, participant_file)
        for participant_id, sessions_group in sessions.groupby("participant_id"):
            cprint(f"Writing sessions for {participant_id}", lvl="info")
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(sessions_filepath, "w") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    cprint(f"Installing {len(scans)} scan(s) into BIDS dataset at {to}", lvl="info")
    for i, (filename, metadata) in enumerate(scans.iterrows()):
        path = dataset_directory / metadata.source_dir
        if _extract_suffix_from_filename(str(filename)) != "nan":
            # Extract subject and session from the BIDS filename for user-facing log.
            parts = Path(filename).parts  # e.g. sub-OAS30001/ses-M000/anat/...
            subject = parts[0] if len(parts) > 0 else "unknown"
            session = parts[1] if len(parts) > 1 else "unknown"
            cprint(
                f"[{i+1}/{len(scans)}] Installing {subject} / {session} :"
                f" {Path(filename).name}",
                lvl="info",
            )
            _install_bids(sourcedata_dir=path, bids_filename=to / filename)
    return scans.index.to_list()


def _extract_suffix_from_filename(filename: str) -> str:
    return filename.split("_")[-1].split(".")[0]
