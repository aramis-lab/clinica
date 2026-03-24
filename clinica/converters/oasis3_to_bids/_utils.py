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


def read_clinical_data(clinical_data_directory: Path) -> dict[str, pd.DataFrame]:
    """Read clinical data from the OASIS3_data_files FTP directory structure."""
    cprint(f"Reading clinical data from {clinical_data_directory}", lvl="info")
    csv_files = list(clinical_data_directory.rglob("*.csv"))
    if not csv_files:
        cprint(f"No CSV files found under {clinical_data_directory}.", lvl="error")
    file_map: dict[str, pd.DataFrame] = {}
    for f in csv_files:
        cprint(f"  Loading {f.name}", lvl="debug")
        file_map[f.stem] = pd.read_csv(f)

    dict_df: dict[str, pd.DataFrame] = {}
    for category, filenames in _CLINICAL_FILES.items():
        if category == "pet":
            # Concatenate all PET files (different tracers → add rows).
            dfs = [file_map[name] for name in filenames if name in file_map]
            if dfs:
                dict_df["pet"] = pd.concat(dfs, ignore_index=True)
        elif category == "clinical":
            # Merge UDS sub-files side-by-side on shared visit keys.
            # Normalise days_to_visit to int first: UDSb1 stores it as a
            # zero-padded string (e.g. "0339") while UDSb4 stores it as int.
            dfs = []
            for name in filenames:
                if name in file_map:
                    df = file_map[name].copy()
                    df["days_to_visit"] = (
                        pd.to_numeric(df["days_to_visit"], errors="coerce")
                        .fillna(0)
                        .astype(int)
                    )
                    dfs.append(df)
            if len(dfs) >= 2:
                merged = dfs[0]
                for df in dfs[1:]:
                    merged = merged.merge(df, on=_CLINICAL_MERGE_KEYS, how="outer")
                dict_df["clinical"] = merged
            elif len(dfs) == 1:
                dict_df["clinical"] = dfs[0]
        else:
            # "mri" and "demo": single file per category.
            for name in filenames:
                if name in file_map:
                    dict_df[category] = file_map[name]
                    break

    return dict_df


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


def _find_imaging_data(path_to_source_data: Path) -> Iterable[Path]:
    for image in path_to_source_data.rglob("*.nii.gz"):
        yield image.relative_to(path_to_source_data)


def intersect_data(
    df_source: pd.DataFrame, dict_df: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Rename new-schema columns to the names expected by downstream functions.
    df_clinical = dict_df["clinical"].rename(
        columns={
            "OASISID": "Subject",
            "MMSE": "mmse",
            "CDRTOT": "cdr",
            "CDRSUM": "sumbox",
        }
    )
    df_demo = dict_df["demo"].rename(
        columns={
            "OASISID": "Subject",
            "AgeatEntry": "ageAtEntry",
            "GENDER": "M/F",
            "HAND": "Hand",
            "EDUC": "Education",
            "APOE": "apoe",
        }
    )

    # Derive a d0000-style session_id from days_to_visit so we can identify the
    # baseline visit and compute session bins consistently with the imaging data.
    df_clinical = df_clinical.assign(
        session_id=lambda df: pd.to_numeric(df["days_to_visit"], errors="coerce")
        .fillna(0)
        .astype(int)
        .apply(lambda x: f"d{x:04d}")
    )

    # Baseline visit (d0000) per subject.
    df_clinical_small = df_clinical[df_clinical.session_id == "d0000"]
    df_clinical_small = df_clinical_small.merge(
        df_source["Subject"], how="inner", on="Subject"
    ).drop_duplicates()

    # Age at entry and APOE come from demographics (one row per subject).
    df_source = df_source.merge(
        df_demo[["Subject", "ageAtEntry", "apoe"]], how="inner", on="Subject"
    )
    df_source = df_source.assign(
        age=lambda df: df["ageAtEntry"] + df["Date"].str[1:].astype("float") / 365.25
    )

    # Subject-level demographics filtered to subjects present in the imaging data.
    df_subject_small = df_demo.merge(
        df_source["Subject"], how="inner", on="Subject"
    ).drop_duplicates()

    # Add MRI acquisition metadata (scanner manufacturer) to imaging sessions.
    mri_scanner = dict_df["mri"][
        ["label", "Manufacturer", "ManufacturersModelName"]
    ].drop_duplicates(subset=["label"])
    df_source = df_source.merge(
        mri_scanner, how="left", left_on="path", right_on="label"
    )

    # Merge demographics baseline + clinical baseline for the participants file.
    df_small = df_subject_small.merge(df_clinical_small, how="inner", on="Subject")
    df_small = df_small.assign(participant_id=lambda df: "sub-" + df.Subject)

    # Convert imaging session date (dXXXX days) to a 6-month BIDS session bin.
    df_source = df_source.assign(
        session=lambda df: (round(df["Date"].str[1:].astype("int") / (365.25 / 2)) * 6)
    )
    df_source = df_source.assign(
        ses=lambda df: df.session.apply(lambda x: f"ses-M{str(int(x)).zfill(3)}")
    )
    df_source = df_source.join(
        df_source.modality.map(
            {
                "dwi_MR": {"datatype": "dwi", "suffix": "dwi"},
                "T1w_MR": {"datatype": "anat", "suffix": "T1w"},
                "T2star_MR": {"datatype": "anat", "suffix": "T2starw"},
                "FLAIR_MR": {"datatype": "anat", "suffix": "FLAIR"},
                "pet_FDG": {"datatype": "pet", "suffix": "pet", "trc_label": "18FFDG"},
                "pet_PIB": {"datatype": "pet", "suffix": "pet", "trc_label": "11CPIB"},
                "pet_AV45": {
                    "datatype": "pet",
                    "suffix": "pet",
                    "trc_label": "18FAV45",
                },
                "pet_AV1451": {
                    "datatype": "pet",
                    "suffix": "pet",
                    "trc_label": "18FAV1451",
                },
            }
        ).apply(pd.Series)
    )
    if "trc_label" in df_source.columns:
        df_source = df_source.assign(
            filename=lambda df: df.apply(
                lambda x: f"{x.participant_id}/{x.ses}/{x.datatype}/"
                f"{x.participant_id}_{x.ses}"
                f"{'_trc-' + x.trc_label if pd.notna(x.trc_label) else ''}"
                f"_{x.run_number}_{x.suffix}.nii.gz",
                axis=1,
            )
        )
    else:
        df_source = df_source.assign(
            filename=lambda df: df.apply(
                lambda x: f"{x.participant_id}/{x.ses}/{x.datatype}/"
                f"{x.participant_id}_{x.ses}"
                f"_{x.run_number}_{x.suffix}.nii.gz",
                axis=1,
            )
        )

    # Convert clinical visit days to session bins and merge onto imaging sessions.
    df_clinical = df_clinical.merge(df_source["Subject"], how="inner", on="Subject")
    df_clinical = df_clinical.assign(
        session=lambda df: round(
            pd.to_numeric(df["days_to_visit"], errors="coerce").fillna(0) / (364.25 / 2)
        )
        * 6
    )
    df_clinical = df_clinical.drop_duplicates().set_index(["Subject", "session_id"])
    df_source = df_source.merge(
        df_clinical[
            [
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
        ],
        how="left",
        on="session",
    )
    return df_source, df_small


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

    n_scans = len(scans)
    cprint(f"Installing {n_scans} scan(s) into BIDS dataset at {to}", lvl="info")
    for i, (filename, metadata) in enumerate(scans.iterrows(), start=1):
        path = dataset_directory / metadata.source_dir
        suffix = _extract_suffix_from_filename(str(filename))
        if suffix != "nan":
            # Extract subject and session from the BIDS filename for user-facing log.
            parts = Path(filename).parts  # e.g. sub-OAS30001/ses-M000/anat/...
            subject = parts[0] if len(parts) > 0 else "unknown"
            session = parts[1] if len(parts) > 1 else "unknown"
            cprint(
                f"[{i}/{n_scans}] Installing {subject} / {session} : {Path(filename).name}",
                lvl="info",
            )
            _install_bids(sourcedata_dir=path, bids_filename=to / filename)
    return scans.index.to_list()


def _extract_suffix_from_filename(filename: str) -> str:
    return filename.split("_")[-1].split(".")[0]
