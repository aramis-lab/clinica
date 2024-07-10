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


def read_clinical_data(clinical_data_directory: Path) -> dict[str, pd.DataFrame]:
    try:
        image_data_files = [f for f in clinical_data_directory.glob("*.csv")]
    except StopIteration:
        raise FileNotFoundError("Imaging collection file not found.")
    image_metadata = [pd.read_csv(image_path) for image_path in image_data_files]
    for i, df in enumerate(image_metadata):
        image_metadata[i] = df.set_index(df.axes[1][0])

    return {
        k: _get_df_based_on_index_name(image_metadata, index_name)
        for k, index_name in zip(
            ["pet", "mri", "subject", "pup", "adrc"],
            [
                "XNAT_PETSESSIONDATA ID",
                "MR ID",
                "Subject",
                "PUP_PUPTIMECOURSEDATA ID",
                "ADRC_ADRCCLINICALDATA ID",
            ],
        )
    }


def _get_df_based_on_index_name(
    image_metadata: list[pd.DataFrame], index_name: str
) -> pd.DataFrame:
    matching_dfs = [df for df in image_metadata if df.index.name == index_name]
    if len(matching_dfs) == 0:
        raise FileNotFoundError(f"Clinical data not found for {index_name}.")
    if len(matching_dfs) > 1:
        raise ValueError(f"Multiple data found for {index_name}.")
    return matching_dfs[0]


def read_imaging_data(imaging_data_directory: Path) -> pd.DataFrame:
    from clinica.iotools.bids_utils import StudyName, bids_id_factory

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
    df_adrc = dict_df["adrc"].assign(
        session_id=lambda df: df.index.map(lambda x: x.split("_")[2])
    )
    df_adrc_small = df_adrc[df_adrc.session_id == "d0000"]
    df_adrc_small = df_adrc_small.merge(
        df_source["Subject"], how="inner", on="Subject"
    ).drop_duplicates()
    df_source = df_source.merge(
        df_adrc_small[["Subject", "ageAtEntry"]], how="inner", on="Subject"
    )
    df_source = df_source.assign(
        age=lambda df: (df["ageAtEntry"]) + df["Date"].str[1:].astype("float") / 365.25
    )
    df_subject_small = (
        dict_df["subject"]
        .merge(df_source["Subject"], how="inner", on="Subject")
        .drop_duplicates()
    )
    df_source = df_source.merge(
        dict_df["mri"]["Scanner"], how="left", left_on="path", right_on="MR ID"
    )
    df_small = df_subject_small.merge(df_adrc_small, how="inner", on="Subject")
    df_small = df_small.assign(participant_id=lambda df: "sub-" + df.Subject)
    df_source = df_source.assign(
        session=lambda df: (round(df["Date"].str[1:].astype("int") / (365.25 / 2)) * 6)
    )
    df_source = df_source.assign(
        ses=lambda df: df.session.apply(
            lambda x: f"ses-M{(5 - (len(str(x)))) * '0' + str(int(x))}"
        )
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
    df_adrc = df_adrc.merge(df_source["Subject"], how="inner", on="Subject")
    df_adrc = df_adrc.assign(
        session=lambda df: round(df["session_id"].str[1:].astype("int") / (364.25 / 2))
        * 6
    )
    df_adrc = df_adrc.drop_duplicates().set_index(["Subject", "session_id"])
    df_source = df_source.merge(
        df_adrc[
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
                "apoe",
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
        df_small[["participant_id", "ageAtEntry", "M/F", "Hand", "Education"]]
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
                "apoe",
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

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription
    from clinica.iotools.bids_utils import StudyName, write_to_tsv

    fs = LocalFileSystem(auto_mkdir=True)

    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name=StudyName.OASIS3).write(
                to=dataset_description_file
            )
        with fs.open(to / "participants.tsv", "w") as participant_file:
            write_to_tsv(participants, participant_file)
        for participant_id, sessions_group in sessions.groupby("participant_id"):
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(sessions_filepath, "w") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    for filename, metadata in scans.iterrows():
        path = dataset_directory / metadata.source_dir
        if _extract_suffix_from_filename(str(filename)) != "nan":
            _install_bids(sourcedata_dir=path, bids_filename=to / filename)
    return scans.index.to_list()


def _extract_suffix_from_filename(filename: str) -> str:
    return filename.split("_")[-1].split(".")[0]
