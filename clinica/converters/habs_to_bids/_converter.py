from pathlib import Path
from typing import Iterator, Optional

import pandas as pd

from clinica.converters.study_models import StudyName, bids_id_factory
from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def _get_protocol_to_bids_df() -> pd.DataFrame:
    return pd.DataFrame.from_dict(
        {
            "3DFLAIR": {"datatype": "anat", "modality": "FLAIR"},
            "ADNI_1X-MPRAGE": {"datatype": "anat", "modality": "T1w"},
            "ADNI_2X-MPRAGE": {"datatype": "anat", "modality": "T1w"},
            "FDG": {
                "datatype": "pet",
                "modality": "pet",
                "trc_label": "18FFDG",
            },
            "FTP_75-105": {
                "datatype": "pet",
                "modality": "pet",
                "trc_label": "18FFTP",
            },
            "FTP_80-100": {
                "datatype": "pet",
                "modality": "pet",
                "trc_label": "18FFTP",
            },
            "PIB_40-60": {
                "datatype": "pet",
                "modality": "pet",
                "trc_label": "11CPIB",
            },
        },
        orient="index",
    )


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    *args,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    from clinica.converters.factory import get_converter_name
    from clinica.converters.study_models import StudyName
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.HABS)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.HABS)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )
    clinical_data = {
        k: _read_clinical_data(path_to_dataset / p, c)
        for k, p, c in _find_clinical_data(path_to_dataset)
    }
    imaging_data = pd.concat(
        [_parse_imaging_data(x) for x in _find_imaging_data(path_to_dataset)]
    )
    _write_bids(
        sourcedata=path_to_dataset,
        rawdata=bids_dir,
        imaging_data=imaging_data,
        clinical_data=clinical_data,
    )


def _source_session_id_to_bids(dataframe: pd.DataFrame) -> pd.Series:
    import re

    # HABS session format as `{years}.{months}`
    years_dot_months_pattern = r".*_(\d+).(\d+)"

    def replace_years_dot_months_session(match: re.Match) -> str:
        years = int(match.group(1)) - 1
        months = int(match.group(2))
        return f"ses-M{12 * years + months:03d}"

    # HABS undocumented format as `{months}m`
    months_only_pattern = r".*_(\d+)m.*"

    def replace_months_only_pattern(match: re.Match) -> str:
        months = int(match.group(1))
        return f"ses-M{months:03d}"

    return dataframe.source_session_id.str.replace(
        years_dot_months_pattern, replace_years_dot_months_session, regex=True
    ).str.replace(months_only_pattern, replace_months_only_pattern, regex=True)


def _find_clinical_data(sourcedata: Path) -> Iterator[tuple[str, Path, dict[str, str]]]:
    from csv import DictReader

    for f in sourcedata.rglob("*HABS_DataRelease*.csv"):
        header = DictReader(open(f, encoding="utf-8-sig")).fieldnames
        yield (
            f.name.split("_")[0],
            f.relative_to(sourcedata),
            {
                header[0]: "source_participant_id",
                header[1]: "source_session_id",
                header[2]: "date",
            },
        )


def _read_clinical_data(path: Path, rename_columns: dict[str, str]) -> pd.DataFrame:
    return (
        pd.read_csv(path)
        .rename(columns=rename_columns)
        .assign(date=lambda df: pd.to_datetime(df.date))
        .assign(
            participant_id=lambda df: df.source_participant_id.apply(
                lambda x: bids_id_factory(StudyName.HABS).from_original_study_id(x)
            )
        )
        .drop(columns="source_participant_id")
        .assign(session_id=_source_session_id_to_bids)
        .drop(columns="source_session_id")
        .convert_dtypes()
        .set_index(["participant_id", "session_id", "date"], verify_integrity=True)
        .sort_index()
    )


def _find_imaging_data(sourcedata: Path) -> Iterator[list[tuple[str, str]]]:
    from zipfile import ZipFile

    for z in sourcedata.rglob("*HABS_DataRelease*.zip"):
        yield [
            (str(z.relative_to(sourcedata)), f)
            for f in ZipFile(z).namelist()
            if f.endswith(".nii.gz")
        ]


def _parse_imaging_data(paths: list[tuple[str, str]]) -> Optional[pd.DataFrame]:
    df = pd.DataFrame.from_records(paths, columns=["source_zipfile", "source_filename"])
    # Parse imaging metadata from file paths.
    pattern = (
        r"(?P<source_participant_id>P_\w{6})"
        r"_(?P<date>[\d-]+)"
        r"_(?P<source_session_id>HAB_\d{1}.\d{1,2})"
        r"_(?P<protocol>[\w-]+)"
    )
    df = pd.concat(
        [df, df["source_filename"].str.extract(pattern)],
        axis="columns",
    ).assign(date=lambda x: pd.to_datetime(x.date))

    # Map protocol to BIDS entities.
    df = (
        df.join(_get_protocol_to_bids_df(), on="protocol", how="left")
        .drop(columns="protocol")
        .dropna(subset=["datatype", "modality"])
    )
    if df.empty:
        return None
    # Compute BIDS participant ID, session ID and filename.
    df = (
        df.assign(
            participant_id=lambda df: df.source_participant_id.apply(
                lambda x: bids_id_factory(StudyName.HABS).from_original_study_id(x)
            )
        )
        .drop(columns="source_participant_id")
        .assign(session_id=_source_session_id_to_bids)
        .drop(columns="source_session_id")
        .set_index(
            ["participant_id", "session_id", "datatype", "modality", "trc_label"],
            verify_integrity=True,
        )
        .sort_index()
    )

    return df


def _write_bids(
    sourcedata: Path,
    rawdata: Path,
    imaging_data: pd.DataFrame,
    clinical_data: dict[str, pd.DataFrame],
):
    import fsspec
    from pandas import notna

    from clinica.converters._utils import (
        install_nifti,
        write_modality_agnostic_files,
        write_to_tsv,
    )
    from clinica.converters.study_models import StudyName
    from clinica.dataset import BIDSDatasetDescription

    participants = (
        clinical_data["Demographics"]
        .drop(columns="NP_Age")
        .rename(
            columns={
                "BiologicalSex": "sex",
                "YrsOfEd": "years_of_education",
            }
        )
        .xs("ses-M000", level="session_id")
        .reset_index(level="date", drop=True)
    ).rename(columns=str.lower)

    with fsspec.open(str(rawdata / "dataset_description.json"), mode="wt") as f:
        BIDSDatasetDescription(name=StudyName.HABS).write(to=f)

    participants_file = rawdata / "participants.tsv"
    with fsspec.open(str(participants_file), mode="wb") as f:
        write_to_tsv(participants, f)

    sessions = (
        clinical_data["Demographics"][["NP_Age"]]
        .reset_index(level="date", drop=True)
        .rename(columns={"NP_Age": "age"})
        .join(
            clinical_data["ClinicalMeasures"]
            .rename(columns={"HABS_DX": "diagnostic"})
            .reset_index(level="date", drop=True)
        )
    ).rename(columns=str.lower)

    for participant_id, dataframe in sessions.groupby("participant_id"):
        sessions_file = rawdata / participant_id / f"{participant_id}_sessions.tsv"
        dataframe = dataframe.droplevel("participant_id")
        with fsspec.open(sessions_file, mode="wb") as f:
            write_to_tsv(dataframe, f)

    for grouped_by, dataframe in imaging_data.groupby(["participant_id", "session_id"]):
        participant_id, session_id = grouped_by
        bids_basedir = rawdata / participant_id / session_id
        bids_prefix = f"{participant_id}_{session_id}"
        dataframe = (
            dataframe.droplevel(level=["participant_id", "session_id"])
            .reset_index()
            .assign(
                filename=lambda df: df.apply(
                    lambda row: (
                        f"{row.datatype}/{bids_prefix}"
                        f"{'_trc-' + row.trc_label if notna(row.trc_label) else ''}"
                        f"_{row.modality}.nii.gz"
                    ),
                    axis="columns",
                )
            )
            .drop(columns=["datatype", "modality", "trc_label"])
            .set_index("filename", verify_integrity=True)
        )

        for filename, row in dataframe.iterrows():
            install_nifti(
                sourcedata_dir=sourcedata / row.source_zipfile,
                bids_filename=bids_basedir / filename,
                source_filename=row.source_filename,
            )

        dataframe = dataframe.drop(columns=["source_zipfile", "source_filename"])

        with fsspec.open(bids_basedir / f"{bids_prefix}_scans.tsv", mode="wb") as f:
            write_to_tsv(dataframe, f)

    readme_data = {
        "link": "https://habs.mgh.harvard.edu",
        "desc": (
            "The overall goal of the Harvard Aging Brain Study (HABS) is to elucidate the earliest changes in "
            "molecular, functional and structural imaging markers that signal the transition from normal cognition to "
            "progressive cognitive decline along the trajectory of preclinical Alzheimer’s Disease."
        ),
    }
    write_modality_agnostic_files(
        study_name=StudyName.HABS,
        readme_data=readme_data,
        bids_dir=rawdata,
    )
