from io import TextIOBase
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from pandas import DataFrame, Series

import clinica.iotools.bids_utils as bids

PROTOCOL_TO_BIDS = {
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
}


def source_participant_id_to_bids(dataframe: DataFrame) -> Series:
    # HABS participant format prefixed with `P_`
    return dataframe.source_participant_id.str.replace("P_", "sub-HABS", regex=False)


def source_session_id_to_bids(dataframe: DataFrame) -> Series:
    import re

    # HABS session format as `{years}.{months}`
    years_dot_months_pattern = r".*_(\d+).(\d+)"

    def replace_years_dot_months_session(match: re.Match) -> str:
        years = int(match.group(1)) - 1
        months = int(match.group(2))
        return f"ses-M{12 * years + months:02d}"

    # HABS undocumented format as `{months}m`
    months_only_pattern = r".*_(\d+)m.*"

    def replace_months_only_pattern(match: re.Match) -> str:
        months = int(match.group(1))
        return f"ses-M{months:02d}"

    return dataframe.source_session_id.str.replace(
        years_dot_months_pattern, replace_years_dot_months_session, regex=True
    ).str.replace(months_only_pattern, replace_months_only_pattern, regex=True)


def source_filename_to_bids(dataframe: DataFrame) -> Series:
    from pandas import notna

    cols = ["participant_id", "session_id", "datatype", "modality", "trc_label"]
    return dataframe[cols].apply(
        lambda row: (
            f"{row.datatype}/{row.participant_id}_{row.session_id}"
            f"{'_trc-'+row.trc_label if notna(row.trc_label) else ''}"
            f"_{row.modality}.nii.gz"
        ),
        axis=1,
    )


def find_clinical_data(sourcedata: Path):
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


def read_clinical_data(path: Path, rename_columns: Dict[str, str]) -> DataFrame:
    from pandas import read_csv

    return (
        read_csv(path)
        .rename(columns=rename_columns)
        .assign(date=lambda df: pd.to_datetime(df.date))
        .assign(participant_id=source_participant_id_to_bids)
        .drop(columns="source_participant_id")
        .assign(session_id=source_session_id_to_bids)
        .drop(columns="source_session_id")
        .convert_dtypes()
        .set_index(["participant_id", "session_id", "date"], verify_integrity=True)
        .sort_index()
    )


def find_imaging_data(sourcedata: Path) -> List[Tuple[str, str]]:
    from zipfile import ZipFile

    ext = ".nii.gz"

    for z in sourcedata.rglob("*HABS_DataRelease*.zip"):
        yield [
            (str(z.relative_to(sourcedata)), f)
            for f in ZipFile(z).namelist()
            if f.endswith(ext)
        ]


def parse_imaging_data(paths: List[Tuple[str, str]]) -> Optional[DataFrame]:
    from pandas import concat, to_datetime

    dataframe = DataFrame.from_records(
        paths, columns=["source_zipfile", "source_filename"]
    )

    # Parse imaging metadata from file paths.
    pattern = (
        r"(?P<source_participant_id>P_\w{6})"
        r"_(?P<date>[\d-]+)"
        r"_(?P<source_session_id>HAB_\d{1}.\d{1,2})"
        r"_(?P<protocol>[\w-]+)"
    )

    dataframe = concat(
        [dataframe, dataframe["source_filename"].str.extract(pattern)],
        axis="columns",
    ).assign(date=lambda df: to_datetime(df.date))

    # Map protocol to BIDS entities.
    protocol_to_bids = DataFrame.from_dict(PROTOCOL_TO_BIDS, orient="index")

    dataframe = (
        dataframe.join(protocol_to_bids, on="protocol", how="left")
        .drop(columns="protocol")
        .dropna(subset=["datatype", "modality"])
    )

    if dataframe.empty:
        return

    # Compute BIDS participant ID, session ID and filename.
    dataframe = (
        dataframe.assign(participant_id=source_participant_id_to_bids)
        .drop(columns="source_participant_id")
        .assign(session_id=source_session_id_to_bids)
        .drop(columns="source_session_id")
        .set_index(
            ["participant_id", "session_id", "datatype", "modality", "trc_label"],
            verify_integrity=True,
        )
        .sort_index()
    )

    return dataframe


def write_to_tsv(dataframe: DataFrame, buffer: TextIOBase) -> None:
    """Save the input dataframe to a BIDS-compliant TSV format."""
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def install_nifti(zipfile: str, filename: str, bids_path: str) -> None:
    """Install a NIfTI file from a source archive to the target BIDS path."""
    import fsspec

    fo = fsspec.open(zipfile)
    fs = fsspec.filesystem("zip", fo=fo)
    with fsspec.open(bids_path, mode="wb") as f:
        f.write(fs.cat(filename))


def write_bids(
    sourcedata: Path,
    rawdata: Path,
    imaging_data: DataFrame,
    clinical_data: Dict[str, DataFrame],
):
    import fsspec
    from pandas import notna

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    participants = (
        clinical_data["Demographics"]
        .drop(columns="NP_Age")
        .rename(
            columns={
                "BiologicalSex": "sex",
                "YrsOfEd": "years_of_education",
            }
        )
        .xs("ses-M00", level="session_id")
        .reset_index(level="date", drop=True)
    ).rename(columns=str.lower)

    with fsspec.open(str(rawdata / "dataset_description.json"), mode="wt") as f:
        BIDSDatasetDescription(name="HABS").write(to=f)

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

    for participant_id, dataframe in sessions.groupby(["participant_id"]):
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
                        f"{'_trc-'+row.trc_label if notna(row.trc_label) else ''}"
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
                zipfile=str(sourcedata / row.source_zipfile),
                filename=row.source_filename,
                bids_path=str(bids_basedir / filename),
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
    bids.write_modality_agnostic_files(
        study_name="HABS", readme_data=readme_data, bids_dir=rawdata
    )
