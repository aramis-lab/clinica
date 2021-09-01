from os import PathLike
from typing import BinaryIO, Dict, Iterable, Optional, Union

import pandas as pd
from pandas import DataFrame, Series

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


def find_clinical_data(sourcedata: PathLike) -> Iterable[PathLike]:
    from csv import DictReader
    from pathlib import Path

    path = Path(sourcedata)

    for f in path.rglob("*HABS_DataRelease*.csv"):
        header = DictReader(open(f, encoding="utf-8-sig")).fieldnames
        yield (
            f.name.split("_")[0],
            str(f.relative_to(path)),
            {
                header[0]: "source_participant_id",
                header[1]: "source_session_id",
                header[2]: "date",
            },
        )


def read_clinical_data(path: PathLike, rename_columns: Dict[str, str]) -> DataFrame:
    from pathlib import Path

    from pandas import read_csv

    return (
        read_csv(Path(path))
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


def find_imaging_data(source_data_dir: PathLike) -> Iterable[PathLike]:
    from pathlib import Path
    from zipfile import ZipFile

    path = Path(source_data_dir)
    ext = ".nii.gz"

    for z in path.rglob("*HABS_DataRelease*.zip"):
        yield [
            str(z.relative_to(path) / f)
            for f in ZipFile(z).namelist()
            if f.endswith(ext)
        ]


def parse_imaging_data(paths: PathLike) -> Optional[DataFrame]:
    from pandas import Series, concat, to_datetime

    dataframe = Series(paths, name="source_location")

    # Parse imaging metadata from file paths.
    pattern = (
        r"(?P<source_participant_id>P_\w{6})"
        r"_(?P<date>[\d-]+)"
        r"_(?P<source_session_id>HAB_\d{1}.\d{1,2})"
        r"_(?P<protocol>[\w-]+)"
    )

    dataframe = concat(
        [dataframe, dataframe.str.extract(pattern)], axis="columns"
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


def write_to_tsv(dataframe: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    # Save dataframe as a BIDS-compliant TSV file.
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def install_nifti(source: PathLike, target: PathLike) -> None:
    from pathlib import Path

    import fsspec

    source = Path(source)
    zip_file = source.parent
    nifti_file = source.name

    fo = fsspec.open(str(zip_file))
    fs = fsspec.filesystem("zip", fo=fo)
    with fsspec.open(str(target), mode="wb") as f:
        f.write(fs.cat(nifti_file))


def write_bids(
    sourcedata: PathLike,
    rawdata: PathLike,
    imaging_data: DataFrame,
    clinical_data: Dict[str, DataFrame],
):
    from pathlib import Path

    import fsspec
    from pandas import notna

    participants = (
        clinical_data["Demographics"][
            [
                "NP_Age",
                "BiologicalSex",
                "YrsOfEd",
                "Race",
                "Ethnicity",
                "Holingshead",
                "APOE_haplotype",
                "E4_Status",
            ]
        ]
        .rename(
            columns={
                "NP_Age": "age",
                "BiologicalSex": "sex",
                "YrsOfEd": "education",
                "Race": "race",
                "Ethnicity": "ethnicity",
                "Holingshead": "holingshead",
                "APOE_haplotype": "apoe_haplotype",
                "E4_Status": "e4_status",
            }
        )
        .xs("ses-M00", level="session_id")
        .reset_index(level="date", drop=True)
    )

    participants_file = Path(rawdata) / "participants.tsv"
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
    )

    for participant_id, dataframe in sessions.groupby(["participant_id"]):
        sessions_file = (
            Path(rawdata) / participant_id / f"{participant_id}_sessions.tsv"
        )
        dataframe = dataframe.droplevel("participant_id")
        with fsspec.open(sessions_file, mode="wb") as f:
            write_to_tsv(dataframe, f)

    for grouped_by, dataframe in imaging_data.groupby(["participant_id", "session_id"]):
        participant_id, session_id = grouped_by

        bids_basedir = Path(rawdata) / participant_id / session_id
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
                source=Path(sourcedata) / row.source_location,
                target=Path(bids_basedir) / filename,
            )

        dataframe = dataframe.drop(columns="source_location")

        with fsspec.open(bids_basedir / "scans.tsv", mode="wb") as f:
            write_to_tsv(dataframe, f)
