from os import PathLike
from typing import BinaryIO, Dict, Iterable, List, Optional, Tuple, Union

from pandas import DataFrame


def find_clinical_data(clinical_data_directory: PathLike) -> Optional[DataFrame]:
    from pathlib import Path

    from pandas import read_excel

    dataframe = None
    for path in Path(clinical_data_directory).rglob("*.xlsx"):
        try:
            dataframe = (
                read_excel(path, index_col=[0, 4])  # noqa
                .rename(columns=lambda x: x.lower().replace(" ", "_"))
                .rename_axis(index=lambda x: x.lower().replace(" ", "_"))
                .convert_dtypes()
                .sort_index()
            )
        except (ValueError, IndexError):
            continue

    return dataframe


def read_clinical_data(clinical_data_directory: PathLike) -> DataFrame:
    import pandas as pd

    dataframe = find_clinical_data(clinical_data_directory)

    if dataframe is None:
        raise FileNotFoundError("Clinical data not found")

    # Compute participant and session IDs.
    dataframe = dataframe.rename_axis(
        index={"loni_id": "participant_id", "visit_number": "session_id"}
    )
    dataframe.index = dataframe.index.map(
        lambda x: (f"sub-NIFD{x[0].replace('_', '')}", f"ses-M{(6 * (x[1] - 1)):02d}")
    )

    # Keep relevant columns and rename them.
    dataframe = (
        dataframe[["dx", "site", "education", "race", "cdr_box_score", "mmse_tot"]]
        .rename(columns={"dx": "diagnosis", "cdr_box_score": "cdr", "mmse_tot": "mmse"})
        .astype(
            dtype={
                "diagnosis": pd.CategoricalDtype(
                    ["BV", "CON", "L_SD", "PATIENT (OTHER)", "PNFA", "SV"]
                ),
                "site": pd.CategoricalDtype(["UCSF", "MAYO", "MGH"]),
                "education": pd.Int64Dtype(),
                "race": pd.Int64Dtype(),
                "cdr": pd.Float64Dtype(),
                "mmse": pd.Float64Dtype(),
            }
        )
        .replace({"education": {99: pd.NA}, "race": {50: pd.NA, 99: pd.NA}})
    )

    # Keep positive MMSE values only.
    dataframe.mmse = dataframe.mmse.mask(dataframe.mmse < 0)

    return dataframe


def find_collection_data(imaging_data_directory: PathLike) -> Optional[DataFrame]:
    from pathlib import Path

    from pandas import read_csv

    dataframe = None
    for path in Path(imaging_data_directory).rglob("*.csv"):
        try:
            dataframe = (
                read_csv(
                    path,
                    index_col="Image Data ID",
                    parse_dates=["Acq Date"],
                )
                .rename(columns=lambda x: x.lower().replace(" ", "_"))
                .rename_axis(index=lambda x: x.lower().replace(" ", "_"))
                .convert_dtypes()
                .sort_index()
            )
        except (ValueError, IndexError):
            continue

    return dataframe


def find_imaging_data(imaging_data_directory: PathLike) -> Iterable[Tuple[str, str]]:
    import re
    from pathlib import Path

    # Pattern for extracting the image data ID from the NIFD files.
    pattern = re.compile(r"(I\d{6})")

    def find_files(in_: PathLike) -> Iterable[Path]:
        return filter(lambda x: x.is_file(), Path(in_).rglob("NIFD*.*"))

    def extract_id_with_dir(files: Iterable[Path]) -> Tuple[str, str]:
        for f in files:
            found = pattern.search(f.name)
            if found:
                yield found.group(1), str(f.parent)

    for image_data_id, source_dir in set(
        sorted(extract_id_with_dir(find_files(imaging_data_directory)))
    ):
        yield image_data_id, source_dir


def read_imaging_data(imaging_data_directory: PathLike) -> DataFrame:
    from pandas import DataFrame

    try:
        imaging_data = DataFrame.from_records(
            data=find_imaging_data(imaging_data_directory),
            columns=["image_data_id", "source_dir"],
            index="image_data_id",
        ).convert_dtypes()
    except TypeError:
        raise FileNotFoundError("No imaging data found")

    collection_data = find_collection_data(imaging_data_directory)

    if collection_data is None:
        raise FileNotFoundError("No collection data found")

    return collection_data.join(imaging_data)


def parse_mri_description(description: str) -> Optional[Dict[str, Optional[str]]]:
    description = description.lower().replace("-", "")

    if "mprage" in description:
        return {
            "datatype": "anat",
            "suffix": "T1w",
            "trc_label": None,
            "rec_label": None,
        }
    elif "flair" in description:
        return {
            "datatype": "anat",
            "suffix": "FLAIR",
            "trc_label": None,
            "rec_label": None,
        }
    elif "t2" in description:
        return {
            "datatype": "anat",
            "suffix": "T2w",
            "trc_label": None,
            "rec_label": None,
        }
    elif "asl" in description:
        return {
            "datatype": "anat",
            "suffix": "PDw",
            "trc_label": None,
            "rec_label": None,
        }
    elif any([x in description for x in ["mt1", "gradwarp", "n3m"]]):
        return {
            "datatype": "anat",
            "suffix": "T1w",
            "trc_label": None,
            "rec_label": None,
        }
    else:
        return None


def parse_pet_description(description: str) -> Optional[Dict[str, str]]:
    import re

    match = re.search(r"3D:(\w+):(\w+)", description)

    if match:
        return {
            "datatype": "pet",
            "suffix": "pet",
            "trc_label": "11CPIB" if "PIB" in match.group(1) else "18FFDG",
            "rec_label": "IR" if "IR" in match.group(2) else "RP",
        }
    else:
        return None


def parse_preprocessing(description: str) -> dict:
    description = description.lower()

    return {
        "gradwarp": any([x in description for x in ["gradwarp", "dis3d"]]),
        "n3": "n3m" in description,
    }


def dataset_to_bids(
    imaging_data: DataFrame,
    clinical_data: Optional[DataFrame] = None,
) -> Tuple[DataFrame, DataFrame, DataFrame]:
    from pandas import Series

    # Parse preprocessing information from scan descriptions.
    preprocessing = imaging_data.description.apply(parse_preprocessing).apply(Series)

    # Parse BIDS entities from scan descriptions.
    bids = (
        imaging_data.apply(
            lambda x: parse_pet_description(x.description)
            if x.modality == "PET"
            else parse_mri_description(x.description),
            axis=1,
        )
        .dropna()
        .apply(Series)
    )

    # Compute quality metric for each scan:
    # - MRI: Applied preprocessing (0: None, 1: GradWarp, 2: N3)
    # - PET: Reconstruction method (0: Fourier, 1: Iterative)
    quality = (
        preprocessing.sum(axis=1)
        + bids.rec_label.apply(lambda x: 1 if x == "IR" else 0)
    ).rename("quality")

    # Select one scan per BIDS modality based on quality metric.
    subset = ["subject", "visit", "datatype", "suffix", "trc_label"]
    scans = (
        bids.join(quality)
        .join(imaging_data)
        .sort_values(by=subset + ["quality"])
        .drop_duplicates(subset=subset, keep="last")
        .drop(columns="quality")
    )

    # Compute the BIDS-compliant participant, session and scan IDs.
    scans = scans.assign(
        participant_id=lambda df: df.subject.apply(
            lambda x: f"sub-NIFD{x.replace('_', '')}"
        ),
        session_id=lambda df: df.visit.apply(lambda x: f"ses-M{(6 * (x - 1)):02d}"),
        filename=lambda df: df.apply(
            lambda x: f"{x.participant_id}/{x.session_id}/{x.datatype}/"
            f"{x.participant_id}_{x.session_id}"
            f"{'_trc-' + x.trc_label if x.trc_label else ''}"
            f"{'_rec-' + x.rec_label if x.rec_label else ''}"
            f"_{x.suffix}.nii.gz",
            axis=1,
        ),
    )

    # Prepare subjects manifest.
    subjects = (
        scans[["participant_id", "session_id", "group", "sex", "age"]]
        .sort_values(by=["participant_id", "session_id"])
        .drop(columns="session_id")
        .drop_duplicates(subset="participant_id")
        .set_index(["participant_id"], verify_integrity=True)
        .sort_index()
    ).join(
        clinical_data.xs("ses-M00", level="session_id")[
            ["diagnosis", "site", "education", "race"]
        ]
    )

    # Prepare sessions manifest.
    sessions = (
        scans[["participant_id", "session_id", "acq_date", "age"]]
        .rename(columns={"acq_date": "date"})
        .drop_duplicates()
        .set_index(["participant_id", "session_id"], verify_integrity=True)
        .sort_index()
    ).join(clinical_data[["cdr", "mmse"]])

    # Prepare scans manifest.
    scans = scans[["filename", "source_dir", "format"]].set_index(
        "filename", verify_integrity=True
    )

    return subjects, sessions, scans


def write_to_tsv(dataframe: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    # Save dataframe as a BIDS-compliant TSV file.
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def convert_dicom(sourcedata_dir: PathLike, bids_filename: PathLike) -> None:
    from pathlib import PurePath

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_utils import run_dcm2niix

    output_fmt = str(PurePath(bids_filename).name).replace(".nii.gz", "")
    output_dir = str(PurePath(bids_filename).parent)

    # Ensure output directory is empty.
    fs = LocalFileSystem()
    if fs.exists(output_dir):
        fs.rm(output_dir, recursive=True)
    fs.makedirs(output_dir)

    # Run conversion with dcm2niix with anonymization and maximum compression.
    run_dcm2niix(
        input_dir=str(sourcedata_dir),
        output_dir=output_dir,
        output_fmt=output_fmt,
        compress=True,
        bids_sidecar=True,
    )


def install_nifti(sourcedata_dir: PathLike, bids_filename: PathLike) -> None:
    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)
    source_file = fs.open(fs.ls(str(sourcedata_dir))[0], mode="rb")
    target_file = fs.open(str(bids_filename), mode="wb", compression="gzip")

    with source_file as sf, target_file as tf:
        tf.write(sf.read())


def write_bids(
    to: PathLike,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
) -> List[PathLike]:
    from pathlib import PurePath

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    to = PurePath(to)
    fs = LocalFileSystem(auto_mkdir=True)

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(
            str(to / "dataset_description.json"), "w"
        ) as dataset_description_file:
            BIDSDatasetDescription(name="NIFD").write(to=dataset_description_file)
        with fs.open(str(to / "participants.tsv"), "w") as participant_file:
            write_to_tsv(participants, participant_file)

        for participant_id, sessions_group in sessions.groupby("participant_id"):
            participant_id = str(participant_id)
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(str(sessions_filepath), "w") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    # Perform import of imaging data next.
    for filename, metadata in scans.iterrows():
        filename = str(filename)
        if metadata.format == "DCM":
            convert_dicom(
                sourcedata_dir=metadata.source_dir, bids_filename=to / filename
            )
        else:
            install_nifti(
                sourcedata_dir=metadata.source_dir, bids_filename=to / filename
            )

    return scans.index.to_list()
