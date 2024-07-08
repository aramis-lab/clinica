from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

__all__ = ["read_clinical_data", "read_imaging_data", "dataset_to_bids", "write_bids"]


def _find_clinical_data(clinical_data_directory: Path) -> Optional[pd.DataFrame]:
    dataframe = None
    for file in clinical_data_directory.rglob("*.xlsx"):
        try:
            dataframe = (
                pd.read_excel(file, index_col=[0, 4])
                .rename(columns=lambda x: x.lower().replace(" ", "_"))
                .rename_axis(index=lambda x: x.lower().replace(" ", "_"))
                .convert_dtypes()
                .sort_index()
            )
        except (ValueError, IndexError):
            continue

    return dataframe


def read_clinical_data(clinical_data_directory: Path) -> pd.DataFrame:
    if (dataframe := _find_clinical_data(clinical_data_directory)) is None:
        raise FileNotFoundError("Clinical data not found")
    # Compute participant and session IDs.
    dataframe = dataframe.rename_axis(
        index={"loni_id": "participant_id", "visit_number": "session_id"}
    )
    dataframe.index = dataframe.index.map(
        lambda x: (f"sub-NIFD{x[0].replace('_', '')}", f"ses-M{(6 * (x[1] - 1)):03d}")
    )

    # todo : use func here

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


def _find_collection_data(imaging_data_directory: Path) -> Optional[pd.DataFrame]:
    dataframe = None
    for file in imaging_data_directory.rglob("*.csv"):
        try:
            dataframe = (
                pd.read_csv(
                    file,
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


def _find_imaging_data(imaging_data_directory: Path) -> Iterable[tuple[str, str]]:
    for image_data_id, source_dir in set(
        sorted(_extract_id_with_dir(_find_files(imaging_data_directory)))
    ):
        yield image_data_id, source_dir


def _extract_id_with_dir(files: Iterable[Path]) -> tuple[str, str]:
    """Extract the image data ID from the NIFD files."""
    import re

    pattern = re.compile(r"(I\d{6})")
    for file in files:
        if (found := pattern.search(file.name)) is not None:
            yield found.group(1), str(file.parent)


def _find_files(in_: Path) -> Iterable[Path]:
    return filter(lambda x: x.is_file(), Path(in_).rglob("NIFD*.*"))


def read_imaging_data(imaging_data_directory: Path) -> pd.DataFrame:
    try:
        imaging_data = pd.DataFrame.from_records(
            data=_find_imaging_data(imaging_data_directory),
            columns=["image_data_id", "source_dir"],
            index="image_data_id",
        ).convert_dtypes()
    except TypeError:
        raise FileNotFoundError("No imaging data found")
    if (collection_data := _find_collection_data(imaging_data_directory)) is None:
        raise FileNotFoundError("No collection data found")

    return collection_data.join(imaging_data)


def _parse_mri_description(description: str) -> Optional[dict[str, Optional[str]]]:
    description = description.lower().replace("-", "")

    if (suffix := _infer_bids_suffix_from_description(description)) is None:
        return None
    return {
        "datatype": "anat",
        "suffix": suffix,
        "trc_label": None,
        "rec_label": None,
    }


def _infer_bids_suffix_from_description(description: str) -> Optional[str]:
    if any([x in description for x in ("mprage", "mt1", "gradwarp", "n3m")]):
        return "T1w"
    if "flair" in description:
        return "FLAIR"
    if "t2" in description:
        return "T2w"
    if "asl" in description:
        return "PDw"
    return None


def _parse_pet_description(description: str) -> Optional[dict[str, str]]:
    import re

    match = re.search(r"3D:(\w+):(\w+)", description)

    if match:
        return {
            "datatype": "pet",
            "suffix": "pet",
            "trc_label": "11CPIB" if "PIB" in match.group(1) else "18FFDG",
            "rec_label": "IR" if "IR" in match.group(2) else "RP",
        }
    return None


def _parse_preprocessing(description: str) -> dict[str, bool]:
    description = description.lower()

    return {
        "gradwarp": any([x in description for x in ["gradwarp", "dis3d"]]),
        "n3": "n3m" in description,
    }


def dataset_to_bids(
    imaging_data: pd.DataFrame,
    clinical_data: Optional[pd.DataFrame] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    from clinica.iotools.bids_utils import StudyName, _id_factory

    # Parse preprocessing information from scan descriptions.
    preprocessing = imaging_data.description.apply(_parse_preprocessing).apply(
        pd.Series
    )
    # Parse BIDS entities from scan descriptions.
    bids = (
        imaging_data.apply(
            lambda x: _parse_pet_description(x.description)
            if x.modality == "PET"
            else _parse_mri_description(x.description),
            axis=1,
        )
        .dropna()
        .apply(pd.Series)
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
            lambda x: _id_factory(StudyName.NIFD).from_original_study_id(x)
        ),
        session_id=lambda df: df.visit.apply(
            lambda x: (
                "ses-M" + x.strip("M").zfill(3)
                if pd.api.types.is_string_dtype(x)
                else f"ses-M{(6 * (x - 1)):03d}"
            )
        ),
        filename=lambda df: df.apply(
            lambda x: (
                f"{x.participant_id}/{x.session_id}/{x.datatype}/"
                f"{x.participant_id}_{x.session_id}"
                f"{'_trc-' + x.trc_label if x.trc_label else ''}"
                f"{'_rec-' + x.rec_label if x.rec_label else ''}"
                f"_{x.suffix}.nii.gz"
            ),
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
        clinical_data.xs("ses-M000", level="session_id")[
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


def _convert_dicom(sourcedata_dir: Path, bids_filename: Path) -> None:
    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_utils import run_dcm2niix

    output_fmt = str(bids_filename.name).replace(".nii.gz", "")
    output_dir = bids_filename.parent

    # Ensure output directory is empty.
    fs = LocalFileSystem()
    if fs.exists(output_dir):
        fs.rm(str(output_dir), recursive=True)
    fs.makedirs(str(output_dir))

    # Run conversion with dcm2niix with anonymization and maximum compression.
    run_dcm2niix(
        input_dir=sourcedata_dir,
        output_dir=output_dir,
        output_fmt=output_fmt,
        compress=True,
        bids_sidecar=True,
    )


def _install_nifti(sourcedata_dir: Path, bids_filename: Path) -> None:
    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)
    source_file = fs.open(fs.ls(str(sourcedata_dir))[0], mode="rb")
    target_file = fs.open(str(bids_filename), mode="wb", compression="gzip")

    with source_file as sf, target_file as tf:
        tf.write(sf.read())


def write_bids(
    to: Path,
    participants: pd.DataFrame,
    sessions: pd.DataFrame,
    scans: pd.DataFrame,
) -> list[Path]:
    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription
    from clinica.iotools.bids_utils import StudyName, write_to_tsv

    fs = LocalFileSystem(auto_mkdir=True)

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(
            str(to / "dataset_description.json"), "w"
        ) as dataset_description_file:
            BIDSDatasetDescription(name=StudyName.NIFD).write(
                to=dataset_description_file
            )
        with fs.open(str(to / "participants.tsv"), "w") as participant_file:
            write_to_tsv(participants, participant_file)

        for participant_id, sessions_group in sessions.groupby("participant_id"):
            participant_id = str(participant_id)
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(str(sessions_filepath), "w") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    for filename, metadata in scans.iterrows():
        filename = str(filename)
        (_convert_dicom if metadata.format == "DCM" else _install_nifti)(
            sourcedata_dir=metadata.source_dir, bids_filename=to / filename
        )

    return scans.index.to_list()
