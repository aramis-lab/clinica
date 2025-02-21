from cmath import nan
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

__all__ = [
    "find_clinical_data",
    "read_imaging_data",
    "merge_imaging_and_clinical_data",
    "prepare_dataset_to_bids_format",
    "write_bids",
]


def find_clinical_data(clinical_data_directory: Path) -> pd.DataFrame:
    from clinica.utils.stream import cprint

    image_data_file = [f for f in clinical_data_directory.glob("*.csv")]
    if len(image_data_file) == 0:
        raise FileNotFoundError(
            f"Clinical data not found or incomplete in {clinical_data_directory}."
        )
    if len(image_data_file) > 1:
        raise ValueError(
            f"Too many data files found in {clinical_data_directory}, expected only one."
        )
    df_clinical = pd.read_csv(
        image_data_file[0],
        usecols=["eid", "31-0.0", "34-0.0", "21022-0.0", "21003-2.0", "21003-3.0"],
    )
    cprint(msg="All clinical data have been found", lvl="info")
    required_columns = {
        "eid",
        "31-0.0",
        "34-0.0",
        "21022-0.0",
        "21003-2.0",
        "21003-3.0",
    }
    if (
        len(missing_columns := required_columns.difference(set(df_clinical.columns)))
        > 0
    ):
        raise ValueError(f"Required columns {missing_columns} not found")
    df_clinical = df_clinical.rename(
        columns={
            "31-0.0": "sex_f31_0_0",
            "34-0.0": "year_of_birth_f34_0_0",
            "21022-0.0": "age_at_recruitment_f21022_0_0",
            "21003-2.0": "age_when_attended_assessment_centre_f21003_2_0",
            "21003-3.0": "age_when_attended_assessment_centre_f21003_3_0",
        }
    )
    return df_clinical


def read_imaging_data(imaging_data_directory: Path) -> pd.DataFrame:
    source_path_series_nifti = pd.Series(
        _find_imaging_data(imaging_data_directory), name="source_path"
    )
    source_path_series_dicom = pd.Series(
        _find_dicom_data(imaging_data_directory), name="source_path"
    )
    # Rationales for this choice of files can be found in the documentation:
    # https://github.com/aramis-lab/clinica/blob/dev/docs/Converters/UKBtoBIDS.md
    # list of the files we want to build the bids for each modality
    file_mod_list = [
        "T1.nii.gz",
        "T2_FLAIR_orig_defaced.nii.gz",
        "AP.nii.gz",
        "PA.nii.gz",
        "rfMRI.dcm",
        "tfMRI.dcm",
        "SWI.nii.gz",
    ]

    dataframe_nifti = pd.DataFrame.from_records(
        source_path_series_nifti, columns=["source_zipfile", "source_filename"]
    )
    dataframe_nifti = dataframe_nifti[
        dataframe_nifti["source_filename"].str.contains("unusable") != True  # noqa
    ]
    dataframe_dicom = pd.DataFrame.from_records(
        source_path_series_dicom, columns=["source_zipfile", "source_filename"]
    )
    dataframe = pd.concat([dataframe_nifti, dataframe_dicom])
    filename = dataframe["source_filename"].apply(lambda x: x.name).rename("filename")
    dataframe = dataframe[filename.isin(file_mod_list)]
    filename = filename[filename.isin(file_mod_list)]
    if dataframe.empty:
        raise ValueError(
            f"No imaging data were found in the provided folder: {imaging_data_directory}, "
            "or they are not handled by Clinica. Please check your data."
        )
    dataframe["source_zipfile"] = dataframe["source_zipfile"].apply(lambda x: str(x))
    split_zipfile = dataframe["source_zipfile"].str.split("_", expand=True)
    split_zipfile = split_zipfile.rename(
        {0: "source_id", 1: "modality_num", 2: "source_sessions_number"}, axis="columns"
    )  # .drop_duplicates()
    split_zipfile = split_zipfile.astype({"source_id": "int64"})

    df_source = pd.concat(
        [
            dataframe,
            filename,
            split_zipfile,
        ],
        axis=1,
    )
    return df_source


def _find_dicom_data(path_to_source_data: Path) -> Iterable[list[Path, Path]]:
    """
    This function finds the paths to the dicoms, only for fMRI, since we use the niftis for the other modalities.
    More information here: https://github.com/aramis-lab/clinica/blob/dev/docs/Converters/UKBtoBIDS.md
    """
    for z in path_to_source_data.rglob("*.zip"):
        if z.name.split("_")[1] == "20217":
            yield [z.relative_to(path_to_source_data), Path("fMRI/tfMRI.dcm")]
        elif z.name.split("_")[1] == "20225":
            yield [z.relative_to(path_to_source_data), Path("fMRI/rfMRI.dcm")]


def _find_imaging_data(path_to_source_data: Path) -> Iterable[list[Path, Path]]:
    from zipfile import ZipFile

    for z in path_to_source_data.rglob("*.zip"):
        for f in ZipFile(z).namelist():
            if f.endswith(".nii.gz"):
                yield [z.relative_to(path_to_source_data), Path(f)]


def merge_imaging_and_clinical_data(
    imaging_data: pd.DataFrame, clinical_data: pd.DataFrame
) -> pd.DataFrame:
    """Merges the imaging and clinical dataframes based on the subject id."""
    merged_data = clinical_data.merge(
        imaging_data, how="inner", right_on="source_id", left_on="eid"
    )
    return merged_data


def _complete_clinical(df_clinical: pd.DataFrame) -> pd.DataFrame:
    """This function uses the existing data to create the columns needed for
    the bids hierarchy (subject_id, ses, age_at _sessions, etc.)"""
    from clinica.converters.bids_utils import StudyName, bids_id_factory

    df_clinical = df_clinical.assign(
        participant_id=lambda df: df.source_id.astype("str").apply(
            lambda x: bids_id_factory(StudyName.UKB).from_original_study_id(x)
        )
    )
    df_clinical = df_clinical.assign(
        sessions=lambda df: "ses-" + df.source_sessions_number.astype("str")
    )
    df_clinical = df_clinical.join(
        df_clinical.filename.map(
            {
                "T2_FLAIR_orig_defaced.nii.gz": {
                    "datatype": "anat",
                    "modality": "FLAIR",
                    "suffix": "FLAIR",
                    "sidecars": ["T2_FLAIR/T2_FLAIR.json"],
                    "task": "",
                    "dir": "",
                },
                "T1.nii.gz": {
                    "datatype": "anat",
                    "modality": "T1W",
                    "suffix": "T1w",
                    "sidecars": ["T1/T1.json"],
                    "task": "",
                    "dir": "",
                },
                "AP.nii.gz": {
                    "datatype": "dwi",
                    "modality": "dwi",
                    "suffix": "dwi",
                    "sidecars": [
                        "dMRI/raw/AP.json",
                        "dMRI/raw/AP.bval",
                        "dMRI/raw/AP.bvec",
                    ],
                    "task": "",
                    "dir": "_dir-AP",
                },
                "PA.nii.gz": {
                    "datatype": "dwi",
                    "modality": "dwi",
                    "suffix": "dwi",
                    "sidecars": [
                        "dMRI/raw/PA.json",
                        "dMRI/raw/PA.bval",
                        "dMRI/raw/PA.bvec",
                    ],
                    "task": "",
                    "dir": "_dir-PA",
                },
                "rfMRI.dcm": {
                    "datatype": "func",
                    "modality": "rsfmri",
                    "suffix": "bold",
                    "sidecars": [],
                    "task": "_task-rest",
                    "dir": "",
                },
                "tfMRI.dcm": {
                    "datatype": "func",
                    "modality": "tfmri",
                    "suffix": "bold",
                    "sidecars": [],
                    "task": "_task-facesshapesemotion",
                    "dir": "",
                },
                "SWI.nii.gz": {
                    "datatype": "swi",
                    "modality": "swi",
                    "suffix": "swi",
                    "sidecars": [],
                    "task": "",
                    "dir": "",
                },
            }
        ).apply(pd.Series)
    )
    df_clinical = df_clinical.assign(year_of_birth=lambda df: df.year_of_birth_f34_0_0)

    df_clinical = df_clinical.assign(age=lambda df: df.age_at_recruitment_f21022_0_0)
    df_clinical = df_clinical.assign(
        age_at_first_session=lambda df: df.age_when_attended_assessment_centre_f21003_2_0
    )
    df_clinical["age_at_session"] = df_clinical.apply(
        lambda df: _select_sessions(df), axis=1
    )
    df_clinical = df_clinical[df_clinical["age_at_session"].notna()]
    df_clinical = df_clinical.assign(
        sessions_month=lambda df: (df.age_at_session - df.age_at_first_session) * 12
    )
    df_clinical = df_clinical.assign(
        sessions=lambda df: df.sessions_month.map(lambda x: f"ses-M{int(x):03d}")
    )
    df_clinical = df_clinical.assign(
        bids_filename=lambda df: (
            "sub-UKB"
            + df.source_id.astype("str")
            + "_"
            + df.sessions.astype("str")
            + df.dir.astype("str")
            + df.task.astype("str")
            + "_"
            + df.suffix
        )
    )
    df_clinical = df_clinical.join(
        df_clinical.sex_f31_0_0.astype("str")
        .map(
            {
                "0": {"sex": "F"},
                "1": {"sex": "M"},
            }
        )
        .apply(pd.Series)
    )
    df_clinical = df_clinical.assign(
        bids_full_path=lambda df: df.participant_id
        + "/"
        + df.sessions
        + "/"
        + df.datatype
        + "/"
        + df.bids_filename
    )
    return df_clinical


def prepare_dataset_to_bids_format(
    df_clinical: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    df_clinical = _complete_clinical(df_clinical)
    df_specifications = pd.read_csv(
        Path(__file__).parent / "specifications.csv", sep=";"
    )
    df_clinical = df_clinical.set_index(
        ["participant_id", "sessions", "modality", "bids_filename"],
        verify_integrity=True,
    )

    return {
        col: df_clinical.filter(items=list(df_specifications[col]))
        for col in ["participants", "sessions", "scans"]
    }


def write_bids(
    to: Path,
    participants: pd.DataFrame,
    sessions: pd.DataFrame,
    scans: pd.DataFrame,
    dataset_directory: Path,
) -> None:
    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters.bids_dataset_description import BIDSDatasetDescription
    from clinica.converters.bids_utils import StudyName, write_to_tsv

    fs = LocalFileSystem(auto_mkdir=True)

    participants = participants.droplevel(
        ["sessions", "modality", "bids_filename"]
    ).drop_duplicates()

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(
            str(to / "dataset_description.json"), "w"
        ) as dataset_description_file:
            BIDSDatasetDescription(name=StudyName.UKB).write(
                to=dataset_description_file
            )
        with fs.open(str(to / "participants.tsv"), "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby("participant_id"):
        sessions = data_frame.droplevel(
            ["participant_id", "modality", "bids_filename"]
        ).drop_duplicates()

        sessions_filepath = to / str(participant_id) / f"{participant_id}_sessions.tsv"
        with fs.open(str(sessions_filepath), "w") as sessions_file:
            write_to_tsv(sessions, sessions_file)
    scans = scans.set_index(["bids_full_path"], verify_integrity=True)
    for bids_full_path, metadata in scans.iterrows():
        if metadata["modality_num"] != "20217" and metadata["modality_num"] != "20225":
            _copy_file_to_bids(
                zipfile=dataset_directory / metadata["source_zipfile"],
                filenames=[metadata["source_filename"]] + metadata["sidecars"],
                bids_path=to / bids_full_path,
            )
        else:
            _convert_dicom_to_nifti(
                zipfiles=dataset_directory / metadata["source_zipfile"],
                bids_path=to / bids_full_path,
            )
            if metadata["modality_num"] == "20217":
                _import_event_tsv(bids_path=to)
    return


def _copy_file_to_bids(zipfile: Path, filenames: List[Path], bids_path: Path) -> None:
    """Install the requested files in the BIDS  dataset."""
    import fsspec

    fo = fsspec.open(str(zipfile))
    fs = fsspec.filesystem("zip", fo=fo)
    for filename in filenames:
        filename = str(filename)
        if fs.exists(filename):
            bids_path_extension = str(bids_path) + "." + (filename.split(".", 1)[1])
            with fsspec.open(bids_path_extension, mode="wb") as f:
                f.write(fs.cat(filename))


def _convert_dicom_to_nifti(zipfiles: Path, bids_path: Path) -> None:
    """Install the requested files in the BIDS  dataset.
    First, the dicom is extracted in a temporary directory
    Second, the dicom extracted is converted in the right place using dcm2niix"""
    import json
    import subprocess
    import tempfile
    import zipfile
    from pathlib import PurePath

    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)

    zf = zipfile.ZipFile(zipfiles)
    try:
        bids_path.parent.mkdir(exist_ok=True, parents=True)
    except OSError:
        # Folder already created with previous instance
        pass
    with tempfile.TemporaryDirectory() as tempdir:
        zf.extractall(tempdir)
        command = [
            "dcm2niix",
            "-w",
            "0",
        ]
        command += ["-9", "-z", "y"]
        command += ["-b", "y", "-ba", "y"]
        command += [tempdir]
        subprocess.run(command)
        fmri_image_path = _find_largest_imaging_data(Path(tempdir))
        fmri_image_path = fmri_image_path or ""
        fs.copy(str(fmri_image_path), str(bids_path) + ".nii.gz")
        fs.copy(
            str(fmri_image_path.parent)
            + "/"
            + PurePath(PurePath(fmri_image_path.name).stem).stem
            + ".json",
            str(bids_path) + ".json",
        )
    # Add the taskname to the json
    with open(str(bids_path) + ".json", "r+") as f:
        json_file = json.load(f)
        json_file["TaskName"] = "facesshapesemotion"
        f.seek(0)
        json.dump(json_file, f, indent=4)


def _select_sessions(x: pd.DataFrame) -> Optional[pd.Series]:
    from clinica.utils.stream import cprint

    if (
        x["source_sessions_number"] == "2"
        and x.age_when_attended_assessment_centre_f21003_2_0 != nan
    ):
        return x.age_when_attended_assessment_centre_f21003_2_0
    elif (
        x["source_sessions_number"] == "2"
        and x.age_when_attended_assessment_centre_f21003_2_0 == nan
    ):
        cprint(
            msg=(
                f"The subject {x.eid} doesn't have the age for the imaging session "
                "number one (age_when_attended_assessment_centre_f21003_2_0)."
                f"It will not be converted. To have it converted, please update your clinical data."
            ),
            lvl="warning",
        )
        return None
    elif (
        x["source_sessions_number"] == "3"
        and x.age_when_attended_assessment_centre_f21003_3_0 != nan
    ):
        return x.age_when_attended_assessment_centre_f21003_3_0
    elif (
        x["source_sessions_number"] == "3"
        and x.age_when_attended_assessment_centre_f21003_3_0 == nan
    ):
        cprint(
            msg=(
                f"The subject {x.eid} doesn't have the age for the imaging session "
                "number two (age_when_attended_assessment_centre_f21003_3_0)."
                f"It will not be converted. To have it converted, please update your clinical data."
            ),
            lvl="warning",
        )
        return None


def _import_event_tsv(bids_path: Path) -> None:
    """Import the csv containing the events' information."""
    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)
    event_tsv = (
        Path(__file__).parents[3]
        / "resources"
        / "fmri"
        / "task-facesshapesemotion_events.tsv"
    )
    if not event_tsv.exists():
        raise FileNotFoundError(
            f"Could not find file {event_tsv} which is an internal file of Clinica."
        )
    bids_path_extension = bids_path / "task-facesshapesemotion_events.tsv"
    fs.copy(str(event_tsv), str(bids_path_extension))


def _find_largest_imaging_data(path_to_source_data: Path) -> Optional[Path]:
    import os

    weight = 0
    file_to_return = None
    for z in path_to_source_data.rglob("*.nii.gz"):
        if os.path.getsize(z) > weight:
            weight = os.path.getsize(z)
            file_to_return = z
    if file_to_return:
        return file_to_return
    return None
