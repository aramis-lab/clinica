from cmath import nan
from os import PathLike
from typing import Dict, Iterable, List, Optional

from pandas import DataFrame, Series


def find_clinical_data(
    clinical_data_directory: PathLike,
) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    # Read the xls
    try:
        image_data_file = list(Path(clinical_data_directory).glob("*.csv"))
    except StopIteration:
        raise FileNotFoundError("Clinical data file not found.")

    # Assign the dfs
    if len(image_data_file) == 0:
        raise FileNotFoundError("Clinical data not found or incomplete. Aborting")
    if len(image_data_file) > 1:
        raise ValueError("Too many data files found, expected one. Aborting.")
    df_clinical = pd.read_csv(
        str(image_data_file[0]),
        usecols=["eid", "31-0.0", "34-0.0", "21022-0.0", "21003-2.0", "21003-3.0"],
    )
    cprint(msg="All clinical data have been found", lvl="info")
    required_columns = [
        "eid",
        "31-0.0",
        "34-0.0",
        "21022-0.0",
        "21003-2.0",
        "21003-3.0",
    ]
    missing_columns = set(required_columns).difference(set(df_clinical.columns))
    if missing_columns:
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


def read_imaging_data(imaging_data_directory: PathLike) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    source_path_series_nifti = pd.Series(
        find_imaging_data(imaging_data_directory), name="source_path"
    )
    source_path_series_dicom = pd.Series(
        find_dicom_data(imaging_data_directory), name="source_path"
    )
    # Rationales for this choice of files can be found in the documentation:
    # https://github.com/aramis-lab/clinica/blob/dev/docs/Converters/UKBtoBIDS.md
    # list of the files we want to build the bids for each modality
    file_mod_list = [
        "T1.nii.gz",
        "T2_FLAIR.nii.gz",
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
        ~dataframe_nifti["source_filename"].str.contains("unusable")
    ]
    dataframe_dicom = pd.DataFrame.from_records(
        source_path_series_dicom, columns=["source_zipfile", "source_filename"]
    )
    dataframe = pd.concat([dataframe_nifti, dataframe_dicom], 0)
    filename = (
        dataframe["source_filename"]
        .apply(lambda x: Path(str(x)).name)
        .rename("filename")
    )
    dataframe = dataframe[filename.isin(file_mod_list)]
    filename = filename[filename.isin(file_mod_list)]
    if dataframe.empty:
        raise ValueError(
            f"No imaging data were found in the provided folder: {imaging_data_directory}, "
            "or they are not handled by Clinica. Please check your data."
        )
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


def find_dicom_data(path_to_source_data: PathLike) -> Iterable[PathLike]:
    """
    This function finds the paths to the dicoms, only for fMRI, since we use the niftis for the other modalities.
    More information here: https://github.com/aramis-lab/clinica/blob/dev/docs/Converters/UKBtoBIDS.md"""
    from pathlib import Path

    for z in Path(path_to_source_data).rglob("*.zip"):
        if str(z.name).split("_")[1] == "20217":
            yield [str(z.relative_to(path_to_source_data)), "fMRI/tfMRI.dcm"]
        elif str(z.name).split("_")[1] == "20225":
            yield [str(z.relative_to(path_to_source_data)), "fMRI/rfMRI.dcm"]


def find_imaging_data(path_to_source_data: PathLike) -> Iterable[PathLike]:
    from pathlib import Path
    from zipfile import ZipFile

    for z in Path(path_to_source_data).rglob("*.zip"):
        for f in ZipFile(z).namelist():
            if f.endswith(".nii.gz"):
                yield [str(z.relative_to(path_to_source_data)), f]


def intersect_data(df_source: DataFrame, df_clinical_data: DataFrame) -> DataFrame:
    """
    This function merges the two dataframes given as inputs based on the subject id
    """

    df_clinical = df_clinical_data.merge(
        df_source, how="inner", right_on="source_id", left_on="eid"
    )
    return df_clinical


def complete_clinical(df_clinical: DataFrame) -> DataFrame:
    """This function uses the existing data to create the columns needed for
    the bids hierarchy (subject_id, ses, age_at _sessions, ect.)"""
    import pandas as pd

    df_clinical = df_clinical.assign(
        participant_id=lambda df: ("sub-UKB" + df.source_id.astype("str"))
    )
    df_clinical = df_clinical.assign(
        sessions=lambda df: "ses-" + df.source_sessions_number.astype("str")
    )
    df_clinical = df_clinical.join(
        df_clinical.filename.map(
            {
                "T2_FLAIR.nii.gz": {
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
        lambda df: select_sessions(df), axis=1
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


def dataset_to_bids(df_clinical: DataFrame) -> Dict[str, DataFrame]:

    import os

    import pandas as pd

    # open the reference for building the tsvs:
    path_to_ref_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "data",
        "ukb_ref.csv",
    )
    df_ref = pd.read_csv(path_to_ref_csv, sep=";")

    df_clinical = df_clinical.set_index(
        ["participant_id", "sessions", "modality", "bids_filename"],
        verify_integrity=True,
    )

    return {
        col: df_clinical.filter(items=list(df_ref[col]))
        for col in ["participants", "sessions", "scans"]
    }


def write_bids(
    to: PathLike,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
    dataset_directory: PathLike,
) -> None:
    from pathlib import Path

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription
    from clinica.iotools.bids_utils import write_to_tsv

    to = Path(to)
    fs = LocalFileSystem(auto_mkdir=True)

    participants = participants.droplevel(
        ["sessions", "modality", "bids_filename"]
    ).drop_duplicates()

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(
            str(to / "dataset_description.json"), "w"
        ) as dataset_description_file:
            BIDSDatasetDescription(name="UKB").write(to=dataset_description_file)
        with fs.open(str(to / "participants.tsv"), "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby(["participant_id"]):
        sessions = data_frame.droplevel(
            ["participant_id", "modality", "bids_filename"]
        ).drop_duplicates()

        sessions_filepath = to / str(participant_id) / f"{participant_id}_sessions.tsv"
        with fs.open(str(sessions_filepath), "w") as sessions_file:
            write_to_tsv(sessions, sessions_file)
    scans = scans.set_index(["bids_full_path"], verify_integrity=True)
    for bids_full_path, metadata in scans.iterrows():
        if metadata["modality_num"] != "20217" and metadata["modality_num"] != "20225":
            copy_file_to_bids(
                zipfile=str(dataset_directory) + "/" + metadata["source_zipfile"],
                filenames=[metadata["source_filename"]] + metadata["sidecars"],
                bids_path=to / bids_full_path,
            )
        else:
            convert_dicom_to_nifti(
                zipfiles=str(dataset_directory) + "/" + metadata["source_zipfile"],
                bids_path=to / bids_full_path,
            )
            if metadata["modality_num"] == "20217":
                import_event_tsv(bids_path=str(to))
    return


def copy_file_to_bids(zipfile: str, filenames: List[str], bids_path: str) -> None:
    """Install the requested files in the BIDS  dataset."""
    import fsspec

    fo = fsspec.open(zipfile)
    fs = fsspec.filesystem("zip", fo=fo)
    for filename in filenames:
        if fs.exists(filename):
            bids_path_extension = str(bids_path) + "." + (filename.split(".", 1)[1])
            with fsspec.open(bids_path_extension, mode="wb") as f:
                f.write(fs.cat(filename))


def convert_dicom_to_nifti(zipfiles: str, bids_path: str) -> None:
    """Install the requested files in the BIDS  dataset.
    First, the dicom is extracted in a temporary directory
    Second, the dicom extracted is converted in the right place using dcm2niix"""
    import json
    import os
    import subprocess
    import tempfile
    import zipfile
    from pathlib import PurePath

    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)

    zf = zipfile.ZipFile(zipfiles)
    try:
        os.makedirs(PurePath(bids_path).parent)
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
        fmri_image_path = PurePath(find_largest_imaging_data(tempdir))
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
    return


def select_sessions(x: DataFrame) -> Optional[Series]:
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
            msg=f"The subject {x.eid} doesn't have the age for the imaging session number one (age_when_attended_assessment_centre_f21003_2_0)."
            f"It will not be converted. To have it converted, please update your clinical data.",
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
            msg=f"The subject {x.eid} doesn't have the age for the imaging session number two (age_when_attended_assessment_centre_f21003_3_0)."
            f"It will not be converted. To have it converted, please update your clinical data.",
            lvl="warning",
        )
        return None


def import_event_tsv(bids_path: str) -> None:
    """Import the csv containing the events' information."""
    import os

    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)
    path_to_event_tsv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
        "resources",
        "fmri",
        "task-facesshapesemotion_events.tsv",
    )

    bids_path_extension = str(bids_path) + "/" + "task-facesshapesemotion_events.tsv"
    fs.copy(path_to_event_tsv, bids_path_extension)
    return


def find_largest_imaging_data(path_to_source_data: PathLike) -> PathLike:
    import os
    from pathlib import Path

    weight = 0
    file_to_return = ""
    for z in Path(path_to_source_data).rglob("*.nii.gz"):
        if os.path.getsize(z) > weight:
            weight = os.path.getsize(z)
            file_to_return = z
    return Path(file_to_return)
