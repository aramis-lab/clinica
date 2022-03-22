from cmath import nan
from os import PathLike
from typing import BinaryIO, Iterable, List, Tuple, Union

from pandas import DataFrame, notna


def find_clinical_data(
    clinical_data_directory: PathLike,
) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    # Read the xls
    try:
        image_data_file = list(Path(clinical_data_directory).glob("*.tsv"))
    except StopIteration:
        raise FileNotFoundError("Clinical data file not found.")

    # Assign the dfs
    if len(image_data_file) == 1:
        df_clinical = pd.read_csv(str(image_data_file[0]), sep="\t")
        cprint(msg="All clinical data have been found", lvl="info")
    elif len(image_data_file) == 0:
        raise FileNotFoundError("Clinical data not found or incomplete")
    elif len(image_data_file) > 1:
        raise FileNotFoundError("Too many data files found")
    return df_clinical


def read_imaging_data(imaging_data_directory: PathLike) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    source_path_series = pd.Series(
        find_imaging_data(imaging_data_directory), name="source_path"
    )
    # list of the files we want to build the bids for each modality
    file_mod_list = [
        "T1_orig_defaced.nii.gz",
        "T2_FLAIR_orig_defaced.nii.gz",
        "AP.nii.gz",
        # "rfMRI.nii.gz",
        # "tfMRI.nii.gz",
    ]

    dataframe = pd.DataFrame.from_records(
        source_path_series, columns=["source_zipfile", "source_filename"]
    )
    filename = (
        dataframe["source_filename"]
        .apply(lambda x: Path(str(x)).name)
        .rename("filename")
    )
    dataframe = dataframe[filename.isin(file_mod_list)]
    filename = filename[filename.isin(file_mod_list)]
    split_zipfile = dataframe["source_zipfile"].str.split("_", expand=True)
    split_zipfile = split_zipfile.rename(
        {0: "source_id", 1: "modality_num", 2: "source_session_number"}, axis="columns"
    ).drop_duplicates()
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
    import pandas as pd

    df_clinical = df_clinical_data.merge(
        df_source, how="inner", right_on="source_id", left_on="eid"
    )
    return df_clinical


def complete_clinical(df_clinical: DataFrame) -> DataFrame:
    """This funbction uses the existing data to create the columns needed for the bids hierarchy (subject_id, ses, age_at _session, ect.)"""
    import pandas as pd

    df_clinical = df_clinical.assign(
        participant_id=lambda df: ("sub-UKB" + df.source_id.astype("str"))
    )
    df_clinical = df_clinical.assign(
        session=lambda df: "ses-" + df.source_session_number.astype("str")
    )
    df_clinical = df_clinical.join(
        df_clinical.modality_num.map(
            {
                "20253": {
                    "datatype": "anat",
                    "suffix": "FLAIR",
                    "sidecars": ["T2_FLAIR/T2_FLAIR.json"],
                },
                "20252": {
                    "datatype": "anat",
                    "suffix": "T1w",
                    "sidecars": ["T1/T1.json"],
                },
                "20250": {
                    "datatype": "dwi",
                    "suffix": "dwi",
                    "sidecars": [
                        "dMRI/raw/AP.json",
                        "dMRI/raw/AP.bval",
                        "dMRI/raw/AP.bvec",
                    ],
                },
                # "20227": {"datatype":"rfMRI", "suffix": "rfMRI","sidecars": []},
                # "20249": {"datatype": "tfMRI","suffix":"tfMRI", "sidecars": []},
                # "20252": {"datatype":"swi", "suffix": "swi", "sidecars":[]}
            }
        ).apply(pd.Series)
    )

    df_clinical = df_clinical.assign(year_of_birth=lambda df: df.year_of_birth_f34_0_0)

    df_clinical = df_clinical.assign(age=lambda df: df.age_at_recruitment_f21022_0_0)
    df_clinical = df_clinical.assign(
        age_at_first_session=lambda df: df.age_when_attended_assessment_centre_f21003_2_0
    )
    df_clinical["age_at_session"] = df_clinical.apply(
        lambda df: select_session(df), axis=1
    )
    df_clinical = df_clinical.assign(
        session_month=lambda df: (df.age_at_session - df.age_at_first_session) * 12
    )

    df_clinical = df_clinical.assign(
        session=lambda df: df.session_month.map(lambda x: f"ses-M{x:03d}")
    )
    df_clinical = df_clinical.assign(
        bids_filename=lambda df: (
            "sub-UKB"
            + df.source_id.astype("str")
            + "_"
            + df.session.astype("str")
            + "_"
            + df.suffix
            # + ".nii.gz"
        )
    )
    df_clinical = df_clinical.assign(
        sidecar_filename=lambda df: (
            "sub-UKB"
            + df.source_id.astype("str")
            + "_"
            + df.session.astype("str")
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
        + df.session
        + "/"
        + df.datatype
        + "/"
        + df.bids_filename
    )
    df_clinical = df_clinical.assign(
        bids_sidecar_path=lambda df: df.participant_id
        + "/"
        + df.session
        + "/"
        + df.datatype
        + "/"
        + df.sidecar_filename
    )
    print("clinical df id: ", df_clinical.source_id)
    print("clinical df session: ", df_clinical.source_session_number)
    return df_clinical


def dataset_to_bids(df_clinical: DataFrame) -> Tuple[DataFrame, DataFrame, DataFrame]:

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
        ["participant_id", "session", "suffix"], verify_integrity=True
    )

    # Build participants dataframe
    df_participants = df_clinical.filter(items=list(df_ref["participants"]))

    # Build sessions dataframe
    df_session = df_clinical.filter(items=list(df_ref["session"]))

    # Build scans dataframe
    df_scan = df_clinical.filter(items=list(df_ref["scan"]))
    return df_participants, df_session, df_scan


def write_to_tsv(dataframe: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    # Save dataframe as a BIDS-compliant TSV file.
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def write_bids(
    to: PathLike,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
    dataset_directory: PathLike,
) -> List[PathLike]:
    from pathlib import Path

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    to = Path(to)
    fs = LocalFileSystem(auto_mkdir=True)

    participants = participants.droplevel(["session", "suffix"]).drop_duplicates()

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name="UKB").write(to=dataset_description_file)

        with fs.open(to / "participants.tsv", "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby(["participant_id"]):
        session = data_frame.droplevel(["participant_id", "suffix"]).drop_duplicates()

        session_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
        print("session_filepath", session_filepath)
        with fs.open(session_filepath, "w") as sessions_file:
            write_to_tsv(session, sessions_file)

    scans = scans.set_index(["bids_full_path"], verify_integrity=True)
    for bids_full_path, metadata in scans.iterrows():
        copy_file_to_bids(
            zipfile=str(dataset_directory) + "/" + metadata["source_zipfile"],
            filenames=[metadata["source_filename"]] + metadata["sidecars"],
            bids_path=to / bids_full_path,
        )
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


def select_session(x):
    if x["source_session_number"] == "2":
        return x.age_when_attended_assessment_centre_f21003_2_0
    elif x["source_session_number"] == "3":
        return x.age_when_attended_assessment_centre_f21003_3_0
