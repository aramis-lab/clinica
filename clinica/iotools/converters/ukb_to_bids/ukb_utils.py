from cmath import nan
from os import PathLike
from typing import BinaryIO, Iterable, List, Tuple, Union

from pandas import DataFrame, notna


def find_clinical_data(
    clinical_data_directory: PathLike,
) -> Tuple[DataFrame, DataFrame, DataFrame, DataFrame, DataFrame]:
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    # Read the xls
    try:
        image_data_file = list(Path(clinical_data_directory).glob("*.tsv"))
    except StopIteration:
        raise FileNotFoundError("Clinical data file not found.")

    # Assign the dfs
    df_list = []
    for i in range(0, len(image_data_file)):
        df_list.append(pd.read_csv(str(image_data_file[i]), sep="\t"))
        df_list[i] = df_list[i].set_index(df_list[i].axes[1][0])
        if df_list[i].index.name == "eid":
            df_clinical = df_list[i]

    # Warn if a dataframe is missing
    if "df_clinical" in locals():
        cprint(msg="All clinical data have been found", lvl="info")
    else:
        raise FileNotFoundError("Clinical data not found or incomplete")

    return df_clinical


def read_clinical_data(
    clinical_data_directory: PathLike,
) -> Tuple[DataFrame, DataFrame, DataFrame, DataFrame, DataFrame]:
    df_clinical = find_clinical_data(clinical_data_directory)
    dictionnary = {
        "clinical": df_clinical,
    }
    return dictionnary


def read_imaging_data(imaging_data_directory: PathLike) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    source_path_series = pd.Series(
        find_imaging_data(imaging_data_directory), name="source_path"
    )
    # list of the files we want to build the bids for each mopdalkty
    file_mod_list = [
        "T1_orig_defaced.nii.gz",
        "T2_FLAIR_orig_defaced.nii.gz",
        "AP.nii.gz",
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
    df_data = dataframe["source_zipfile"].str.split("_", expand=True)
    df_data = df_data.rename(
        {0: "Subject", 1: "modality_alt", 2: "session_number"}, axis="columns"
    ).drop_duplicates()
    df_data = df_data.astype({"Subject": "int64"})

    df_source = pd.concat(
        [
            dataframe,
            filename,
            df_data,
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


def intersect_data(df_source: DataFrame, dict_df: dict) -> Tuple[DataFrame, DataFrame]:
    import pandas as pd

    df_clinical = dict_df["clinical"]
    df_clinical_ev = df_clinical.merge(
        df_source, how="inner", right_on="Subject", left_on="eid"
    )
    df_clinical_alt = df_clinical_ev.assign(
        participant_id=lambda df: ("sub-" + df.Subject.astype("str"))
    )
    df_clinical_alt = df_clinical_alt.assign(
        session=lambda df: "ses-" + df.session_number.astype("str")
    )
    df_clinical_alt = df_clinical_alt.join(
        df_clinical_alt.modality_alt.map(
            {
                "20253": {
                    "datatype": "anat",
                    "suffix": "FLAIR",
                    "json": "T2_FLAIR/T2_FLAIR.json",
                    "sidecars": [],
                },
                "20252": {
                    "datatype": "anat",
                    "suffix": "T1w",
                    "json": "T1/T1.json",
                    "sidecars": [],
                },
                "20250": {
                    "datatype": "dwi",
                    "suffix": "dwi",
                    "json": "dMRI/raw/AP.json",
                    "sidecars": ["dMRI/raw/AP.bval", "dMRI/raw/AP.bvec"],
                },
            }
        ).apply(pd.Series)
    )

    df_clinical_alt = df_clinical_alt.assign(
        year_of_birth=lambda df: df.year_of_birth_f34_0_0
    )

    df_clinical_alt = df_clinical_alt.assign(
        age=lambda df: df.age_at_recruitment_f21022_0_0
    )
    df_clinical_alt["age_at_session"] = df_clinical_alt.apply(
        lambda df: select_session(df), axis=1
    )
    df_clinical_alt = df_clinical_alt.assign(
        session_month=lambda df: (df.age_at_session - df.age) * 12
    )

    df_clinical_alt = df_clinical_alt.assign(
        session=lambda df: df.session_month.map(lambda x: f"ses-M{x}")
    )
    df_clinical_alt = df_clinical_alt.assign(
        bids_filename=lambda df: (
            "sub-UKB"
            + df.Subject.astype("str")
            + "_"
            + df.session.astype("str")
            + "_"
            + df.suffix
            + ".nii.gz"
        )
    )
    df_clinical_alt = df_clinical_alt.assign(
        sidecar_filename=lambda df: (
            "sub-UKB"
            + df.Subject.astype("str")
            + "_"
            + df.session.astype("str")
            + "_"
            + df.suffix
        )
    )

    print("df_clinical alt sessions: ", df_clinical_alt.session, "\n")
    df_clinical_alt = df_clinical_alt.join(
        df_clinical_alt.sex_f31_0_0.astype("str")
        .map(
            {
                "0": {"sex": "F"},
                "1": {"sex": "M"},
            }
        )
        .apply(pd.Series)
    )
    df_clinical_alt = df_clinical_alt.assign(
        bids_full_path=lambda df: df.participant_id
        + "/"
        + df.session
        + "/"
        + df.datatype
        + "/"
        + df.bids_filename
    )
    df_clinical_alt = df_clinical_alt.assign(
        bids_sidecar_path=lambda df: df.participant_id
        + "/"
        + df.session
        + "/"
        + df.datatype
        + "/"
        + df.sidecar_filename
    )
    return df_source, df_clinical_alt


def dataset_to_bids(
    df_source: DataFrame, df_clinical_ev: DataFrame
) -> Tuple[DataFrame, DataFrame, DataFrame]:

    import os

    import pandas as pd

    # open the reference for building the tsvs:
    path_to_ref_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "data",
        "ukb_ref.csv",
    )
    df_ref = pd.read_csv(path_to_ref_csv, sep=";")

    df_clinical_ev = df_clinical_ev.set_index(
        ["participant_id", "session", "suffix"], verify_integrity=True
    )

    # Build participants dataframe
    df_participants = df_clinical_ev.filter(items=list(df_ref["participants"]))

    # Build sessions dataframe
    df_session = df_clinical_ev.filter(items=list(df_ref["session"]))

    # Build scans dataframe
    df_scan = df_clinical_ev.filter(items=list(df_ref["scan"]))
    return df_participants, df_session, df_scan


def write_to_tsv(dataframe: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    # Save dataframe as a BIDS-compliant TSV file.
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def install_bids(sourcedata_dir: PathLike, bids_filename: PathLike) -> None:
    from pathlib import Path

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
    sidecar_dir = sourcedata_dir.parent
    for source_sidecar in sidecar_dir.rglob(f"{source_basename}*"):
        target_sidecar = Path.joinpath(bids_filename.parent, target_basename).with_name(
            f"{target_basename}{source_sidecar.suffix}"
        )
        source_file = fs.open(source_sidecar, mode="rb")
        target_file = fs.open(target_sidecar, mode="wb")
        with source_file as sf, target_file as tf:
            tf.write(sf.read())


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

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name="UKB").write(to=dataset_description_file)

        with fs.open(to / "participants.tsv", "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby(["participant_id"]):
        session = data_frame.droplevel("participant_id")
        session_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
        with fs.open(session_filepath, "w") as sessions_file:
            write_to_tsv(session, sessions_file)

    scans = scans.set_index(["bids_full_path"], verify_integrity=True)
    for bids_full_path, metadata in scans.iterrows():
        install_nifti(
            zipfile=str(dataset_directory) + "/" + metadata["source_zipfile"],
            filename=metadata["source_filename"],
            bids_path=to / bids_full_path,
        )
        install_json(
            zipfile=str(dataset_directory) + "/" + metadata["source_zipfile"],
            filename=metadata["json"],
            bids_path=to / metadata["bids_sidecar_path"],
        )
        install_sidecars(
            zipfile=str(dataset_directory) + "/" + metadata["source_zipfile"],
            filename=metadata["sidecars"],
            bids_path=to / metadata["bids_sidecar_path"],
        )
    return


def install_nifti(zipfile: str, filename: str, bids_path: str) -> None:
    """Install a NIfTI file from a source archive to the target BIDS path."""
    import fsspec

    fo = fsspec.open(zipfile)
    fs = fsspec.filesystem("zip", fo=fo)
    with fsspec.open(bids_path, mode="wb") as f:
        f.write(fs.cat(filename))


def install_json(zipfile: str, filename: str, bids_path: str) -> None:
    """Install a NIfTI file from a source archive to the target BIDS path."""
    import fsspec

    fo = fsspec.open(zipfile)
    fs = fsspec.filesystem("zip", fo=fo)
    bids_path_json = str(bids_path) + ".json"
    if fs.exists(filename):
        with fsspec.open(bids_path_json, mode="wb") as f:
            f.write(fs.cat(filename))


def install_sidecars(zipfile: str, filename: str, bids_path: str) -> None:
    """Install a NIfTI file from a source archive to the target BIDS path."""
    import fsspec

    fo = fsspec.open(zipfile)
    fs = fsspec.filesystem("zip", fo=fo)
    for i in range(0, len(filename)):
        print("metadata !:", filename[i])
        if fs.exists(filename[i]):
            bids_path_extension = str(bids_path) + "." + (filename[i].split(".")[1])
            print("bids_path_extension: ", bids_path_extension)
            with fsspec.open(bids_path_extension, mode="wb") as f:
                f.write(fs.cat(filename[i]))


def select_session(x):
    if x["session_number"] == "2":
        return x.age_when_attended_assessment_centre_f21003_2_0
    elif x["session_number"] == "3":
        return x.age_when_attended_assessment_centre_f21003_3_0
