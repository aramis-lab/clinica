from os import PathLike
from typing import BinaryIO, Iterable, List, Tuple, Union

from pandas import DataFrame


def find_clinical_data(
    clinical_data_directory: PathLike,
) -> Tuple[DataFrame, DataFrame, DataFrame, DataFrame, DataFrame]:
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    # Read the xls
    try:
        image_data_file = list(Path(clinical_data_directory).glob("*.csv"))
    except StopIteration:
        raise FileNotFoundError("Imaging collection file not found.")

    # Assign the dfs
    df_list = []
    for i in range(0, len(image_data_file)):
        df_list.append(pd.read_csv(str(image_data_file[i])))
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
    file_mod_list = ["T1_orig_defaced.nii.gz", "T2_FLAIR_orig_defaced.nii.gz"]

    subject_id = source_path_series.apply(lambda x: Path(str(x)).name).rename(
        "filename"
    )
    source_path_series = source_path_series[subject_id.isin(file_mod_list)]
    subject_id = subject_id[subject_id.isin(file_mod_list)]

    source_dir_series = (
        source_path_series.apply(lambda x: Path(Path(str(x)).parent).name)
        .rename("eid")
        .astype("int64")
    )

    df_source = pd.concat(
        [
            source_path_series,
            source_dir_series,
            subject_id,
        ],
        axis=1,
    )
    # df_source = df_source.rename(
    #     {0: "Subject", 1: "modality", 2: "Date"}, axis="columns"
    # ).drop_duplicates()

    # df_source = df_source.assign(participant_id=lambda df: "sub-" + df.Subject)
    return df_source


def find_imaging_data(path_to_source_data: PathLike) -> Iterable[PathLike]:
    from pathlib import Path

    for path in Path(path_to_source_data).rglob("*.nii.gz"):
        yield str(path.relative_to(path_to_source_data))


def intersect_data(df_source: DataFrame, dict_df: dict) -> Tuple[DataFrame, DataFrame]:
    import pandas as pd

    df_clinical = dict_df["clinical"]
    df_clinical_ev = df_clinical.merge(df_source["eid"], how="inner", on="eid")

    return df_source, df_clinical_ev


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

    # Build participants dataframe
    df_participants = df_clinical_ev.filter(items=list(df_ref["participants"]))
    name_dict = dict(
        zip(list(df_ref["participants"]), list(df_ref["participants_name"]))
    )
    type_dict = dict(
        zip(list(df_ref["participants_name"]), list(df_ref["participants_type"]))
    )
    df_participants = df_participants.rename(columns=name_dict)
    df_participants = df_participants.astype(type_dict)

    # Build sessions dataframe
    df_session = df_clinical_ev.filter(items=list(df_ref["session"]))
    name_dict = dict(zip(list(df_ref["session"]), list(df_ref["session_name"])))
    type_dict = dict(zip(list(df_ref["session_name"]), list(df_ref["session_type"])))
    df_session = df_session.rename(columns=name_dict)
    df_session = df_session.astype(type_dict)

    # Build scans dataframe
    df_scan = df_clinical_ev.filter(items=list(df_ref["scan"]))
    name_dict = dict(zip(list(df_ref["scan"]), list(df_ref["scan_name"])))
    type_dict = dict(zip(list(df_ref["scan_name"]), list(df_ref["scan_type"])))
    df_scan = df_scan.rename(columns=name_dict)
    df_scan = df_scan.astype(type_dict)

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

        for participant_id, sessions_group in sessions.groupby("participant_id"):
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(sessions_filepath, "w") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    # Perform import of imaging data next.
    for filename, metadata in scans.iterrows():
        path = Path(dataset_directory) / metadata.source_dir
        install_bids(sourcedata_dir=path, bids_filename=to / filename)
    return scans.index.to_list()
