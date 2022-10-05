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
        if df_list[i].index.name == "XNAT_PETSESSIONDATA ID":
            df_pet = df_list[i]
        if df_list[i].index.name == "MR ID":
            df_mri = df_list[i]
        if df_list[i].index.name == "Subject":
            df_subject = df_list[i]
        if df_list[i].index.name == "PUP_PUPTIMECOURSEDATA ID":
            df_pup = df_list[i]
        if df_list[i].index.name == "ADRC_ADRCCLINICALDATA ID":
            df_adrc = df_list[i]

    # Warn if a dataframe is missing
    if (
        "df_pup" in locals()
        and "df_mri" in locals()
        and "df_adrc" in locals()
        and "df_pet" in locals()
    ):
        cprint(msg="All clinical data have been found", lvl="info")
    else:
        raise FileNotFoundError("Clinical data not found or incomplete")

    return df_pet, df_mri, df_subject, df_pup, df_adrc


def read_clinical_data(
    clinical_data_directory: PathLike,
) -> Tuple[DataFrame, DataFrame, DataFrame, DataFrame, DataFrame]:
    df_pet, df_mri, df_subject, df_pup, df_adrc = find_clinical_data(
        clinical_data_directory
    )
    dictionnary = {
        "pet": df_pet,
        "mri": df_mri,
        "subject": df_subject,
        "pup": df_pup,
        "adrc": df_adrc,
    }
    return dictionnary


def read_imaging_data(imaging_data_directory: PathLike) -> DataFrame:
    from pathlib import Path

    import pandas as pd

    source_path_series = pd.Series(
        find_imaging_data(imaging_data_directory), name="source_path"
    )

    source_dir_series = source_path_series.apply(lambda x: Path(str(x)).parent).rename(
        "source_dir"
    )

    file_spec_series = source_path_series.apply(lambda x: Path(str(x)).parts[0]).rename(
        "path"
    )

    df_source = pd.concat(
        [
            source_path_series,
            file_spec_series,
            source_dir_series,
            file_spec_series.str.split("_", expand=True),
        ],
        axis=1,
    )
    df_source = df_source.rename(
        {0: "Subject", 1: "modality", 2: "Date"}, axis="columns"
    ).drop_duplicates()

    df_source = df_source.assign(participant_id=lambda df: "sub-" + df.Subject)
    return df_source


def find_imaging_data(path_to_source_data: PathLike) -> Iterable[PathLike]:
    from pathlib import Path

    for path in Path(path_to_source_data).rglob("*.nii.gz"):
        yield str(path.relative_to(path_to_source_data))


def intersect_data(df_source: DataFrame, dict_df: dict) -> Tuple[DataFrame, DataFrame]:
    import pandas as pd

    df_adrc = dict_df["adrc"].assign(
        session_id=lambda df: df.index.map(lambda x: x.split("_")[2])
    )
    df_adrc_small = df_adrc[df_adrc.session_id == "d0000"]
    df_adrc_small = df_adrc_small.merge(
        df_source["Subject"], how="inner", on="Subject"
    ).drop_duplicates()
    df_source = df_source.merge(
        df_adrc_small[["Subject", "ageAtEntry"]], how="inner", on="Subject"
    )
    df_source = df_source.assign(
        age=lambda df: (df["ageAtEntry"]) + df["Date"].str[1:].astype("float") / 365.25
    )

    df_subject_small = (
        dict_df["subject"]
        .merge(df_source["Subject"], how="inner", on="Subject")
        .drop_duplicates()
    )
    df_source = df_source.merge(
        dict_df["mri"]["Scanner"], how="left", left_on="path", right_on="MR ID"
    )

    df_small = df_subject_small.merge(df_adrc_small, how="inner", on="Subject")
    df_small = df_small.assign(participant_id=lambda df: "sub-" + df.Subject)
    df_source = df_source.assign(
        session=lambda df: (round(df["Date"].str[1:].astype("int") / (365.25 / 2)) * 6)
    )

    df_source = df_source.assign(
        ses=lambda df: df.session.apply(
            lambda x: f"ses-M{ (5-(len(str(x)))) * '0' + str(int(x))}"
        )
    )

    df_source = df_source.join(
        df_source.modality.map(
            {
                "MR": {"datatype": "anat", "suffix": "T1w"},
                "FDG": {"datatype": "pet", "suffix": "pet", "trc_label": "18FFDG"},
                "PIB": {"datatype": "pet", "suffix": "pet", "trc_label": "11CPIB"},
                "AV45": {"datatype": "pet", "suffix": "pet", "trc_label": "18FAV45"},
            }
        ).apply(pd.Series)
    )
    if "trc_label" in df_source.columns:
        df_source = df_source.assign(
            filename=lambda df: df.apply(
                lambda x: f"{x.participant_id}/{x.ses}/{x.datatype}/"
                f"{x.participant_id}_{x.ses}"
                f"{'_trc-'+x.trc_label if pd.notna(x.trc_label) else ''}"
                f"_{x.suffix}.nii.gz",
                axis=1,
            )
        )
    else:
        df_source = df_source.assign(
            filename=lambda df: df.apply(
                lambda x: f"{x.participant_id}/{x.ses}/{x.datatype}/"
                f"{x.participant_id}_{x.ses}"
                f"_{x.suffix}.nii.gz",
                axis=1,
            )
        )
    df_adrc = df_adrc.merge(df_source["Subject"], how="inner", on="Subject")
    df_adrc = df_adrc.assign(
        session=lambda df: round(df["session_id"].str[1:].astype("int") / (364.25 / 2))
        * 6
    )
    df_adrc = df_adrc.drop_duplicates().set_index(["Subject", "session_id"])
    df_source = df_source.merge(
        df_adrc[
            [
                "session",
                "mmse",
                "cdr",
                "commun",
                "dx1",
                "homehobb",
                "judgment",
                "memory",
                "orient",
                "perscare",
                "sumbox",
                "apoe",
            ]
        ],
        how="left",
        on="session",
    )
    return df_source, df_small


def dataset_to_bids(
    df_source: DataFrame, df_small: DataFrame
) -> Tuple[DataFrame, DataFrame, DataFrame]:
    # Build participants dataframe
    df_participants = (
        df_small[["participant_id", "ageAtEntry", "M/F", "Hand", "Education"]]
        .rename(
            columns={
                "ageAtEntry": "age",
                "M/F": "sex",
                "Hand": "handedness",
            }
        )
        .set_index("participant_id", verify_integrity=True)
    )

    # Build sessions dataframe
    df_session = (
        df_source[
            [
                "participant_id",
                "ses",
                "Date",
                "age",
                "mmse",
                "cdr",
                "commun",
                "dx1",
                "homehobb",
                "judgment",
                "memory",
                "orient",
                "perscare",
                "sumbox",
                "apoe",
            ]
        ]
        .rename(
            columns={
                "ses": "session_id",
                "Date": "source_session_id",
            }
        )
        .astype({"age": "int"})
        .drop_duplicates(subset=["participant_id", "session_id"])
        .set_index(["participant_id", "session_id"], verify_integrity=True)
    )

    # Build scans dataframe
    df_scan = (
        df_source[["filename", "source_dir"]]
        .drop_duplicates(subset=["filename"])
        .set_index("filename", verify_integrity=True)
    )

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
    sidecar_dir = sourcedata_dir.parent / "BIDS"
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

    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name="OASIS-3").write(to=dataset_description_file)
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
