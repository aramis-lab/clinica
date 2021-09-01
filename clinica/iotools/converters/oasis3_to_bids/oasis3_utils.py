from os import PathLike
from typing import BinaryIO, Iterable, List, Optional, Tuple, Union

from pandas import DataFrame


def find_clinical_data(clinical_data_directory: PathLike) -> Optional[DataFrame]:
    from pathlib import Path

    import pandas as pd

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
        print("All the clinical data have been found.")
    else:
        raise FileNotFoundError("Clinical data not found or incomplete")

    return df_pet, df_mri, df_subject, df_pup, df_adrc


def read_clinical_data(clinical_data_directory: PathLike) -> DataFrame:
    from pandas import NA, CategoricalDtype

    df_pet, df_mri, df_subject, df_pup, df_adrc = find_clinical_data(
        clinical_data_directory
    )
    return df_pet, df_mri, df_subject, df_pup, df_adrc


def read_imaging_data(imaging_data_directory: PathLike) -> Iterable[Tuple[str, str]]:
    import pandas as pd

    s = pd.Series(find_imaging_data(imaging_data_directory), name="source_path")

    a = s.str.split(pat="/", n=1, expand=True).drop(columns=1)[0].rename("path")
    df_source = pd.concat([s, a, a.str.split("_", expand=True)], axis=1)
    df_source = df_source.rename(
        {0: "Subject", 1: "modality", 2: "Date"}, axis="columns"
    ).drop_duplicates()

    return df_source


def find_imaging_data(path_to_source_data):
    from pathlib import Path

    for path in Path(path_to_source_data).rglob("*.nii.gz"):
        yield str(path.relative_to(path_to_source_data))


def intersect_data(df_source, df_mri, df_subject, df_adrc):
    import pandas as pd

    df_adrc = df_adrc.assign(
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
        age=lambda df: (df["ageAtEntry"]) + df["Date"].str[1:].astype("float") / 365
    )

    df_subject_small = df_subject.merge(
        df_source["Subject"], how="inner", on="Subject"
    ).drop_duplicates()
    df_source = df_source.merge(
        df_mri["Scanner"], how="left", left_on="path", right_on="MR ID"
    )

    df_small = df_subject_small.merge(df_adrc_small, how="inner", on="Subject")

    df_source = df_source.assign(
        session=lambda df: (round(df["Date"].str[1:].astype("int") / (364.25 / 2)) * 6)
    )

    df_source = df_source.assign(
        ses=lambda df: df.session.apply(
            lambda x: f"ses-M{ (5-(len(str(x)))) * '0' + str(int(x))}"
        )
    )

    df_source = df_source.join(
        df_source.modality.map(
            {
                "MR": {"datatype": "anat", "suffix": "T1w", "trc_label": ""},
                "FDG": {"datatype": "pet", "suffix": "pet", "trc_label": "acq-fdg_"},
                "PIB": {"datatype": "pet", "suffix": "pet", "trc_label": "acq-pib_"},
                "AV45": {"datatype": "pet", "suffix": "pet", "trc_label": "acq-av45_"},
            }
        ).apply(pd.Series)
    )

    df_source = df_source.assign(
        Output_filename=lambda df: df.datatype
        + "/sub-"
        + df.Subject
        + "_"
        + df.ses
        + "_"
        + df.trc_label
        + df.suffix
        + ".nii.gz"
    )
    return df_source, df_small


def dataset_to_bids(df_source, df_small):
    import pandas as pd

    # build participants .tsv (subjects)
    df_participants = pd.DataFrame(
        {
            "participant_id": df_small["Subject"],
            "age": df_small["ageAtEntry"],
            "sex": df_small["M/F"],
            "handedness": df_small["Hand"],
        }
    )
    # build sessions .tsv
    df_session = pd.DataFrame(
        {
            "session_id": df_source.ses,
            "date_from_bl": df_source["Date"],
            "age": df_source["age"].astype("int"),
        }
    )
    # df_session = df_session.merge(df_adrc)
    df_scan = pd.DataFrame(
        {"filename": df_source.Output_filename, "source_dir": df_source.source_path}
    )
    return df_participants, df_session, df_scan


def write_to_tsv(dataframe: DataFrame, buffer: Union[PathLike, BinaryIO]) -> None:
    # Save dataframe as a BIDS-compliant TSV file.
    dataframe.to_csv(buffer, sep="\t", na_rep="n/a", date_format="%Y-%m-%d")


def install_nifti(sourcedata_dir: PathLike, bids_filename: PathLike) -> None:
    from fsspec.implementations.local import LocalFileSystem

    fs = LocalFileSystem(auto_mkdir=True)
    source_file = fs.open(fs.ls(sourcedata_dir)[0], mode="rb")
    target_file = fs.open(bids_filename, mode="wb")

    with source_file as sf, target_file as tf:
        tf.write(sf.read())


def write_bids(
    to: PathLike, participants: DataFrame, sessions: DataFrame, scans: DataFrame
) -> List[PathLike]:
    from pathlib import Path

    from fsspec.implementations.local import LocalFileSystem

    to = Path(to)
    fs = LocalFileSystem(auto_mkdir=True)

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(to / "participants.tsv", "wb") as participant_file:
            write_to_tsv(participants, participant_file)

        for participant_id, sessions_group in sessions.groupby("participant_id"):
            sessions_group = sessions_group.droplevel("participant_id")
            sessions_filepath = to / participant_id / f"{participant_id}_sessions.tsv"
            with fs.open(sessions_filepath, "wb") as sessions_file:
                write_to_tsv(sessions_group, sessions_file)

    # Perform import of imaging data next.
    for filename, metadata in scans.iterrows():
        install_nifti(sourcedata_dir=metadata.source_dir, bids_filename=to / filename)
    return scans.index.to_list()
