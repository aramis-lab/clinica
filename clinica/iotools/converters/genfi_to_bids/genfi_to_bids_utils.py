from cmath import nan
from os import PathLike
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
import pydicom as pdcm
from pandas import DataFrame, Series


def find_dicoms(path_to_source_data):

    for z in Path(path_to_source_data).rglob("*.dcm"):
        yield str(z), str(Path(z).parent)


def filter_dicoms(df_dicom):
    not_genfi_list = [
        "t2_2d_axial",
        "dwi_trace",
        "dwi_adc",
        "dwi_fa",
        "localizer",
        "localiser",
        "nonimage",
        "dual_pd_t2_axial",
        "t1_1.25mm_iso",
        "MoCo",
    ]
    not_genfi2_list = [
        "ASL",
        "motioncorrected",
        "localiser",
        "localizer",
    ]

    df_dicom = df_dicom.drop_duplicates(subset=["source"])  # .set_index("source_path")
    df_dicom = df_dicom.assign(
        series_desc=lambda df: df.source_path.apply(
            lambda x: pdcm.dcmread(x).SeriesDescription
        )
    )
    df_dicom = df_dicom.assign(
        acq_date=lambda df: df.source_path.apply(lambda x: pdcm.dcmread(x).StudyDate)
    )
    df_dicom = df_dicom.set_index(["source_path"], verify_integrity=True)

    df_dicom = df_dicom[~df_dicom["source"].str.contains("secondary", case=False)]
    for file_mod in not_genfi2_list:
        df_dicom = df_dicom[~df_dicom["series_desc"].str.contains(file_mod, case=False)]
    for file_mod in not_genfi_list:
        df_dicom = df_dicom[~df_dicom["series_desc"].str.contains(file_mod, case=False)]
    return df_dicom


def identify_modality(x):
    x = x.lower()
    if "dwi" in x:
        return "dwi"
    if "t1" in x:
        return "T1"
    if "t2" in x:
        return "T2w"
    if "fieldmap" in x:
        return "fieldmap"
    if "fmri" in x:
        return "rsfmri"
    else:
        return None


def convert_dicom_to_nifti(in_file: str, bids_path: str, bids_filename: str):
    import subprocess

    command = [
        "dcm2niix",
        "-w",
        "0",
    ]
    command += ["-9", "-z", "y"]
    command += ["-b", "y", "-ba", "y"]
    command += ["-o", bids_path, "-f", bids_filename]
    command += [in_file]
    subprocess.run(command, capture_output=True)
    return


def find_clinical_data():
    return


def complete_clinical():
    return


def dataset_to_bids(complete_data_df):
    # generates participants, sessions and scans tsv
    complete_data_df = complete_data_df.set_index(
        ["participant_id", "session_id", "modality", "bids_filename"],
        verify_integrity=True,
    )
    participants = complete_data_df.filter(items=["participant_id", "source"])
    sessions = complete_data_df.filter(
        items=["participant_id", "session_id", "genfi_version"]
    )
    scans = complete_data_df
    return participants, sessions, scans


def intersect_data():
    return


def read_imaging_data(source_path):
    df_dicom = pd.DataFrame(find_dicoms(source_path), columns=["source_path", "source"])
    df_dicom = filter_dicoms(df_dicom)
    return df_dicom


def complete_data(df_dicom):
    df_dicom = df_dicom.assign(
        source_id=lambda df: df.source.apply(
            lambda x: Path(Path(Path(x).parent).parent).name.split("-")[0]
        )
    )
    df_dicom = df_dicom.assign(
        source_ses_id=lambda df: df.source.apply(
            lambda x: Path(Path(Path(x).parent).parent).name.split("-")[1]
        )
    )
    df_dicom = df_dicom.assign(
        origin=lambda df: df.source.apply(
            lambda x: Path(Path(Path(Path(Path(x).parent).parent).parent).parent)
        )
    )

    df_dicom = df_dicom.reset_index()

    df_dicom = df_dicom.assign(source_ses_id=lambda df: df.source_ses_id.astype("int"))
    df_dicom = df_dicom.assign(
        genfi_version=lambda df: df.source_ses_id.apply(lambda x: f"GENFI{len(str(x))}")
    )
    df_1 = (
        df_dicom[["source_id", "source_ses_id", "acq_date"]]
        .groupby(["source_id", "source_ses_id"])
        .min()
    )
    df_2 = df_dicom[["source_id", "acq_date"]].groupby(["source_id"]).min()
    df_1 = df_1.join(df_2.rename(columns={"acq_date": "baseline"}))
    df_1 = df_1.assign(
        ses_month=lambda df: (
            (
                df["acq_date"].str[:4].astype("int")
                - df["baseline"].str[:4].astype("int")
            )
            * 12
            + (
                df["acq_date"].str[4:6].astype("int")
                - df["baseline"].str[4:6].astype("int")
            )
        )
    )
    df_1 = df_1.assign(
        session_id=lambda df: df.ses_month.map(lambda x: f"ses-M{x:03d}")
    )

    df_sub_ses = df_dicom.merge(
        df_1[["baseline", "session_id"]], how="inner", on=["source_id", "source_ses_id"]
    )
    df_sub_ses = df_sub_ses.assign(
        participant_id=lambda df: df.source_id.apply(lambda x: f"sub-{x}")
    )

    df_sub_ses = df_sub_ses.assign(
        modality=lambda df: df.series_desc.apply(lambda x: identify_modality(x))
    )
    df_sub_ses = df_sub_ses.join(
        df_sub_ses.modality.map(
            {
                "T2w": {
                    "datatype": "anat",
                    "suffix": "T2w",
                    "sidecars": ["T2w.json"],
                    "task": "",
                },
                "T1": {
                    "datatype": "anat",
                    "suffix": "T1w",
                    "sidecars": ["T1w.json"],
                    "task": "",
                },
                "dwi": {
                    "datatype": "dwi",
                    "suffix": "dwi",
                    "sidecars": [
                        "dwi.json",
                        "dwi.bval",
                        "dwi.bvec",
                    ],
                    "task": "",
                },
                "rsfmri": {
                    "datatype": "func",
                    "suffix": "bold",
                    "sidecars": ["rsfmri.json"],
                    "task": "_task-rest",
                },
                "fieldmap": {
                    "datatype": "fmap",
                    "suffix": "fmap",
                    "sidecars": ["fmap.json"],
                    "task": "",
                },
            }
        ).apply(pd.Series)
    )

    # take into account fieldmaps -> same method as for baseline (extract number of directory)
    df_sub_ses_fmap = df_sub_ses.assign(
        dir_num=lambda df: df.source.apply(lambda x: int(Path(Path(x).parent).name))
    )
    df_map_1 = (
        df_sub_ses_fmap[
            ["source_id", "source_ses_id", "modality", "dir_num", "suffix"]
        ][df_sub_ses["modality"].str.contains("fieldmap")]
        .groupby(["source_id", "source_ses_id", "modality", "dir_num"])
        .min()
    )
    df_map_2 = (
        df_sub_ses_fmap[["source_id", "source_ses_id", "modality", "dir_num"]][
            df_sub_ses["modality"].str.contains("fieldmap")
        ]
        .groupby(["source_id", "source_ses_id", "modality"])
        .min()
    )
    df_map_1 = df_map_1.join(df_map_2.rename(columns={"dir_num": "run_01_dir_num"}))
    df_alt_map = df_map_1.reset_index().assign(
        fmap_type_p=lambda df: (df.run_01_dir_num != df.dir_num)
    )
    df_alt_map = df_alt_map.assign(
        fmap=lambda df: df.fmap_type_p.apply(
            lambda x: "phasediff" if x else "magnitude"
        )
    )

    # merge fieldmap suffix
    df_sub_ses_map_1 = df_sub_ses_fmap.merge(
        df_alt_map[["source_id", "source_ses_id", "modality", "dir_num", "fmap"]],
        how="inner",
        on=["source_id", "source_ses_id", "modality", "dir_num"],
    )
    df_sub_ses_map_2 = df_sub_ses_fmap.merge(
        df_alt_map[["source_id", "source_ses_id", "modality", "dir_num", "fmap"]],
        how="left",
        on=["source_id", "source_ses_id", "modality", "dir_num"],
    )
    df_sub_ses_map_2 = df_sub_ses_map_2[
        ~df_sub_ses_map_2["modality"].str.contains("fieldmap")
    ]
    df_sub_ses_map_1 = df_sub_ses_map_1.assign(suffix=lambda df: df.fmap)
    df_suf = pd.concat([df_sub_ses_map_2, df_sub_ses_map_1], ignore_index=True)

    # take into account runs -> same method as for baseline (extract number of directory)
    df_suf_dir = df_suf.assign(
        dir_num=lambda df: df.source.apply(lambda x: int(Path(Path(x).parent).name))
    )
    df_1 = (
        df_suf_dir[["source_id", "source_ses_id", "suffix", "dir_num"]]
        .groupby(["source_id", "source_ses_id", "suffix", "dir_num"])
        .min()
    )
    df_2 = (
        df_suf_dir[["source_id", "source_ses_id", "suffix", "dir_num"]]
        .groupby(["source_id", "source_ses_id", "suffix"])
        .min()
    )
    df_1 = df_1.join(df_2.rename(columns={"dir_num": "run_01_dir_num"}))
    df_alt = df_1.reset_index().assign(run=lambda df: (df.run_01_dir_num != df.dir_num))
    df_alt = df_alt.assign(
        run_num=lambda df: df.run.apply(lambda x: f"run-0{int(x)+1}")
    )

    df_sub_ses_run = df_suf_dir.merge(
        df_alt[["source_id", "source_ses_id", "suffix", "dir_num", "run_num"]],
        how="left",
        on=["source_id", "source_ses_id", "suffix", "dir_num"],
    )

    ##builds path to bids
    df_sub_ses_run = df_sub_ses_run.assign(
        bids_filename=lambda df: (
            df.participant_id + "_" + df.session_id + "_" + df.run_num + "_" + df.suffix
        )
    )
    df_sub_ses_run = df_sub_ses_run.assign(
        bids_full_path=lambda df: (
            df.participant_id
            + "/"
            + df.session_id
            + "/"
            + df.datatype
            + "/"
            + df.bids_filename
        )
    )

    return df_sub_ses_run


def write_bids(
    to: PathLike,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
    dataset_directory: PathLike,
) -> None:
    import os
    from pathlib import Path

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription
    from clinica.iotools.bids_utils import write_to_tsv

    to = Path(to)
    fs = LocalFileSystem(auto_mkdir=True)

    # Ensure BIDS hierarchy is written first.
    with fs.transaction:
        with fs.open(
            str(to / "dataset_description.json"), "w"
        ) as dataset_description_file:
            BIDSDatasetDescription(name="GENFI").write(to=dataset_description_file)
        with fs.open(str(to / "participants.tsv"), "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby(["participant_id"]):
        sessions = data_frame.droplevel(
            ["participant_id", "modality", "bids_filename"]
        ).drop_duplicates()

        sessions_filepath = to / str(participant_id) / f"{participant_id}_sessions.tsv"
        with fs.open(str(sessions_filepath), "w") as sessions_file:
            write_to_tsv(sessions, sessions_file)
    scans = scans.reset_index().set_index(["bids_full_path"], verify_integrity=True)

    for bids_full_path, metadata in scans.iterrows():
        try:
            os.makedirs(to / (Path(bids_full_path).parent))
        except OSError:
            pass
        convert_dicom_to_nifti(
            Path(metadata["source_path"]).parent,
            to / str(Path(bids_full_path).parent),
            metadata["bids_filename"],
        )

    return
