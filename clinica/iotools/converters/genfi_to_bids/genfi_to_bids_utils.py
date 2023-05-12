from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import pydicom as pdcm
from pandas import DataFrame


def find_dicoms(path_to_source_data: PathLike) -> Iterable[Tuple[PathLike, PathLike]]:
    """Find the dicoms in the given directory.

    Parameters
    ----------
    path_to_source_data: PathLike
        Path to the source data

    Returns
    -------
    Iterable[Tuple[PathLike, PathLike]]
        Path to found files and parent directory
    """
    from clinica.utils.stream import cprint

    cprint("Looking for imaging data.", lvl="info")
    for z in Path(path_to_source_data).rglob("*.dcm"):
        yield str(z), str(Path(z).parent)


def handle_manufacturer(x: str) -> str:
    try:
        manufacturer = pdcm.dcmread(x).Manufacturer
    except:
        manufacturer = "Unknown"
    return manufacturer


def filter_dicoms(df: DataFrame) -> DataFrame:
    """Filters modalities handled by the converter.

    Parameters
    ----------
    df: DataFrame
        Dataframe containing all of the dicoms available in the source directory.

    Returns
    -------
    df: DataFrame
        Dataframe with only the modalities handled.
    """
    to_filter = [
        # from GENFI 1
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
        # from GENFI 2
        "ASL",
        "motioncorrected",
        "localiser",
        "localizer",
    ]

    df = df.drop_duplicates(subset=["source"])
    df = df.assign(
        series_desc=lambda df: df.source_path.apply(
            lambda x: pdcm.dcmread(x).SeriesDescription
        ),
        acq_date=lambda df: df.source_path.apply(lambda x: pdcm.dcmread(x).StudyDate),
        manufacturer=lambda df: df.source_path.apply(lambda x: handle_manufacturer(x)),
    )
    df = df.set_index(["source_path"], verify_integrity=True)

    df = df[~df["source"].str.contains("secondary", case=False)]
    for file_mod in to_filter:
        df = df[~df["series_desc"].str.contains(file_mod, case=False)]
    return df


def _check_file(directory: PathLike, pattern: str) -> PathLike:
    from pathlib import Path

    try:
        data_file = list(Path(directory).glob(pattern))
    except StopIteration:
        raise FileNotFoundError("Clinical data file not found.")
    if len(data_file) == 0:
        raise FileNotFoundError("Clinical data not found or incomplete. Aborting")
    if len(data_file) > 1:
        raise ValueError("Too many data files found, expected one. Aborting.")
    return data_file[0]


def _read_file(data_file: PathLike) -> pd.DataFrame:

    return (
        pd.concat(
            [
                pd.read_excel(str(data_file)),
                pd.read_excel(str(data_file), sheet_name=1),
            ]
        )
        .convert_dtypes()
        .rename(columns=lambda x: x.lower().replace(" ", "_"))
    )


def find_clinical_data(
    clinical_data_directory: PathLike,
) -> List[DataFrame]:
    """Finds the clinical data associated with the dataset.

    Parameters
    ----------
    clinical_data_directory: PathLike
        Path to the clinical data.

    Returns
    -------
    List[DataFrame]
        Dataframes containing the clinical data
    """
    from clinica.utils.stream import cprint

    cprint("Looking for clinical data.", lvl="info")

    return (
        _read_file(_check_file(clinical_data_directory, pattern))
        for pattern in (
            "FINAL*DEMOGRAPHICS*.xlsx",
            "FINAL*IMAGING*.xlsx",
            "FINAL*CLINICAL*.xlsx",
        )
    )


def complete_clinical_data(
    df_demographics: DataFrame, df_imaging: DataFrame, df_clinical: DataFrame
) -> DataFrame:
    """Merges the different clincal dataframes into one.

    Parameters
    ----------
    df_demographics: DataFrame
        Dataframe containing the demographic data

    df_imaging: DataFrame
        Dataframe containing the imaging data

    df_clinical: DataFrame
        Dataframe containing the clinical data

    Returns
    -------
    df_clinical_complete: DataFrame
        Dataframe with the data of the 3 input dataframes
    """
    merge_key = ["blinded_code", "blinded_site", "visit"]
    df_clinical_complete = df_imaging.merge(
        df_demographics, how="inner", on=merge_key
    ).drop(columns="diagnosis")
    df_clinical = df_clinical.dropna(subset=merge_key)
    return df_clinical_complete.merge(
        df_clinical[
            [
                "blinded_code",
                "blinded_site",
                "visit",
                "diagnosis",
                "ftld-cdr-global",
                "cdr-sob",
            ]
        ],
        how="inner",
        on=merge_key,
    )


def dataset_to_bids(complete_data_df: DataFrame, gif: bool) -> Dict[str, DataFrame]:
    """Selects the data needed to write the participants, sessions, and scans tsvs.

    Parameters
    ----------
    complete_data_df: DataFrame
        Dataframe containing the merged data extracted from the raw images and the clinical data

    gif: bool
        If True, indicates the user wants to have the values of the gif parcellation

    Returns
    -------
    Dict[str, DataFrame]
        Dictionary containing as key participants, sessions and scans, and the values wanted for each tsv
    """
    import os

    # generates participants, sessions and scans tsv
    complete_data_df = complete_data_df.drop_duplicates(
        subset=["participant_id", "session_id", "modality", "run_num", "bids_filename"]
    ).set_index(
        ["participant_id", "session_id", "modality", "run_num", "bids_filename"],
        verify_integrity=True,
    )
    # open the reference for building the tsvs:
    path_to_ref_csv = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "data",
        "genfi_ref.csv",
    )
    df_ref = pd.read_csv(path_to_ref_csv, sep=";")

    if not gif:
        df_ref = df_ref.head(8)
    return {
        col: complete_data_df.filter(items=list(df_ref[col]))
        for col in ["participants", "sessions", "scans"]
    }


def intersect_data(
    imaging_data: DataFrame, df_clinical_complete: DataFrame
) -> DataFrame:
    """This function merges the dataframe containing the data extracted
    from the raw images and from the clinical data.

    Parameters
    ----------
    imaging_data: DataFrame
        Dataframe containing the data extracted from the raw images

    df_clinical_complete: DataFrame
        Dataframe containing the clinical data

    Returns
    -------
    df_complete: DataFrame
        Dataframe containing the merged data
    """
    df_complete = imaging_data.merge(
        df_clinical_complete,
        how="inner",
        left_on=["source_id", "source_ses_id"],
        right_on=["blinded_code", "visit"],
    )
    return df_complete.loc[:, ~df_complete.columns.duplicated()]


def read_imaging_data(source_path: PathLike) -> DataFrame:
    """This function finds the imaging data and filters it.

    Parameters
    ----------
    source_path: PathLike
        Path to the raw data

    Returns
    -------
    df_dicom: DataFrame
        Dataframe containing the data extracted.
    """

    return filter_dicoms(
        pd.DataFrame(find_dicoms(source_path), columns=["source_path", "source"])
    )


def compute_baseline_date(df: DataFrame) -> DataFrame:
    """Computes the baseline date by taking the minimum
    acq_date for each subject.

    Parameters
    ----------
    df: Dataframe

    Returns
    -------
    Dataframe
        Contains the baseline date.
    """
    df_1 = (
        df[["source_id", "source_ses_id", "acq_date"]]
        .groupby(["source_id", "source_ses_id"])
        .min()
    )
    df_2 = df[["source_id", "acq_date"]].groupby(["source_id"]).min()
    return df_1.join(df_2.rename(columns={"acq_date": "baseline"}))


def compute_session_numbers(df: DataFrame) -> DataFrame:
    """Computes the session IDs obtained from the number of months between an acquisition and the baseline acquisition.

    Parameters
    ----------
    df: DataFrame
        Dataframe containing the timestamp of the acquisition.

    Returns
    -------
    DataFrame
        Dataframe containing the session_id computed from the timestamps.
    """
    from datetime import datetime

    for col in ("acq_date", "baseline"):
        df[col] = df[col].apply(lambda x: datetime.strptime(x, "%Y%m%d"))
    return df.assign(
        ses_month=lambda x: 12
        * (x.acq_date.apply(lambda y: y.year) - x.baseline.apply(lambda y: y.year))
        + (x.acq_date.apply(lambda y: y.month) - x.baseline.apply(lambda y: y.month)),
        session_id=lambda x: x.ses_month.map(lambda y: f"ses-M{y:03d}"),
    )


def compute_modality(df: DataFrame) -> DataFrame:
    """Parses the series_desc column to identify the modality, and map
    this modality to a metadata dictionary.

    Parameters
    ----------
    df: Dataframe
        DataFrame on which to compute the modality.

    Returns
    -------
    Dataframe
        Contains the datatype, suffix, sidecars and task that may be needed to form the complete file name.
    """
    from clinica.iotools.bids_utils import identify_modality

    modality_mapping = {
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
            "task": "task-rest",
        },
        "fieldmap": {
            "datatype": "fmap",
            "suffix": "fmap",
            "sidecars": ["fmap.json"],
            "task": "",
        },
    }
    df = (
        df.assign(
            modality=lambda x: x.series_desc.apply(lambda y: identify_modality(y))
        )
        .dropna(subset=["modality"])
        .drop_duplicates()
    )
    return df.join(df.modality.map(modality_mapping).apply(pd.Series))


def identify_fieldmaps(df: DataFrame) -> DataFrame:
    """Identifies the fieldmaps imaging type: magnitude or phase difference.

    Parameters
    ----------
    df:Dataframe
        Dataframe without the fieldmaps identified.

    Returns
    -------
    Dataframe
        Dataframe with the fieldmaps identified
    """
    filter = ["source_id", "source_ses_id", "modality", "dir_num", "suffix"]
    df1 = df[filter][df["modality"].str.contains("fieldmap")].groupby(filter[:-1]).min()
    df2 = (
        df[filter[:-1]][df["modality"].str.contains("fieldmap")]
        .groupby(filter[:-2])
        .min()
    )
    df1 = df1.join(df2.rename(columns={"dir_num": "run_01_dir_num"}))
    df_alt = df1.reset_index().assign(
        fmap_type_p=lambda x: (x.run_01_dir_num != x.dir_num)
    )
    return df_alt.assign(
        fmap=lambda x: x.fmap_type_p.apply(lambda y: "phasediff" if y else "magnitude")
    )


def merge_fieldmaps(df: DataFrame, df_fmap: DataFrame) -> DataFrame:
    """Merges the dataframe containing the fieldmaps names

    Parameters
    ----------
    df: Dataframe
        Initial dataframe. Fieldmaps are identified on information extracted here.

    df_fmap: Dataframe
        Dataframe in which fieldmaps are identified.

    Returns
    -------
    Dataframe
        Dataframe containing the correct fieldmap for each fieldmap acquisition.
    """
    filter = ["source_id", "source_ses_id", "modality", "dir_num", "fmap"]
    df1 = df.merge(df_fmap[filter], how="inner", on=filter[:-1])
    df2 = df.merge(df_fmap[filter], how="left", on=filter[:-1])
    df2 = df2[~df2["modality"].str.contains("fieldmap")]
    df1 = df1.assign(suffix=lambda x: x.fmap)
    return pd.concat([df2, df1], ignore_index=True)


def compute_runs(df: DataFrame) -> DataFrame:
    """This functions computes the run numbers.

    Parameters
    ----------
    df: DataFrame
        Dataframe without runs.

    Returns
    -------
    DataFrame
        Dataframe containing the correct run for each acquisition.
    """

    filter = ["source_id", "source_ses_id", "suffix", "number_of_parts", "dir_num"]
    df1 = df[filter].groupby(filter).min()
    df2 = df[filter].groupby(filter[:-1]).min()
    df1 = df1.join(df2.rename(columns={"dir_num": "run_01_dir_num"}))
    df_alt = df1.reset_index().assign(run=lambda x: (x.run_01_dir_num != x.dir_num))

    df_run = pd.concat(
        [
            df_alt,
            pd.DataFrame(
                _compute_scan_sequence_numbers(df_alt.run.tolist()),
                columns=["run_number"],
            ),
        ],
        axis=1,
    )
    df_run = df_run.assign(
        run_num=lambda df: df.run_number.apply(lambda x: f"run-{x:02d}")
    )
    return df_run


def compute_philips_parts(df: DataFrame) -> DataFrame:
    """Compute the parts numbers for philips dwi acquisitions.
    The amount of dwi acquisitions linked together is indicated.
    For example, if a dwi acquisition is split in 9,
    the `number_of_parts` column will have a value of 9 for all of these acquisitions.
    Two columns are added:
        - part_number, which contains the number for each part of a DTI acquisition.
        - number_of_parts, which contains the amount of parts for each DTI acquisition.

    Parameters
    ----------
    df: DataFrame
        Dataframe without runs.

    Returns
    -------
    DataFrame
        Dataframe containing the correct dwi part number for each acquisition. It also contains
        the total amount of dwi parts for each subjects-session.
    """
    df_with_duplicate_flags = _find_duplicate_run(df)
    df_with_part_numbers = _compute_part_numbers(df_with_duplicate_flags)
    df_with_number_of_parts = _compute_number_of_parts(df_with_part_numbers)
    return pd.concat(
        [df_with_part_numbers, df_with_number_of_parts["number_of_parts"]], axis=1
    )


def _find_duplicate_run(df: DataFrame) -> DataFrame:
    """Create a column that contains the information of whether a run is a duplicate or not."""
    filter = ["source_id", "source_ses_id", "suffix", "dir_num"]
    df = df[df["suffix"].str.contains("dwi", case=False)]
    df1 = df[filter].groupby(filter).min()
    df2 = df[filter].groupby(filter[:-1]).min()
    df1 = df1.join(df2.rename(columns={"dir_num": "part_01_dir_num"}))
    df_alt = df1.reset_index().assign(run=lambda x: (x.part_01_dir_num != x.dir_num))
    return df_alt


def _compute_part_numbers(df: DataFrame) -> DataFrame:
    """Compute the sequence number of each part."""
    return pd.concat(
        [
            df,
            pd.DataFrame(
                _compute_scan_sequence_numbers(df.run.tolist()), columns=["part_number"]
            ),
        ],
        axis=1,
    )


def _compute_number_of_parts(df: pd.DataFrame) -> DataFrame:
    """Add the number of parts (the max value of the part_number column) to each part."""
    filter2 = ["source_id", "source_ses_id", "suffix", "part_number"]
    df_parts_1 = df[filter2].groupby(filter2).max()
    df_parts_2 = df[filter2].groupby(filter2[:-1]).max()
    df_max_nb_parts = df_parts_1.join(
        df_parts_2.rename(columns={"part_number": "number_of_parts"})
    ).reset_index()
    return df_max_nb_parts


def _compute_scan_sequence_numbers(duplicate_flags: Iterable[bool]) -> List[int]:
    """Return the run number from an iterable of booleans indicating
    whether each scan is a duplicate or not.

    Parameters
    ---------
    duplicate_flags : Iterable[bool]
        If the element at index k is True, then the scan k is a duplicate.
        Otherwise, it is the first scan of the sequence.

    Returns
    ------
    ses_numbers : List[int]
        The list of scan sequence numbers.

    Examples
    ---------
    >>> _compute_scan_sequence_numbers([False, True, False, True, False, False, False, True, False, True, True])
    [1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 3]

    Raises
    -----
    ValueError :
        If the input list is empty.
    """
    if len(duplicate_flags) == 0:
        raise ValueError("Provided list is empty.")
    ses_numbers = [1]
    for k in range(1, len(duplicate_flags)):
        ses_numbers.append(1 if not duplicate_flags[k] else ses_numbers[k - 1] + 1)
    return ses_numbers


def merge_imaging_data(df_dicom: DataFrame) -> DataFrame:
    """This function uses the raw information extracted from the images,
    to obtain all the information necessary for the BIDS conversion.

    Parameters
    ----------
    df_dicom: DataFrame
        Dataframe containing the data extracted from the images

    Returns
    -------
    df_sub_ses_run: DataFrame
        Dataframe with the data necessary for the BIDS
    """
    from clinica.utils.filemanip import get_parent

    df_dicom = df_dicom.assign(
        source_id=lambda df: df.source.apply(
            lambda x: get_parent(x, 2).name.split("-")[0]
        ),
        source_ses_id=lambda df: df.source.apply(
            lambda x: get_parent(x, 2).name.split("-")[1]
        ),
        origin=lambda df: df.source.apply(lambda x: get_parent(x, 4)),
    )

    df_dicom = df_dicom.reset_index()

    df_dicom = df_dicom.assign(
        source_ses_id=lambda df: df.source_ses_id.astype("int"),
        genfi_version=lambda df: df.source_ses_id.apply(
            lambda x: f"GENFI{len(str(x))}"
        ),
    )
    df_1 = compute_baseline_date(df_dicom)

    df_1 = compute_session_numbers(df_1)

    df_sub_ses = df_dicom.merge(
        df_1[["baseline", "session_id"]], how="inner", on=["source_id", "source_ses_id"]
    )
    df_sub_ses = df_sub_ses.assign(
        participant_id=lambda df: df.source_id.apply(lambda x: f"sub-{x}")
    )
    df_sub_ses = compute_modality(df_sub_ses)

    df_fmap = df_sub_ses.assign(
        dir_num=lambda x: x.source.apply(
            lambda y: int(get_parent(y).name.split("-")[0])
        )
    )

    df_suf = merge_fieldmaps(df_fmap, identify_fieldmaps(df_fmap))

    df_suf_dir = df_suf.assign(
        dir_num=lambda x: x.source.apply(
            lambda y: int(get_parent(y).name.split("-")[0])
        )
    )
    df_alt = compute_philips_parts(df_suf_dir)
    df_parts = df_suf_dir.merge(
        df_alt[["source_id", "source_ses_id", "suffix", "dir_num", "number_of_parts"]],
        how="left",
        on=["source_id", "source_ses_id", "suffix", "dir_num"],
    )
    df_parts[["number_of_parts"]] = df_parts[["number_of_parts"]].fillna(value="1")

    df_alt = compute_runs(df_parts)

    df_sub_ses_run = df_suf_dir.merge(
        df_alt[
            [
                "source_id",
                "source_ses_id",
                "suffix",
                "dir_num",
                "run_num",
                "number_of_parts",
            ]
        ],
        how="left",
        on=["source_id", "source_ses_id", "suffix", "dir_num"],
    )
    df_sub_ses_run = df_sub_ses_run.assign(
        bids_filename=lambda df: df[
            ["participant_id", "session_id", "run_num", "task", "suffix"]
        ].agg("_".join, axis=1),
    )
    df_sub_ses_run = df_sub_ses_run.assign(
        bids_filename=lambda df: df.bids_filename.apply(lambda x: x.replace("__", "_"))
    )
    df_sub_ses_run = df_sub_ses_run.assign(
        bids_full_path=lambda df: df[
            ["participant_id", "session_id", "datatype", "bids_filename"]
        ].agg("/".join, axis=1)
    )
    return df_sub_ses_run


def write_bids(
    to: PathLike,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
) -> None:
    """This function writes the BIDS

    Parameters
    ----------
    to: PathLike
        Path where the BIDS should be written

    participants: DataFrame
        DataFrame containing the data for the participants.tsv

    sessions: DataFrame
        DataFrame containing the data for the sessions.tsv

    scans: DataFrame
        DataFrame containing the data for the scans.tsv
    """
    import os
    from pathlib import Path

    from fsspec.implementations.local import LocalFileSystem

    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription
    from clinica.iotools.bids_utils import run_dcm2niix, write_to_tsv
    from clinica.utils.stream import cprint

    cprint("Starting to write the BIDS.", lvl="info")
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
        dcm2niix_success = run_dcm2niix(
            Path(metadata["source_path"]).parent,
            to / str(Path(bids_full_path).parent),
            metadata["bids_filename"],
            True,
        )
        if (
            "dwi" in metadata["bids_filename"]
            and "Philips" in metadata.manufacturer
            and dcm2niix_success
        ):
            merge_philips_diffusion(
                to / Path(bids_full_path).with_suffix(".json"),
                metadata.number_of_parts,
                metadata.run_num,
            )
    correct_fieldmaps_name(to)
    return


def correct_fieldmaps_name(to: PathLike) -> None:
    """This function scans the BIDS after it has been written to correct the nameds of the fieldmap files.

    Parameters
    ----------
    to: PathLike
        Path to the BIDS

    Examples
    --------
    bids
    ├── README
    ├── dataset_description.json
    ├── participants.tsv
    └── sub-GRN001
        ├── ses-M000
        │   ├── fmap
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude_e1.json
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude_e1.nii.gz
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude_e2.json
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude_e2.nii.gz
        │   │   ├── sub-GRN001_ses-M000_run-01_phasediff_e2_ph.json
        │   │   └── sub-GRN001_ses-M000_run-01_phasediff_e2_ph.nii.gz
        │   └── func
        │       ├── sub-GRN001_ses-M000_run-01_bold.json
        │       └── sub-GRN001_ses-M000_run-01_bold.nii.gz
        └── sub-GRN001_sessions.tsv
    >>> correct_fieldmaps_name("path/to/bids)
    bids
    ├── README
    ├── dataset_description.json
    ├── participants.tsv
    └── sub-GRN001
        ├── ses-M000
        │   ├── fmap
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude1.json
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude1.nii.gz
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude2.json
        │   │   ├── sub-GRN001_ses-M000_run-01_magnitude2.nii.gz
        │   │   ├── sub-GRN001_ses-M000_run-01_phasediff.json
        │   │   └── sub-GRN001_ses-M000_run-01_phasediff.nii.gz
        │   └── func
        │       ├── sub-GRN001_ses-M000_run-01_bold.json
        │       └── sub-GRN001_ses-M000_run-01_bold.nii.gz
        └── sub-GRN001_sessions.tsv
    """
    import os
    import re

    for z in Path(to).rglob("*magnitude_e*"):
        os.rename(z, z.parent / re.sub(r"magnitude_e", "magnitude", z.name))

    for z in Path(to).rglob("*phasediff_e*_ph*"):
        os.rename(z, z.parent / re.sub(r"phasediff_e[1-9]_ph", "phasediff", z.name))


def merge_philips_diffusion(
    json_file: PathLike,
    number_of_parts: float,
    run_num: str,
) -> None:
    """Add the dwi number in the provided JSON file for each run of Philips images.
    The 'MultipartID' field of the input JSON file is set to 'dwi_1' or 'dwi_2' depending
    on the run number.
    """
    import json

    multipart_id = _get_multipart_id(
        PhilipsNumberOfParts.from_int(int(number_of_parts)), run_num
    )
    with open(json_file, "r+") as f:
        if multipart_id is not None:
            data = json.load(f)
            data["MultipartID"] = multipart_id
            json.dump(data, f, indent=4)


class PhilipsNumberOfParts(Enum):
    """DWI scans obtained with a Philips scanner might have
    been divided in either 5 or 9 parts. This distinction is important
    because we will link these different parts together.
    If the number of parts is not 5 or 9, nothing will be done.
    """

    FIVE = 5
    NINE = 9
    OTHER = None

    @classmethod
    def from_int(cls, nb_parts: int):
        import warnings

        for member in cls:
            if member.value == nb_parts:
                return member
        warnings.warn(
            f"Unexpected number of splits {nb_parts}. "
            f"Should be one of {[c.value for c in cls]}."
        )
        return cls.OTHER


def _get_multipart_id(nb_parts: PhilipsNumberOfParts, run_num: str) -> Optional[str]:
    """Return the MultiPartID for the provided run number depending on the number of parts."""
    if nb_parts == PhilipsNumberOfParts.NINE:
        if run_num in (f"run-0{k}" for k in range(1, 5)):
            return "dwi_2"
        if run_num in (f"run-0{k}" for k in range(5, 10)):
            return "dwi_1"
        raise ValueError(f"{run_num} is outside of the scope.")
    if nb_parts == PhilipsNumberOfParts.FIVE:
        if run_num in (f"run-0{k}" for k in range(1, 6)):
            return "dwi_1"
        raise ValueError(f"{run_num} is outside of the scope.")
    if nb_parts == PhilipsNumberOfParts.OTHER:
        return None
