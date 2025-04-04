from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import pydicom as pdcm
from pandas import DataFrame

__all__ = [
    "merge_imaging_and_clinical_data",
    "parse_clinical_data",
    "parse_imaging_data",
    "prepare_dataset_to_bids_format",
    "write_bids",
]


def _find_dicoms(path_to_source_data: Path) -> Iterable[Tuple[Path, Path]]:
    """Find the dicoms in the given directory.

    Parameters
    ----------
    path_to_source_data: Path
        The path to the source data.

    Returns
    -------
    Iterable[Tuple[Path, Path]]
        Path to found files and parent directory
    """
    from clinica.utils.stream import cprint

    cprint("Looking for imaging data.", lvl="info")
    for z in path_to_source_data.rglob("*.dcm"):
        yield z, z.parent


def _filter_dicoms(df: DataFrame) -> DataFrame:
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
    df = df.sort_values(by="source_path").reset_index(drop=True)
    df = df.drop_duplicates(subset=["source"])
    df = df.assign(
        series_desc=lambda x: x.source_path.apply(
            lambda y: _handle_series_description(y)
        ),
        acq_date=lambda x: x.source_path.apply(lambda y: pdcm.dcmread(y).StudyDate),
        manufacturer=lambda x: x.source_path.apply(lambda y: _handle_manufacturer(y)),
    )
    df = df.set_index(["source_path"], verify_integrity=True)
    df = df[df["source"].str.contains("secondary", case=False) != True]  # noqa
    for file_mod in to_filter:
        df = df[~df["series_desc"].str.contains(file_mod, case=False)]
    return df


def _handle_series_description(x: str) -> str:
    """This function handles the series description used to identify the modality later on.

    In the case of T1 and T2, some dicoms are misnamed. The echo time and repetition time
    are supposed to be very different between a T1 and T2. Therefore, we use this proxy to check that the series description is coherent.
    The values used are in milliseconds, and have been chosen arbitrarily.
    """
    import warnings

    MAX_ECHO_TIME_FOR_A_T1 = 60
    MAX_REPETITION_TIME_FOR_A_T1 = 2100
    series_description = pdcm.dcmread(x).SeriesDescription
    if any([modality in series_description.lower() for modality in ["t1", "t2"]]):
        try:
            echo_time = pdcm.dcmread(x).EchoTime
            repetition_time = pdcm.dcmread(x).RepetitionTime
            if (
                "t2" in series_description.lower()
                and echo_time < MAX_ECHO_TIME_FOR_A_T1
                and repetition_time < MAX_REPETITION_TIME_FOR_A_T1
            ):
                return "t1"
        except AttributeError:
            warnings.warn(
                message=(
                    f"The subject from DICOM {Path(x).parent} has no Echo Time "
                    "or Repetition Time, it might be wrongly converted."
                )
            )
    return series_description


def _handle_manufacturer(x: str) -> str:
    try:
        return pdcm.dcmread(x).Manufacturer
    except Exception:
        return "Unknown"


def _check_file(directory: Path, pattern: str) -> Path:
    try:
        data_file = [f for f in directory.glob(pattern)]
    except StopIteration:
        raise FileNotFoundError("Clinical data file not found.")
    if len(data_file) == 0:
        raise FileNotFoundError("Clinical data not found or incomplete. Aborting")
    if len(data_file) > 1:
        raise ValueError("Too many data files found, expected one. Aborting.")
    return data_file[0]


def parse_clinical_data(clinical_data_directory: Path) -> DataFrame:
    return _complete_clinical_data(*_find_clinical_data(clinical_data_directory))


def _find_clinical_data(
    clinical_data_directory: Path,
) -> tuple[DataFrame, DataFrame, DataFrame, DataFrame, DataFrame]:
    """Finds the clinical data associated with the dataset.

    Parameters
    ----------
    clinical_data_directory: Path
        The path to the clinical data.

    Returns
    -------
    List[DataFrame]
        Dataframes containing the clinical data
    """
    from clinica.utils.stream import cprint

    cprint("Looking for clinical data.", lvl="info")

    return tuple(
        _read_file(_check_file(clinical_data_directory, pattern))
        for pattern in (
            "FINAL*DEMOGRAPHICS*.xlsx",
            "FINAL*IMAGING*.xlsx",
            "FINAL*CLINICAL*.xlsx",
            "FINAL*BIOSAMPLES*.xlsx",
            "FINAL*NEUROPSYCH*.xlsx",
        )
    )


def _read_file(data_file: Path) -> pd.DataFrame:
    return (
        pd.concat(
            [
                pd.read_excel(data_file),
                pd.read_excel(data_file, sheet_name=1),
            ]
        )
        .convert_dtypes()
        .rename(columns=lambda x: x.lower().replace(" ", "_"))
    )


def _complete_clinical_data(
    df_demographics: DataFrame,
    df_imaging: DataFrame,
    df_clinical: DataFrame,
    df_biosamples: DataFrame,
    df_neuropsych: DataFrame,
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

    df_biosamples: DataFrame
        Dataframe containing the biosample data

    df_neuropsych: DataFrame
        Dataframe containing the neuropsych data

    Returns
    -------
    df_clinical_complete: DataFrame
        Dataframe with the data of the 3 input dataframes
    """
    merge_key = ["blinded_code", "blinded_site", "visit"]
    df_clinical_complete = df_imaging.merge(
        df_demographics, how="inner", on=merge_key
    ).drop(columns="diagnosis")
    df_clinical_complete = df_clinical_complete.merge(
        df_biosamples, how="inner", on=merge_key
    )
    df_clinical_complete = df_clinical_complete.merge(
        df_neuropsych, how="inner", on=merge_key
    )
    df_clinical = df_clinical.dropna(subset=merge_key)
    return df_clinical_complete.merge(df_clinical, how="inner", on=merge_key)


def prepare_dataset_to_bids_format(
    complete_data_df: DataFrame,
    gif: bool,
    path_to_clinical_tsv: Path,
) -> Dict[str, DataFrame]:
    """Selects the data needed to write the participants, sessions, and scans tsvs.

    Parameters
    ----------
    complete_data_df: DataFrame
        Dataframe containing the merged data extracted from the raw images and the clinical data

    gif: bool
        If True, indicates the user wants to have the values of the gif parcellation

    path_to_clinical_tsv: Path
        TSV file containing the data fields the user wishes to have from the excel spreadsheets

    Returns
    -------
    Dict[str, DataFrame]
        Dictionary containing as key participants, sessions and scans, and the values wanted for each tsv
    """
    complete_data_df = complete_data_df.drop_duplicates(
        subset=["participant_id", "session_id", "modality", "run_num", "bids_filename"]
    ).set_index(
        ["participant_id", "session_id", "modality", "run_num", "bids_filename"],
        verify_integrity=True,
    )
    specifications = pd.read_csv(Path(__file__).parent / "specifications.csv", sep=";")
    if not gif:
        specifications = specifications.head(8)
    # add additional data through csv
    if path_to_clinical_tsv:
        additional_data_df = pd.read_csv(path_to_clinical_tsv, sep="\t")
        data_mapping = pd.read_csv(Path(__file__) / "data_mapping.tsv", sep="\t")
        pre_addi_df = data_mapping.merge(additional_data_df, how="inner", on="data")
        addi_df = pd.DataFrame(
            [
                pre_addi_df["data"][pre_addi_df["dest"] == x].values.tolist()
                for x in ("participants", "sessions", "scans")
            ]
        ).transpose()
        addi_df.columns = ["participants", "sessions", "scans"]
        df_to_write = pd.concat([specifications, addi_df])
    else:
        df_to_write = specifications
    return {
        col: complete_data_df.filter(items=list(df_to_write[col]))
        for col in ["participants", "sessions", "scans"]
    }


def merge_imaging_and_clinical_data(
    imaging_data: DataFrame, clinical_data: DataFrame
) -> DataFrame:
    """This function merges the dataframe containing the data extracted
    from the raw images and from the clinical data.

    Parameters
    ----------
    imaging_data : DataFrame
        Dataframe containing the data extracted from the raw images

    clinical_data : DataFrame
        Dataframe containing the clinical data

    Returns
    -------
    df_complete : DataFrame
        Dataframe containing the merged data
    """
    df_complete = imaging_data.merge(
        clinical_data,
        how="inner",
        left_on=["source_id", "source_ses_id"],
        right_on=["blinded_code", "visit"],
    )
    return df_complete.loc[:, ~df_complete.columns.duplicated()]


def parse_imaging_data(source_path: Path) -> DataFrame:
    return _merge_imaging_data(_read_imaging_data(source_path))


def _read_imaging_data(source_path: Path) -> DataFrame:
    """This function finds the imaging data and filters it.

    Parameters
    ----------
    source_path: Path
        The path to the raw data.

    Returns
    -------
    df_dicom: DataFrame
        Dataframe containing the data extracted.
    """
    return _filter_dicoms(
        pd.DataFrame(_find_dicoms(source_path), columns=["source_path", "source"])
    )


def _merge_imaging_data(df: DataFrame) -> DataFrame:
    """This function uses the raw information extracted from the images,
    to obtain all the information necessary for the BIDS conversion.

    Parameters
    ----------
    df: DataFrame
        Dataframe containing the data extracted from the images

    Returns
    -------
    DataFrame :
        Dataframe with the data necessary for the BIDS
    """
    df = _compute_source_id_and_source_ses_id(df).reset_index()
    df = _compute_genfi_version(df)
    df = _compute_baseline_and_session_numbers(df)
    df = _compute_participant_id(df)
    df = _compute_modality(df)
    df = _compute_fieldmaps(df)
    df = _compute_runs(df)
    df = _compute_bids_full_path(df)

    return df


def _compute_source_id_and_source_ses_id(df: DataFrame) -> DataFrame:
    """Adds two columns built from the column 'source'.

    - 'source': subject ID and session ID joined by a dash
    - 'source_id': subject ID in the raw dataset
    - 'source_ses_id': Session ID in the raw dataset

    Example
    -------
    source = "'C9ORF001-11'
    source_id = 'C9ORF001'
    source_ses_id = '11
    """
    from clinica.utils.filemanip import get_parent

    return df.assign(
        source_id=lambda x: x.source.apply(
            lambda y: get_parent(y, 2).name.split("-")[0]
        ),
        source_ses_id=lambda x: x.source.apply(
            lambda y: get_parent(y, 2).name.split("-")[1]
        ),
    )


def _compute_genfi_version(df: DataFrame) -> DataFrame:
    """Compute the genfi_version from the souce_ses_id column.

    The column `source_ses_id` is casted to integer values.
    """
    return df.assign(
        source_ses_id=lambda x: x.source_ses_id.astype("int"),
        genfi_version=lambda x: x.source_ses_id.apply(lambda y: f"GENFI{len(str(y))}"),
    )


def _compute_baseline_and_session_numbers(df: DataFrame) -> DataFrame:
    """Compute the baseline date and session numbers
    and merge on 'source_id` and `source_ses_id`.
    """
    df_with_acq_date = _compute_baseline_date(df)
    df_with_acq_date_and_session_numbers = _compute_session_numbers(df_with_acq_date)
    return df.merge(
        df_with_acq_date_and_session_numbers[["baseline", "session_id"]],
        how="inner",
        on=["source_id", "source_ses_id"],
    )


def _compute_baseline_date(df: DataFrame) -> DataFrame:
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
    df_2 = df[["source_id", "acq_date"]].groupby("source_id").min()
    return df_1.join(df_2.rename(columns={"acq_date": "baseline"}))


def _compute_session_numbers(df: DataFrame) -> DataFrame:
    """Computes the session IDs obtained from the number of months between
    an acquisition and the baseline acquisition.

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


def _compute_participant_id(df: DataFrame) -> DataFrame:
    """Compute the 'participant_id' column from the 'source_id' column."""
    from clinica.converters.study_models import StudyName, bids_id_factory

    return df.assign(
        participant_id=df.source_id.apply(
            lambda x: bids_id_factory(StudyName.GENFI).from_original_study_id(x)
        )
    )


def _compute_modality(df: DataFrame) -> DataFrame:
    """Parses the series_desc column to identify the modality, and map
    this modality to a metadata dictionary.

    Parameters
    ----------
    df: Dataframe
        DataFrame on which to compute the modality.

    Returns
    -------
    Dataframe
        Contains the datatype, suffix, sidecars and task that may be
        needed to form the complete file name.
    """
    from clinica.converters._utils import identify_modality

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


def _compute_fieldmaps(df: DataFrame) -> DataFrame:
    """Compute the fieldmaps. Please add details...
    For example, why do we need to compute the dir_num before and after ?
    """
    from clinica.utils.filemanip import get_parent

    df_with_dir_num = df.assign(
        dir_num=lambda x: x.source.apply(
            lambda y: int(get_parent(y).name.split("-")[0])
        )
    )
    df_suf = _merge_fieldmaps(df_with_dir_num, _identify_fieldmaps(df_with_dir_num))
    return df_suf


def _identify_fieldmaps(df: DataFrame) -> DataFrame:
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
    column_filter = ["source_id", "source_ses_id", "modality", "dir_num", "suffix"]
    df1 = (
        df[column_filter][df["modality"].str.contains("fieldmap")]
        .groupby(column_filter[:-1])
        .min()
    )
    df2 = (
        df[column_filter[:-1]][df["modality"].str.contains("fieldmap")]
        .groupby(column_filter[:-2])
        .min()
    )
    df1 = df1.join(df2.rename(columns={"dir_num": "run_01_dir_num"}))
    df_alt = df1.reset_index().assign(
        fmap_type_p=lambda x: (x.run_01_dir_num != x.dir_num)
    )
    return df_alt.assign(
        fmap=lambda x: x.fmap_type_p.apply(lambda y: "phasediff" if y else "magnitude")
    )


def _merge_fieldmaps(df: DataFrame, df_fmap: DataFrame) -> DataFrame:
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
    column_filter = ["source_id", "source_ses_id", "modality", "dir_num", "fmap"]
    df1 = df.merge(df_fmap[column_filter], how="inner", on=column_filter[:-1])
    df2 = df.merge(df_fmap[column_filter], how="left", on=column_filter[:-1])
    df2 = df2[~df2["modality"].str.contains("fieldmap")]
    df1 = df1.assign(suffix=lambda x: x.fmap)
    return pd.concat([df2, df1], ignore_index=True)


def _compute_runs(df: DataFrame) -> DataFrame:
    """Compute the run numbers."""
    df_with_number_of_parts = _compute_number_of_scan_parts(df)
    df_with_run_numbers = _compute_run_numbers_from_parts(df_with_number_of_parts)
    return df.merge(
        df_with_run_numbers[
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


def _compute_number_of_scan_parts(df: DataFrame) -> DataFrame:
    """Compute the number of parts. Explain me please..."""
    df_alt = _compute_philips_parts(df)
    df_parts = df.merge(
        df_alt[["source_id", "source_ses_id", "suffix", "dir_num", "number_of_parts"]],
        how="left",
        on=["source_id", "source_ses_id", "suffix", "dir_num"],
    )
    df_parts[["number_of_parts"]] = df_parts[["number_of_parts"]].fillna(value="1")
    return df_parts


def _compute_philips_parts(df: DataFrame) -> DataFrame:
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
    column_filter = ["source_id", "source_ses_id", "suffix", "dir_num"]
    df = df[df["suffix"].str.contains("dwi", case=False)]
    df1 = df[column_filter].groupby(column_filter).min()
    df2 = df[column_filter].groupby(column_filter[:-1]).min()
    df1 = df1.join(df2.rename(columns={"dir_num": "part_01_dir_num"}))
    return df1.reset_index().assign(run=lambda x: (x.part_01_dir_num != x.dir_num))


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
    column_filter = ["source_id", "source_ses_id", "suffix", "part_number"]
    df_parts_1 = df[column_filter].groupby(column_filter).max()
    df_parts_2 = df[column_filter].groupby(column_filter[:-1]).max()
    return df_parts_1.join(
        df_parts_2.rename(columns={"part_number": "number_of_parts"})
    ).reset_index()


def _compute_run_numbers_from_parts(df: DataFrame) -> DataFrame:
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
    column_filter = [
        "source_id",
        "source_ses_id",
        "suffix",
        "number_of_parts",
        "dir_num",
    ]
    df1 = df[column_filter].groupby(column_filter).min()
    df2 = df[column_filter].groupby(column_filter[:-1]).min()
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
    return df_run.assign(run_num=lambda x: x.run_number.apply(lambda y: f"run-{y:02d}"))


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


def _compute_bids_full_path(df: DataFrame) -> DataFrame:
    """Compute the BIDS full path."""
    df_with_bids_filename = df.assign(
        bids_filename=lambda x: x[
            ["participant_id", "session_id", "task", "run_num", "suffix"]
        ].agg("_".join, axis=1),
    )
    df_with_bids_filename = df_with_bids_filename.assign(
        bids_filename=lambda x: x.bids_filename.apply(lambda y: y.replace("__", "_"))
    )
    return df_with_bids_filename.assign(
        bids_full_path=lambda x: x[
            ["participant_id", "session_id", "datatype", "bids_filename"]
        ].agg("/".join, axis=1)
    )


def _drop_duplicate_line_with_nans(participants: pd.DataFrame) -> pd.DataFrame:
    """Performs an operation specific to participants metadata file in GENFI :
    There can subsist two lines per participant in case one holds the information and one holds some nans.
    It is not recommended to just use dropna() since it is possible there is no information for one subject, which would be lost
    The proposed solution here drops the line with nans if there is another line for the participant.

    Parameters
    ----------
    participants :
        The participants metadata frame

    Returns
    -------
        The same frame with dropped duplicate lines with nan values

    """
    from clinica.utils.exceptions import ClinicaBIDSError

    if "participant_id" not in participants.columns:
        raise ClinicaBIDSError(
            "Column participant_id was not found in the participants tsv while it is required by BIDS specifications."
        )

    for participant, size in participants.groupby(["participant_id"]).size().items():
        if size > 1:
            participants.drop(
                participants[
                    (participants["participant_id"] == participant)
                    * (participants.isna().any(axis=1))
                ].index,
                inplace=True,
            )
    return participants


def write_bids(
    to: Path,
    participants: DataFrame,
    sessions: DataFrame,
    scans: DataFrame,
    source: Path,
) -> None:
    """This function writes the BIDS

    Parameters
    ----------
    to: Path
        The path where the BIDS should be written

    participants: DataFrame
        DataFrame containing the data for the participants.tsv

    sessions: DataFrame
        DataFrame containing the data for the sessions.tsv

    scans: DataFrame
        DataFrame containing the data for the scans.tsv

    source : Path
        Path to the source imaging data
    """
    import os

    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters._utils import run_dcm2niix, write_to_tsv
    from clinica.dataset import BIDSDatasetDescription
    from clinica.utils.stream import cprint

    cprint("Starting to write the BIDS.", lvl="info")
    fs = LocalFileSystem(auto_mkdir=True)
    # Ensure BIDS hierarchy is written first.

    participants = (
        participants.reset_index()
        .drop(["session_id", "modality", "run_num", "bids_filename", "source"], axis=1)
        .drop_duplicates()
    )
    participants = _drop_duplicate_line_with_nans(participants).set_index(
        "participant_id"
    )

    with fs.transaction:
        with fs.open(to / "dataset_description.json", "w") as dataset_description_file:
            BIDSDatasetDescription(name="GENFI").write(to=dataset_description_file)
        with fs.open(to / "participants.tsv", "w") as participant_file:
            write_to_tsv(participants, participant_file)

    for participant_id, data_frame in sessions.groupby("participant_id"):
        sessions = data_frame.droplevel(
            ["participant_id", "modality", "bids_filename", "run_num"]
        ).drop_duplicates()

        sessions_filepath = to / str(participant_id) / f"{participant_id}_sessions.tsv"
        with fs.open(sessions_filepath, "w") as sessions_file:
            write_to_tsv(sessions, sessions_file)

    scans = scans.reset_index().set_index(["bids_full_path"], verify_integrity=True)

    for bids_full_path, metadata in scans.iterrows():
        metadata.rename({"bids_filename": "filename"}, inplace=True)
        bids_full_path = Path(bids_full_path)
        try:
            os.makedirs(to / bids_full_path.parent)
        except OSError:
            pass
        dcm2niix_success = run_dcm2niix(
            Path(metadata["source_path"]).parent,
            to / Path(bids_full_path).parent,
            metadata["filename"],
            True,
        )
        if dcm2niix_success:
            scans_filepath = (
                to
                / str(metadata.participant_id)
                / str(metadata.session_id)
                / f"{metadata.participant_id}_{metadata.session_id}_scans.tsv"
            )
            metadata.source_path = metadata.source_path.relative_to(source)
            row_to_write = _serialize_row(
                metadata.drop(["participant_id", "session_id"]),
                write_column_names=not scans_filepath.exists(),
            )
            with open(scans_filepath, "a") as scans_file:
                scans_file.write(f"{row_to_write}\n")
            if "dwi" in metadata["filename"] and "Philips" in metadata.manufacturer:
                _merge_philips_diffusion(
                    to / bids_full_path.with_suffix(".json"),
                    metadata.number_of_parts,
                    metadata.run_num,
                )
    _correct_fieldmaps_name(to)
    _delete_real_and_imaginary_files(to)


def _serialize_row(row: pd.Series, write_column_names: bool) -> str:
    row_dict = row.to_dict()
    to_write = (
        [row_dict.keys(), row_dict.values()]
        if write_column_names
        else [row_dict.values()]
    )
    return "\n".join([_serialize_list(list(_)) for _ in to_write])


def _serialize_list(data: list, sep="\t") -> str:
    return sep.join([str(value) for value in data])


def _correct_fieldmaps_name(to: Path) -> None:
    """This function scans the BIDS after it has been written to correct the nameds of the fieldmap files.

    Parameters
    ----------
    to: Path
        The path to the BIDS

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
    >>> _correct_fieldmaps_name("path/to/bids")
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

    for z in to.rglob("*magnitude_e*"):
        os.rename(z, z.parent / re.sub(r"magnitude_e", "magnitude", z.name))

    for z in to.rglob("*phasediff_e*_ph*"):
        os.rename(z, z.parent / re.sub(r"phasediff_e[1-9]_ph", "phasediff", z.name))


def _merge_philips_diffusion(
    json_file: Path,
    number_of_parts: float,
    run_num: str,
) -> None:
    """Add the dwi number in the provided JSON file for each run of Philips images.
    The 'MultipartID' field of the input JSON file is set to 'dwi_1' or 'dwi_2' depending
    on the run number.
    """
    import json

    from clinica.utils.stream import cprint

    multipart_id = _get_multipart_id(
        PhilipsNumberOfParts.from_int(int(number_of_parts)), run_num
    )
    try:
        with open(json_file, "r+") as f:
            if multipart_id is not None:
                data = json.load(f)
                data["MultipartID"] = multipart_id
                json.dump(data, f, indent=4)
    except FileNotFoundError:
        cprint(msg=f"the file {json_file} does not exist.", lvl="warning")


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


def _delete_real_and_imaginary_files(bids_folder: Path) -> None:
    """Delete all files with a `_real` or `_imaginary` suffix from the BIDS folder.
    These files are generated by dcm2niix when it is not able to correctly convert the full image.
    These files break the BIDS specification, and are therefore removed.

    Parameters
    ----------
    bids_folder : PathLike
        Path to the BIDS folder.
    """
    from clinica.utils.stream import cprint

    for file_type in ("real", "imaginary"):
        for f in bids_folder.rglob(f"*_{file_type}.*"):
            f.unlink()
            cprint(f"file {f} was deleted as it is not BIDS compliant.")
