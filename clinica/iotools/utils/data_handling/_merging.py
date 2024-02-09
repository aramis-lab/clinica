from os import PathLike
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


def create_merge_file(
    bids_dir: PathLike,
    out_tsv: PathLike,
    caps_dir: Optional[PathLike] = None,
    tsv_file: Optional[PathLike] = None,
    pipelines: Optional[List[str]] = None,
    ignore_scan_files: bool = False,
    ignore_sessions_files: bool = False,
    **kwargs,
):
    """Merge all the TSV files containing clinical data of a BIDS compliant
    dataset and store the result inside a TSV file.

    Parameters
    ----------
    bids_dir : PathLike
        Path to the BIDS folder.

    out_tsv : PathLike
        Path to the output tsv file.

    caps_dir : PathLike, optional
        Path to the CAPS folder.

    tsv_file : PathLike, optional
        Path to a TSV file containing the subjects with their sessions.

    pipelines : list of str, optional
        When adding CAPS information, indicates the pipelines that will be merged.

    ignore_scan_files : bool, optional
        If True the information related to scans is not read.

    ignore_sessions_files : bool, optional
        If True the information related to sessions and scans is not read.
    """
    from os import path

    from clinica.utils.stream import cprint

    bids_dir = Path(bids_dir)
    caps_dir = _validate_caps_dir(caps_dir)
    out_path = _validate_output_tsv_path(out_tsv)
    participants_df, sub_ses_df = _get_participants_and_subjects_sessions_df(
        bids_dir, tsv_file, ignore_sessions_files
    )
    merged_df = _create_merge_file_from_bids(
        bids_dir, sub_ses_df, participants_df, ignore_scan_files, ignore_sessions_files
    )
    merged_df.to_csv(out_path, sep="\t", index=False)
    cprint("End of BIDS information merge.", lvl="debug")
    merged_df.reset_index(drop=True, inplace=True)
    if caps_dir is not None:
        merged_df, merged_summary_df = _add_data_to_merge_file_from_caps(
            caps_dir, merged_df, pipelines, **kwargs
        )
        summary_path = path.splitext(str(out_path))[0] + "_summary.tsv"
        merged_summary_df.to_csv(summary_path, sep="\t", index=False)
        merged_df.to_csv(out_path, sep="\t", index=False)
        cprint("End of CAPS information merge.", lvl="debug")


def _validate_caps_dir(caps_dir: Optional[PathLike] = None) -> Optional[Path]:
    if caps_dir is not None:
        caps_dir = Path(caps_dir)
        if not caps_dir.is_dir():
            raise IOError("The path to the CAPS directory is wrong")
    return caps_dir


def _validate_output_tsv_path(out_path: PathLike) -> Path:
    """Validate that provided file path is a TSV file.

    If folders do not exist, this function will create them.
    If provided path is a directory, this function will return
    a file named 'merge.tsv' within this directory.
    """
    import warnings

    from clinica.utils.stream import cprint

    out_path = Path(out_path).resolve()
    if out_path.is_dir():
        out_path = out_path / "merge.tsv"
    elif "." not in out_path.name:
        out_path = out_path.with_suffix(".tsv")
    elif out_path.suffix != ".tsv":
        raise TypeError("Output path extension must be tsv.")
    if out_path.exists():
        msg = (
            f"Path to TSV file {out_path} already exists. The file will be overwritten."
        )
        warnings.warn(msg)
        cprint(msg=msg, lvl="warning")
    out_dir = out_path.parent
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    return out_path


def _get_participants_and_subjects_sessions_df(
    bids_dir: Path,
    tsv_file: Optional[PathLike] = None,
    ignore_sessions_files: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    from clinica.utils.participant import get_subject_session_list
    from clinica.utils.stream import cprint

    index_cols = ["participant_id", "session_id"]
    subjects, sessions = get_subject_session_list(
        bids_dir, ss_file=tsv_file, use_session_tsv=(not ignore_sessions_files)
    )
    if (bids_dir / "participants.tsv").is_file():
        participants_df = pd.read_csv(bids_dir / "participants.tsv", sep="\t")
    else:
        participants_df = pd.DataFrame(
            sorted(list(set(subjects))), columns=["participant_id"]
        )
    sub_ses_df = pd.DataFrame(
        [[subject, session] for subject, session in zip(subjects, sessions)],
        columns=index_cols,
    )
    try:
        sub_ses_df.set_index(index_cols, inplace=True, verify_integrity=True)
    except ValueError:
        cprint(
            "Found duplicate subject/session pair. Keeping first occurrence.",
            lvl="warning",
        )
        sub_ses_df = sub_ses_df.drop_duplicates(subset=index_cols)
        sub_ses_df.set_index(index_cols, inplace=True)

    return participants_df, sub_ses_df


def _create_merge_file_from_bids(
    bids_dir: Path,
    sub_ses_df: pd.DataFrame,
    participants_df: pd.DataFrame,
    ignore_scan_files: bool = False,
    ignore_sessions_files: bool = False,
) -> pd.DataFrame:
    """Create a merge pandas dataframe for a given BIDS dataset."""
    df = pd.DataFrame(columns=participants_df.columns.values)
    for subject, subject_df in sub_ses_df.groupby(level=0):
        row_subject_data_only = _get_row_data_from_participant_df(
            participants_df, subject
        )
        if ignore_sessions_files:
            rows_session_data_only = (
                _compute_session_rows_for_subject_without_session_files(subject_df)
            )
        else:
            rows_session_data_only = (
                _compute_session_rows_for_subject_with_session_files(
                    bids_dir, subject, subject_df, ignore_scan_files
                )
            )
        for row_session_data_only in rows_session_data_only:
            full_row = pd.concat([row_subject_data_only, row_session_data_only], axis=1)
            df = pd.concat([df, full_row], axis=0)
    return _post_process_merge_file_from_bids(df)


def _get_row_data_from_participant_df(
    participants_df: pd.DataFrame, subject: str
) -> pd.DataFrame:
    """Return, for a given subject, part of a merged dataframe row.
    This part of the full row only contains information found in the participants dataframe.
    In order to make a full row entry, it will be concatenated with information from the sessions and scans.
    """
    from clinica.utils.stream import cprint

    row_subject_data = participants_df[participants_df["participant_id"] == subject]
    row_subject_data.reset_index(inplace=True, drop=True)
    if len(row_subject_data) == 0:
        cprint(
            msg=f"Participant {subject} does not exist in participants.tsv",
            lvl="warning",
        )
        row_subject_data = pd.DataFrame([[subject]], columns=["participant_id"])

    return row_subject_data


def _compute_session_rows_for_subject_without_session_files(
    subject_df: pd.DataFrame,
) -> List[pd.DataFrame]:
    """Compute all rows with session information for a given subject.
    This will be concatenated with information from the participants TSV file
    to build a full dataframe row with all relevant information.
    This is a default implementation when sessions and scans are ignored.
    """
    return [
        pd.DataFrame([[session]], columns=["session_id"])
        for _, session in subject_df.index.values
    ]


def _compute_session_rows_for_subject_with_session_files(
    bids_dir: Path,
    subject: str,
    subject_df: pd.DataFrame,
    ignore_scan_files: bool,
) -> List[pd.DataFrame]:
    """Compute all rows with session information for a given subject.
    This will be concatenated with information from the participants TSV file
    to build a full dataframe row with all relevant information.
    This function will try to incorporate data from the sessions and scans
    if available.
    """
    from ..pipeline_handling import DatasetError

    rows = []
    sub_path = bids_dir / subject
    sessions_df = pd.read_csv(sub_path / f"{subject}_sessions.tsv", sep="\t")
    for _, session in subject_df.index.values:
        row_session_df = sessions_df[sessions_df.session_id == session]
        row_session_df.reset_index(inplace=True, drop=True)
        if len(row_session_df) == 0:
            raise DatasetError(sessions_df.loc[0, "session_id"] + " / " + session)
        scan_path = bids_dir / subject / session / f"{subject}_{session}_scans.tsv"
        row_scans_df = _get_row_scan_df(scan_path, ignore_scan_files)
        row_df = pd.concat([row_session_df, row_scans_df], axis=1)
        rows.append(row_df.loc[:, ~row_df.columns.duplicated(keep="last")])

    return rows


def _get_row_scan_df(scan_path: Path, ignore_scan_files: bool) -> pd.DataFrame:
    """Return a dataframe for all scans associated with a given subject and session.
    The scan data come from the scan TSV file. If it doesn't exist, an empty dataframe
    is returned.
    """
    row_scans_df = pd.DataFrame()
    if ignore_scan_files or not scan_path.is_file():
        return row_scans_df
    scans_dict = dict()
    scans_df = pd.read_csv(scan_path, sep="\t")
    for idx in scans_df.index.values:
        filepath = scans_df.loc[idx, "filename"]
        if not filepath.endswith(".nii.gz"):
            continue
        filename = Path(filepath).name.split(".")[0]
        modality = "_".join(filename.split("_")[2::])
        for col in scans_df.columns.values:
            if col != "filename":
                value = scans_df.loc[idx, col]
                new_col_name = f"{modality}_{col}"
                scans_dict.update({new_col_name: value})
        json_path = scan_path.parent / f"{filepath.split('.')[0]}.json"
        scans_dict = _add_metadata_from_json(json_path, scans_dict, modality)
    scans_dict = {str(key): str(value) for key, value in scans_dict.items()}
    row_scans_df = pd.DataFrame(scans_dict, index=[0])

    return row_scans_df


def _add_metadata_from_json(json_path: Path, scans_dict: dict, modality: str) -> dict:
    """Add metadata from the provided JSON sidecar file to the scan dictionary."""
    import json
    import warnings

    if json_path.exists():
        try:
            with open(json_path, "r") as f:
                json_dict = json.load(f)
        except json.JSONDecodeError as e:
            warnings.warn(
                f"Couldn't parse the JSON file {json_path}. Ignoring metadata.\nJson error is:\n{e}"
            )
            json_dict = {}
        for key, value in json_dict.items():
            new_col_name = f"{modality}_{key}"
            scans_dict.update({new_col_name: value})
    return scans_dict


def _post_process_merge_file_from_bids(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Put participant_id and session_id first, and round numbers."""
    col_list = merged_df.columns.values.tolist()
    col_list.insert(0, col_list.pop(col_list.index("participant_id")))
    col_list.insert(1, col_list.pop(col_list.index("session_id")))
    merged_df = merged_df[col_list]
    tmp = merged_df.select_dtypes(include=[np.number])
    merged_df.loc[:, tmp.columns] = np.round(tmp, 6)

    return merged_df


def _add_data_to_merge_file_from_caps(
    caps_dir: Path,
    merged_df: pd.DataFrame,
    pipelines: Optional[List[str]] = None,
    **kwargs,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    from clinica.utils.stream import cprint

    from ..pipeline_handling import (
        PipelineNameForMetricExtraction,
        pipeline_metric_extractor_factory,
    )

    merged_summary_df = pd.DataFrame()
    if "group_selection" in kwargs and kwargs["group_selection"] is None:
        kwargs.pop("group_selection")
    pipelines = pipelines or PipelineNameForMetricExtraction
    for pipeline in pipelines:
        metric_extractor = pipeline_metric_extractor_factory(pipeline)
        cprint(f"Extracting from CAPS pipeline output: {pipeline}...", lvl="info")
        merged_df, summary_df = metric_extractor(caps_dir, merged_df, **kwargs)
        if summary_df is not None and not summary_df.empty:
            merged_summary_df = pd.concat([merged_summary_df, summary_df])
        if summary_df is None or summary_df.empty:
            cprint(f"{pipeline} outputs were not found in the CAPS folder.", lvl="info")

    return _post_process_merged_df_from_caps(merged_df, merged_summary_df)


def _post_process_merged_df_from_caps(
    merged_df: pd.DataFrame,
    merged_summary_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if len(merged_summary_df) == 0:
        raise FileNotFoundError(
            "No outputs were found for any pipeline in the CAPS folder. "
            "The output only contains BIDS information."
        )
    columns = merged_df.columns.values.tolist()
    merged_summary_df.reset_index(inplace=True, drop=True)
    for idx in merged_summary_df.index:
        first_column_name = merged_summary_df.loc[idx, "first_column_name"]
        last_column_name = merged_summary_df.loc[idx, "last_column_name"]
        merged_summary_df.loc[idx, "first_column_index"] = columns.index(
            first_column_name
        )
        merged_summary_df.loc[idx, "last_column_index"] = columns.index(
            last_column_name
        )
    tmp = merged_df.select_dtypes(include=[np.number])
    merged_df.loc[:, tmp.columns] = np.round(tmp, 12)

    return merged_df, merged_summary_df
