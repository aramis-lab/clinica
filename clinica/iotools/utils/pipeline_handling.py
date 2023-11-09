"""Methods to find information in the different pipelines of Clinica."""
import functools
import os
from os import PathLike
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


def _get_atlas_name(atlas_path: Path, pipeline: str) -> str:
    """Helper function for _extract_metrics_from_pipeline."""
    if pipeline == "dwi_dti":
        splitter = "_dwi_space-"
    elif pipeline in ("t1_freesurfer_longitudinal", "t1_freesurfer"):
        splitter = "_parcellation-"
    elif pipeline in ("t1-volume", "pet-volume"):
        splitter = "_space-"
    else:
        raise ValueError(f"Not supported pipeline {pipeline}.")
    if splitter in atlas_path.stem:
        try:
            return atlas_path.stem.split(splitter)[-1].split("_")[0]
        except Exception:
            pass
    raise ValueError(
        f"Unable to infer the atlas name from {atlas_path} for pipeline {pipeline}."
    )


def _get_mod_path(ses_path: Path, pipeline: str) -> Optional[Path]:
    """Helper function for _extract_metrics_from_pipeline.
    Returns the path to the modality of interest depending
    on the pipeline considered.
    """
    if pipeline == "dwi_dti":
        return ses_path / "dwi" / "dti_based_processing" / "atlas_statistics"
    if pipeline == "t1_freesurfer_longitudinal":
        mod_path = ses_path / "t1"
        long_ids = list(mod_path.glob("long*"))
        if len(long_ids) == 0:
            return None
        return (
            mod_path
            / long_ids[0].name
            / "freesurfer_longitudinal"
            / "regional_measures"
        )
    if pipeline == "t1_freesurfer":
        return ses_path / "t1" / "freesurfer_cross_sectional" / "regional_measures"
    if pipeline == "t1-volume":
        return ses_path / "t1" / "spm" / "dartel"
    if pipeline == "pet-volume":
        return ses_path / "pet" / "preprocessing"
    raise ValueError(f"Not supported pipeline {pipeline}.")


def _get_label_list(
    atlas_path: Path, metric: str, pipeline: str, group: str
) -> List[str]:
    """Helper function for _extract_metrics_from_pipeline.
    Returns the list of labels to use in the session df depending on the
    pipeline, the atlas, and the metric considered.
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
    )

    atlas_name = _get_atlas_name(atlas_path, pipeline)
    atlas_df = pd.read_csv(atlas_path, sep="\t")
    if pipeline == "t1-freesurfer":
        return [
            f"t1-freesurfer_atlas-{atlas_name}_ROI-{replace_sequence_chars(roi_name)}_thickness"
            for roi_name in atlas_df.label_name.values
        ]
    if pipeline in ("t1-volume", "pet-volume"):
        additional_desc = ""
        if "trc" in str(atlas_path):
            tracer = str(atlas_path).split("_trc-")[1].split("_")[0]
            additional_desc += f"_trc-{tracer}"
        if "pvc-rbv" in str(atlas_path):
            additional_desc += f"_pvc-rbv"
        return [
            f"{pipeline}_{group}_atlas-{atlas_name}{additional_desc}_ROI-{replace_sequence_chars(roi_name)}_intensity"
            for roi_name in atlas_df.label_name.values
        ]
    if pipeline == "dwi_dti":
        prefix = "dwi-dti_"
        metric = metric.rstrip("_statistics")
    elif pipeline == "t1-freesurfer_longitudinal":
        prefix = "t1-fs-long_"
    else:
        raise ValueError(f"Not supported pipeline {pipeline}.")
    return [
        prefix + metric + "_atlas-" + atlas_name + "_" + x
        for x in atlas_df.label_name.values
    ]


def _skip_atlas(
    atlas_path: Path,
    pipeline: str,
    pvc_restriction: Optional[bool] = None,
    tracers_selection: Optional[List[str]] = None,
) -> bool:
    """Helper function for _extract_metrics_from_pipeline.
    Returns whether the atlas provided through its path should be
    skipped for the provided pipeline.
    """
    if pipeline == "t1-freesurfer_longitudinal":
        return "-wm_" in str(atlas_path) or "-ba_" in str(atlas_path.stem)
    if pipeline in ("t1-volume", "pet-volume"):
        skip = []
        if pvc_restriction is not None:
            if pvc_restriction:
                skip.append("pvc-rbv" not in str(atlas_path))
            else:
                skip.append("pvc-rbv" in str(atlas_path))
        if tracers_selection:
            skip.append(
                all([tracer not in str(atlas_path) for tracer in tracers_selection])
            )
        return any(skip)
    return False


def _extract_metrics_from_pipeline(
    caps_dir: PathLike,
    df: pd.DataFrame,
    metrics: List[str],
    pipeline: str,
    atlas_selection: Optional[List[str]] = None,
    group_selection: Optional[List[str]] = None,
    pvc_restriction: Optional[bool] = None,
    tracers_selection: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract and merge the data of the provided pipeline into the
    merged dataframe already containing the BIDS information.

    Parameters
    ----------
    caps_dir : PathLike
        Path to the CAPS directory.

    df : pd.DataFrame
        DataFrame containing the BIDS information.
        New data will be merged into this dataframe.

    metrics : list of str
        List of metrics to extract from the pipeline's caps folder.

    pipeline : str
        Name of the pipeline for which to extract information.

    atlas_selection : list of str, optional
        Allows to choose the atlas to merge.
        If None all atlases are selected.

    group_selection : list of str, optional
        Allows to choose the DARTEL groups to merge.
        If None all groups are selected.

    pvc_restriction : bool, optional
        If True, only the atlases containing the label 'pvc-rbv' will be used.
        If False, the atlases containing the label won't be used.
        If None, all the atlases will be used.

    tracers_selection : list of str, optional
        Allows to choose the PET tracer to merge.
        If None, all tracers available will be used.

    Returns
    -------
    final_df : pd.DataFrame
        DataFrame containing the information from the bids dataset as well
        as the information extracted from the caps pipeline's folder.

    summary_df : pd.DataFrame
        Summary DataFrame generated by function `generate_summary`.
    """
    from clinica.utils.stream import cprint

    caps_dir = Path(caps_dir)
    if df.index.names != ["participant_id", "session_id"]:
        try:
            df.set_index(
                ["participant_id", "session_id"], inplace=True, verify_integrity=True
            )
        except KeyError:
            raise KeyError("Fields `participant_id` and `session_id` are required.")

    if group_selection is None:
        try:
            group_selection = [f.name for f in (caps_dir / "groups").iterdir()]
        except FileNotFoundError:
            return df, None
    else:
        group_selection = [f"group-{group}" for group in group_selection]
    ignore_groups = group_selection == [""]

    subjects_dir = caps_dir / "subjects"
    records = []
    for participant_id, session_id in df.index.values:
        ses_path = subjects_dir / participant_id / session_id
        mod_path = _get_mod_path(ses_path, pipeline)
        records.append({"participant_id": participant_id, "session_id": session_id})
        if mod_path is None:
            cprint(
                f"Could not find a longitudinal dataset for participant {participant_id} {session_id}",
                lvl="warning",
            )
            continue
        if mod_path.exists():
            for group in group_selection:
                group_path = mod_path / group
                if group_path.exists():
                    for metric in metrics:
                        atlas_paths = sorted(
                            group_path.glob(
                                f"{participant_id}_{session_id}_*{metric}.tsv"
                            )
                        )
                        if len(atlas_paths) == 0:
                            atlas_paths = sorted(
                                (group_path / "atlas_statistics").glob(
                                    f"{participant_id}_{session_id}_*{metric}.tsv"
                                )
                            )
                        for atlas_path in atlas_paths:
                            if not _skip_atlas(
                                atlas_path, pipeline, pvc_restriction, tracers_selection
                            ):
                                atlas_name = _get_atlas_name(atlas_path, pipeline)
                                if atlas_path.exists():
                                    if not (
                                        atlas_selection
                                        or (
                                            atlas_selection
                                            and atlas_name in atlas_selection
                                        )
                                    ):
                                        atlas_df = pd.read_csv(atlas_path, sep="\t")
                                        label_list = _get_label_list(
                                            atlas_path, metric, pipeline, group
                                        )
                                        key = (
                                            "label_value"
                                            if "freesurfer" in pipeline
                                            else "mean_scalar"
                                        )
                                        values = atlas_df[key].to_numpy()
                                        for label, value in zip(label_list, values):
                                            records[-1][label] = value
    pipeline_df = pd.DataFrame.from_records(
        records, index=["participant_id", "session_id"]
    )
    summary_df = generate_summary(pipeline_df, pipeline, ignore_groups=ignore_groups)
    final_df = pd.concat([df, pipeline_df], axis=1)
    final_df.reset_index(inplace=True)

    return final_df, summary_df


dwi_dti_pipeline = functools.partial(
    _extract_metrics_from_pipeline,
    metrics=["FA_statistics", "MD_statistics", "RD_statistics", "AD_statistics"],
    pipeline="dwi_dti",
    group_selection=[""],
)


t1_freesurfer_longitudinal_pipeline = functools.partial(
    _extract_metrics_from_pipeline,
    metrics=["volume", "thickness", "segmentationVolumes"],
    pipeline="t1_freesurfer_longitudinal",
    group_selection=[""],
)

t1_freesurfer_pipeline = functools.partial(
    _extract_metrics_from_pipeline,
    metrics=["thickness", "segmentationVolumes"],
    pipeline="t1_freesurfer",
    group_selection=[""],
)

t1_volume_pipeline = functools.partial(
    _extract_metrics_from_pipeline,
    metrics=["statistics"],
    pipeline="t1-volume",
)

pet_volume_pipeline = functools.partial(
    _extract_metrics_from_pipeline,
    metrics=["statistics"],
    pipeline="pet-volume",
)


def generate_summary(
    pipeline_df: pd.DataFrame, pipeline_name: str, ignore_groups: bool = False
):

    columns = [
        "pipeline_name",
        "group_id",
        "atlas_id",
        "tracer",
        "pvc",
        "regions_number",
        "first_column_name",
        "last_column_name",
    ]
    summary_df = pd.DataFrame(columns=columns)

    if ignore_groups:
        groups = ["_"]
        atlases = list({column.split("_")[1] for column in pipeline_df.columns.values})
    else:
        groups = list({column.split("_")[1] for column in pipeline_df.columns.values})
        atlases = list({column.split("_")[2] for column in pipeline_df.columns.values})
    pvc_restrictions = list(
        {"pvc-rbv" in column for column in pipeline_df.columns.values}
    )
    tracers = list(
        {
            column.split("_trc-")[1].split("_")[0]
            for column in pipeline_df.columns.values
            if "trc" in column
        }
    )
    if len(tracers) == 0:
        tracers.append("_")
        pvc_restrictions = ["_"]
    groups.sort()
    atlases.sort()

    for group in groups:
        group_id = group.split("-")[-1]
        for atlas in atlases:
            atlas_id = atlas.split("-")[-1]
            for tracer in tracers:
                for pvc_restriction in pvc_restrictions:
                    if pvc_restriction and pvc_restriction != "_":
                        regions = [
                            column
                            for column in pipeline_df.columns.values
                            if group in column
                            and atlas in column
                            and "pvc-rbv" in column
                            and tracer in column
                        ]
                    else:
                        regions = [
                            column
                            for column in pipeline_df.columns.values
                            if group in column
                            and atlas in column
                            and "pvc-rbv" not in column
                            and tracer in column
                        ]

                    if len(regions) > 0:
                        row_df = pd.DataFrame(
                            [
                                [
                                    pipeline_name,
                                    group_id,
                                    atlas_id,
                                    tracer,
                                    pvc_restriction,
                                    len(regions),
                                    regions[0],
                                    regions[-1],
                                ]
                            ],
                            columns=columns,
                        )
                        summary_df = pd.concat([summary_df, row_df])

    summary_df = summary_df.replace("_", "n/a")
    return summary_df


class DatasetError(Exception):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return repr("Bad format for the sessions: " + self.name)
