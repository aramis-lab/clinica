"""Methods to find information in the different pipelines of Clinica."""


import os
from glob import glob
from os import path
from typing import Tuple

import pandas as pd


def pet_volume_pipeline(
    caps_dir,
    df,
    group_selection=None,
    volume_atlas_selection=None,
    pvc_restriction=None,
    tracers_selection=None,
    **kwargs,
):
    """Merge the data of the PET-Volume pipeline to the merged file containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        group_selection: allows to choose the DARTEL groups to merge. If None all groups are selected.
        volume_atlas_selection: allows to choose the atlas to merge (default = 'all')
        pvc_restriction: gives the restriction on the inclusion or not of the file with the label 'pvc-rbv'
            1       --> only the atlases containing the label will be used
            0       --> the atlases containing the label won't be used
            None    --> all the atlases will be used
        tracers_selection: allows to choose the PET tracer to merge (default = 'all')

    Returns:
         final_df: a DataFrame containing the information of the bids and the pipeline
    """
    pet_path = path.join("pet", "preprocessing")

    return volume_pipeline(
        caps_dir,
        df,
        pet_path,
        group_selection=group_selection,
        atlas_selection=volume_atlas_selection,
        pvc_restriction=pvc_restriction,
        pipeline_name="pet-volume",
        tracers_selection=tracers_selection,
    )


def dwi_dti_pipeline(
    caps_dir: str, df: pd.DataFrame, dti_atlas_selection=None, **kwargs
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    from pathlib import Path

    # Ensures that df is correctly indexed
    try:
        df.set_index(
            ["participant_id", "session_id"], inplace=True, verify_integrity=True
        )
    except KeyError:
        raise KeyError("Fields `participant_id` and `session_id` are required.")

    subjects_dir = Path(caps_dir) / "subjects"

    pipeline_df = pd.DataFrame()

    metrics = ["FA", "MD", "RD", "AD"]

    for participant_id, session_id in df.index.values:

        ses_path = subjects_dir / participant_id / session_id
        mod_path = ses_path / "dwi" / "dti_based_processing" / "atlas_statistics"

        ses_df = pd.DataFrame(
            [[participant_id, session_id]], columns=["participant_id", "session_id"]
        )
        ses_df.set_index(["participant_id", "session_id"], inplace=True, drop=True)

        if mod_path.exists():

            for metric in metrics:

                atlas_paths = sorted(
                    (mod_path).glob(
                        f"{participant_id}_{session_id}_*{metric}_statistics.tsv"
                    )
                )
                for atlas_path in atlas_paths:

                    atlas_name = atlas_path.stem.split("_dwi_space-")[1].split("_")[0]
                    if atlas_path.exists() and (
                        not (
                            dti_atlas_selection
                            or (
                                dti_atlas_selection
                                and atlas_name in dti_atlas_selection
                            )
                        )
                    ):
                        atlas_df = pd.read_csv(atlas_path, sep="\t")
                        label_list = [
                            "dwi-dti_" + metric + "_atlas-" + atlas_name + "_" + x
                            for x in atlas_df.label_name.values
                        ]
                        ses_df[label_list] = atlas_df["mean_scalar"].to_numpy()

        pipeline_df = pd.concat([pipeline_df, ses_df])

        summary_df = generate_summary(pipeline_df, "dwi-dti", ignore_groups=True)

    final_df = pd.concat([df, pipeline_df], axis=1)
    final_df.reset_index(inplace=True)
    return final_df, summary_df


def t1_freesurfer_longitudinal_pipeline(
    caps_dir: str, df: pd.DataFrame, freesurfer_atlas_selection=None, **kwargs
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Merge the data of the T1FreesurferLongitudinal pipeline to the merged file containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        freesurfer_atlas_selection: allows to choose the atlas to merge (default = 'all')

    Returns:
         final_df: a DataFrame containing the information of the bids and the pipeline
    """

    from pathlib import Path

    from clinica.utils.stream import cprint

    # Ensures that df is correctly indexed
    try:
        df.set_index(
            ["participant_id", "session_id"], inplace=True, verify_integrity=True
        )
    except KeyError:
        raise KeyError("Fields `participant_id` and `session_id` are required.")

    subjects_dir = Path(caps_dir) / "subjects"

    pipeline_df = pd.DataFrame()

    metrics = ["volume", "thickness"]

    for participant_id, session_id in df.index.values:

        ses_path = subjects_dir / participant_id / session_id
        mod_path = ses_path / "t1"

        try:
            long_id = next(mod_path.glob("long*")).name
        except StopIteration:
            cprint(
                f"Could not find a longitudinal dataset for participant {participant_id} {session_id}",
                lvl="warning",
            )
            continue

        mod_path = mod_path / long_id / "freesurfer_longitudinal" / "regional_measures"

        ses_df = pd.DataFrame(
            [[participant_id, session_id]], columns=["participant_id", "session_id"]
        )
        ses_df.set_index(["participant_id", "session_id"], inplace=True, drop=True)

        if mod_path.exists():

            for metric in metrics:

                atlas_paths = sorted(
                    (mod_path).glob(f"{participant_id}_{session_id}_*{metric}.tsv")
                )

                for atlas_path in atlas_paths:
                    if "-wm_" not in str(atlas_path) and "-ba_" not in str(
                        atlas_path.stem
                    ):
                        atlas_name = atlas_path.stem.split("_parcellation-")[1].split(
                            "_"
                        )[0]
                        if atlas_path.exists() and (
                            not freesurfer_atlas_selection
                            or (
                                freesurfer_atlas_selection
                                and atlas_name in freesurfer_atlas_selection
                            )
                        ):
                            atlas_df = pd.read_csv(atlas_path, sep="\t")
                            label_list = [
                                "t1-fs-long_"
                                + metric
                                + "_atlas-"
                                + atlas_name
                                + "_"
                                + x
                                for x in atlas_df.label_name.values
                            ]
                            ses_df[label_list] = atlas_df["label_value"].to_numpy()

            # Always retrieve subcortical volumes

            try:
                atlas_path = next(
                    (mod_path).glob(
                        f"{participant_id}_{session_id}_*segmentationVolumes.tsv"
                    )
                )
                atlas_df = pd.read_csv(atlas_path, sep="\t")
                label_list = [
                    "t1-fs-long_segmentationVolumes_" + x
                    for x in atlas_df.label_name.values
                ]

                ses_df[label_list] = atlas_df["label_value"].to_numpy()

            except StopIteration:
                cprint("Segmentation volume file not found. Skipping.", lvl="warning")

        pipeline_df = pd.concat([pipeline_df, ses_df])

    summary_df = generate_summary(
        pipeline_df, "t1-freesurfer-longitudinal", ignore_groups=True
    )

    final_df = pd.concat([df, pipeline_df], axis=1)
    final_df.reset_index(inplace=True)

    return final_df, summary_df


def t1_freesurfer_pipeline(caps_dir, df, freesurfer_atlas_selection=None, **kwargs):
    """Merge the data of the T1Freesurfer pipeline to the merged file containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        freesurfer_atlas_selection: allows to choose the atlas to merge (default = 'all')

    Returns:
         final_df: a DataFrame containing the information of the bids and the pipeline
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
    )

    # Ensures that df is correctly indexed
    try:
        df.set_index(
            ["participant_id", "session_id"], inplace=True, verify_integrity=True
        )
    except KeyError:
        raise KeyError("Fields `participant_id` and `session_id` are required.")

    subjects_dir = path.join(caps_dir, "subjects")

    pipeline_df = pd.DataFrame()

    for participant_id, session_id in df.index.values:

        ses_path = path.join(subjects_dir, participant_id, session_id)
        mod_path = path.join(
            ses_path, "t1", "freesurfer_cross_sectional", "regional_measures"
        )
        ses_df = pd.DataFrame(
            [[participant_id, session_id]], columns=["participant_id", "session_id"]
        )
        ses_df.set_index(["participant_id", "session_id"], inplace=True, drop=True)

        if os.path.exists(mod_path):
            # Looking for atlases
            atlas_paths = sorted(
                glob(
                    path.join(mod_path, f"{participant_id}_{session_id}_*thickness.tsv")
                )
            )
            for atlas_path in atlas_paths:
                atlas_name = atlas_path.split("_parcellation-")[1].split("_")[0]
                if path.exists(atlas_path) and (
                    not freesurfer_atlas_selection
                    or (
                        freesurfer_atlas_selection
                        and atlas_name in freesurfer_atlas_selection
                    )
                ):
                    atlas_df = pd.read_csv(atlas_path, sep="\t")
                    label_list = [
                        f"t1-freesurfer_atlas-{atlas_name}_ROI-{replace_sequence_chars(roi_name)}_thickness"
                        for roi_name in atlas_df.label_name.values
                    ]
                    ses_df[label_list] = atlas_df["label_value"].to_numpy()

            # Always retrieve subcortical volumes
            atlas_path = path.join(
                mod_path, f"{participant_id}_{session_id}_segmentationVolumes.tsv"
            )
            atlas_df = pd.read_csv(atlas_path, sep="\t")
            label_list = [
                f"t1-freesurfer_segmentation-volumes_ROI-{replace_sequence_chars(roi_name)}_volume"
                for roi_name in atlas_df.label_name.values
            ]
            ses_df[label_list] = atlas_df["label_value"].to_numpy()

        pipeline_df = pd.concat([pipeline_df, ses_df])

    summary_df = generate_summary(pipeline_df, "t1-freesurfer", ignore_groups=True)
    final_df = pd.concat([df, pipeline_df], axis=1)
    final_df.reset_index(inplace=True)

    return final_df, summary_df


def t1_volume_pipeline(
    caps_dir, df, group_selection=None, volume_atlas_selection=None, **kwargs
):
    """Merge data of the t1-volume pipeline to the merged file containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        group_selection: allows to choose the DARTEL groups to merge. If None all groups are selected.
        volume_atlas_selection: allows to choose the atlas to merge. If None all atlases are selected.

    Returns:
        final_df: a DataFrame containing the information of the bids and the pipeline
    """
    t1_spm_path = path.join("t1", "spm", "dartel")

    return volume_pipeline(
        caps_dir,
        df,
        t1_spm_path,
        group_selection=group_selection,
        atlas_selection=volume_atlas_selection,
        pvc_restriction=None,
        pipeline_name="t1-volume",
    )


def volume_pipeline(
    caps_dir,
    df,
    pipeline_path,
    pipeline_name,
    group_selection=None,
    atlas_selection=None,
    pvc_restriction=None,
    tracers_selection=None,
):
    """Merge data of the t1-volume and pet-volume pipelines to the merged file containing the BIDS information.

    Args:
        caps_dir: the path to the CAPS directory
        df: the DataFrame containing the BIDS information
        pipeline_path: path between the session folder and the group folder
        pipeline_name: name of the pipeline
        group_selection: allows to choose the DARTEL groups to merge. If None all groups are selected.
        atlas_selection: allows to choose the atlas to merge. If None all atlases are selected.
        pvc_restriction: gives the restriction on the inclusion or not of the file with the label 'pvc-rbv'
            1       --> only the atlases containing the label will be used
            0       --> the atlases containing the label won't be used
            None    --> all the atlases will be used
        tracers_selection: allows to choose the PET tracer to merge (default = 'all')

    Returns:
        final_df: a DataFrame containing the information of the bids and the pipeline
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
    )

    # Ensures that df is correctly indexed
    try:
        df.set_index(
            ["participant_id", "session_id"], inplace=True, verify_integrity=True
        )
    except KeyError:
        raise KeyError("Fields `participant_id` and `session_id` are required.")

    if not group_selection:
        try:
            group_selection = os.listdir(path.join(caps_dir, "groups"))
        except FileNotFoundError:
            return df, None
    else:
        group_selection = [f"group-{group}" for group in group_selection]

    subjects_dir = path.join(caps_dir, "subjects")

    pipeline_df = pd.DataFrame()

    for participant_id, session_id in df.index.values:

        ses_path = path.join(subjects_dir, participant_id, session_id)
        mod_path = path.join(ses_path, pipeline_path)
        ses_df = pd.DataFrame(
            [[participant_id, session_id]], columns=["participant_id", "session_id"]
        )

        if os.path.exists(mod_path):
            # Looking for groups
            for group in group_selection:
                group_path = path.join(mod_path, group)
                if os.path.exists(group_path):
                    # Looking for atlases
                    if not atlas_selection:
                        atlas_paths = sorted(
                            glob(
                                path.join(
                                    group_path,
                                    "atlas_statistics",
                                    f"{participant_id}_{session_id}_*_statistics.tsv",
                                )
                            )
                        )
                    else:
                        atlas_paths = list()
                        for atlas in atlas_selection:
                            atlas_paths += sorted(
                                glob(
                                    path.join(
                                        group_path,
                                        "atlas_statistics",
                                        f"{participant_id}_{session_id}_*{atlas}*_statistics.tsv",
                                    )
                                )
                            )

                    # Filter pvc_restriction
                    if pvc_restriction is not None:
                        if pvc_restriction == 1:
                            atlas_paths = [
                                atlas_path
                                for atlas_path in atlas_paths
                                if "pvc-rbv" in atlas_path
                            ]
                        else:
                            atlas_paths = [
                                atlas_path
                                for atlas_path in atlas_paths
                                if "pvc-rbv" not in atlas_path
                            ]

                    # Filter tracers
                    if tracers_selection:
                        atlas_paths = [
                            atlas_path
                            for atlas_path in atlas_paths
                            for tracer in tracers_selection
                            if tracer in atlas_path
                        ]

                    atlas_df_list = list()
                    for atlas_path in atlas_paths:
                        atlas_name = atlas_path.split("_space-")[-1].split("_")[0]
                        if path.exists(atlas_path):
                            atlas_df = pd.read_csv(atlas_path, sep="\t")
                            additional_desc = ""
                            if "trc" in atlas_path:
                                tracer = atlas_path.split("_trc-")[1].split("_")[0]
                                additional_desc += f"_trc-{tracer}"
                            if "pvc-rbv" in atlas_path:
                                additional_desc += f"_pvc-rbv"

                            label_list = [
                                f"{pipeline_name}_{group}_atlas-{atlas_name}{additional_desc}_ROI-{replace_sequence_chars(roi_name)}_intensity"
                                for roi_name in atlas_df.label_name.values
                            ]
                            atlas_df_list.append(
                                pd.DataFrame(
                                    atlas_df["mean_scalar"].to_numpy().reshape(1, -1),
                                    columns=label_list,
                                )
                            )
                    ses_df = pd.concat([ses_df, *atlas_df_list], axis=1)

        pipeline_df = pd.concat([pipeline_df, ses_df])

    pipeline_df.set_index(["participant_id", "session_id"], inplace=True, drop=True)
    summary_df = generate_summary(pipeline_df, pipeline_name)
    final_df = pd.concat([df, pipeline_df], axis=1)
    final_df.reset_index(inplace=True)

    return final_df, summary_df


def generate_summary(pipeline_df, pipeline_name, ignore_groups=False):

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
    pvc_rectrictions = list(
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
        pvc_rectrictions = ["_"]
    groups.sort()
    atlases.sort()

    for group in groups:
        group_id = group.split("-")[-1]
        for atlas in atlases:
            atlas_id = atlas.split("-")[-1]
            for tracer in tracers:
                for pvc_rectriction in pvc_rectrictions:
                    if pvc_rectriction and pvc_rectriction != "_":
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
                                    pvc_rectriction,
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
