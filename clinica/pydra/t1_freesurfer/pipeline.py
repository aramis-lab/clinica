import os
import typing as ty

import pydra


def input_workflow(name: str = "input_workflow", **inputs) -> pydra.Workflow:
    from pydra.tasks import bids

    workflow = pydra.Workflow(
        name=name,
        input_spec=["bids_dir"],
        **inputs,
    )

    workflow.add(
        bids.read_bids_dataset(
            output_queries={
                "out": {
                    "suffix": "T1w",
                    "extension": ["nii", "nii.gz"],
                }
            },
            dataset_path=workflow.lzin.bids_dir,
            name="read_t1w",
        )
    )

    workflow.set_output([("t1w", workflow.read_t1w.lzout.out)])

    return workflow


def core_workflow(name: str = "core_workflow", **inputs) -> pydra.Workflow:
    from pydra.tasks import bids, freesurfer

    from .tasks import compute_freesurfer_subject_id

    workflow = pydra.Workflow(
        name=name,
        input_spec=["caps_dir", "t1w"],
        **inputs,
    )

    workflow.add(
        bids.parse_bids_name(
            output_entities={
                "participant_id": "sub",
                "session_id": "ses",
            },
            file_path=workflow.lzin.t1w,
            name="parse_t1w",
        )
    )

    workflow.add(
        compute_freesurfer_subject_id(
            participant_id=workflow.parse_t1w.lzout.participant_id,
            session_id=workflow.parse_t1w.lzout.session_id,
            name="compute_freesurfer_subject_id",
        )
    )

    workflow.add(
        freesurfer.ReconAll(
            name="recon_all",
            t1_volume_file=workflow.lzin.t1w,
            subject_id=workflow.compute_freesurfer_subject_id.lzout.subject_id,
            subjects_dir=workflow.lzin.caps_dir,
        )
    )

    workflow.set_output(
        [
            ("subjects_dir", workflow.recon_all.lzout.subjects_dir),
            ("subject_id", workflow.recon_all.lzout.subject_id),
        ]
    )

    return workflow


def build_workflow(
    bids_dir: os.PathLike,
    caps_dir: os.PathLike,
    atlas_path: os.PathLike,
    name: str = "t1_freesurfer",
) -> pydra.Workflow:

    workflow = pydra.Workflow(
        name=name,
        input_spec=["bids_dir", "caps_dir"],
        bids_dir=bids_dir,
    )

    workflow.add(
        input_workflow(
            name="input_workflow",
            bids_dir=workflow.lzin.bids_dir,
        )
    )

    workflow.add(
        core_workflow(
            name="core_workflow",
            caps_dir=caps_dir,
            t1w=workflow.input_workflow.lzout.t1w,
        ).split("t1w")
    )

    workflow.set_output(
        [
            ("t1w", workflow.input_workflow.lzout.t1w),
            ("subjects_dir", workflow.core_workflow.lzout.subjects_dir),
            ("subject_id", workflow.core_workflow.lzout.subject_id),
        ]
    )

    return workflow
