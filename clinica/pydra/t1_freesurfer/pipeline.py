import os
import typing as ty

import pydra


def build_workflow(
    name: str,
    bids_dir: os.PathLike,
    caps_dir: os.PathLike,
    parameters: ty.Optional[dict] = None,
) -> pydra.Workflow:
    from pydra.tasks import bids
    from pydra.tasks import freesurfer

    from .tasks import compute_freesurfer_subject_id

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("bids_dir", str, {"mandatory": True}),
            ("caps_dir", str, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    workflow = pydra.Workflow(
        input_spec=input_spec,
        bids_dir=bids_dir,
        caps_dir=caps_dir,
        name=name,
    )

    workflow.add(
        bids.read_bids_dataset(
            output_queries={
                "t1w": {"suffix": "T1w", "extension": ["nii", "nii.gz"]}
            },
            dataset_path = workflow.lzin.bids_dir,
            name = "read_t1w_files",
        )
    )

    workflow.add(
        bids.parse_bids_name(
            output_entities={
                "participant_id": "sub",
                "session_id": "ses",
            },
            file_path=workflow.read_t1w_files.lzout.t1w,
            name="parse_t1w_file",
        ).split("file_path")
    )

    workflow.add(
        compute_freesurfer_subject_id(
            name="compute_freesurfer_subject_id",
            participant_id=workflow.parse_t1w_file.lzout.participant_id,
            session_id=workflow.parse_t1w_file.lzout.session_id,
        )
    )

    workflow.add(
        freesurfer.ReconAll(
            name="recon_all",
            t1_volume_file=workflow.read_t1w_files.lzout.t1w,
            subject_id=workflow.compute_freesurfer_subject_id.lzout.subject_id,
            subjects_dir=workflow.lzin.caps_dir,
        ).split("t1_volume_file")
    )

    workflow.set_output(
        [
            ("subjects_dir", workflow.recon_all.lzout.subjects_dir),
            ("subject_id", workflow.recon_all.lzout.subject_id),
        ]
    )

    return workflow
