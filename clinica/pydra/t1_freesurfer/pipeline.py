import pydra

from clinica.pydra.engine import clinica_io


@clinica_io
def build_core_workflow(name: str, parameters: dict) -> pydra.Workflow:
    from pydra.tasks.bids import BIDSFileInfo
    from pydra.tasks.freesurfer import ReconAll

    from .tasks import compute_freesurfer_subject_id

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[("T1w", str, {"mandatory": True})],
        bases=(pydra.specs.BaseSpec,),
    )

    workflow = pydra.Workflow(name=name, input_spec=input_spec)

    workflow.add(
        BIDSFileInfo(
            output_entities={
                "participant_id": "sub",
                "session_id": "ses",
            }
        ).to_task(
            name="bids_file_info",
            file_path=workflow.lzin.T1w,
        ).split("file_path")
    )

    workflow.add(
        compute_freesurfer_subject_id(
            name="compute_freesurfer_subject_id",
            participant_id=workflow.bids_file_info.lzout.participant_id,
            session_id=workflow.bids_file_info.lzout.session_id,
        )
    )

    workflow.add(
        ReconAll(
            name="recon_all",
            t1_volume_file=workflow.lzin.T1w,
            subject_id=workflow.compute_freesurfer_subject_id.lzout.subject_id,
            subjects_dir=".",
        )
    )

    workflow.set_output(
        [
            ("subjects_dir", workflow.recon_all.lzout.subjects_dir),
            ("subject_id", workflow.recon_all.lzout.subject_id),
        ]
    )

    return workflow
