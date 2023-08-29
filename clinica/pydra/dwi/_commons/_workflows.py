from typing import List

from pydra import Workflow


def build_eddy_fsl_workflow(
    compute_mask: bool = True,
    image_id: bool = False,
    field: bool = False,
    output_dir=None,
    name: str = "eddy_fsl",
) -> Workflow:
    """TODO"""
    from pydra.tasks import fsl

    from ._tasks import generate_acq_file_task, generate_index_file_task

    wf = Workflow(
        name, input_spec=_build_eddy_fsl_input_specs(image_id, field, compute_mask)
    )

    bet_config = {
        "name": "mask_reference_b0",
        "fractional_intensity_threshold": 0.3,
        "save_brain_mask": True,
        "with_robust_brain_center_estimation": True,
    }
    if compute_mask:
        bet_config["input_image"] = wf.lzin.reference_b0

    wf.add(fsl.BET(**bet_config))

    wf.add(
        generate_acq_file_task(
            dwi_filename=wf.lzin.dwi_filename,
            fsl_phase_encoding_direction=wf.lzin.phase_encoding_direction,
            total_readout_time=wf.lzin.total_readout_time,
            image_id=wf.lzin.image_id if image_id else None,
        )
    )

    wf.add(
        generate_index_file_task(
            b_values_filename=wf.lzin.b_values_filename,
            image_id=wf.lzin.image_id if image_id else None,
        )
    )

    eddy_config = {
        "name": "eddy_fsl",
        "replace_outliers": True,
        "bvec_file": wf.lzin.b_vectors_filename,
        "bval_file": wf.lzin.b_values_filename,
        "input_image": wf.lzin.dwi_filename,
        "encoding_file": wf.generate_acq_file_task.lzout.out_file,
        "index_file": wf.generate_index_file_task.lzout.out_file,
        "brain_mask": wf.mask_reference_b0.lzout.brain_mask
        if compute_mask
        else wf.lzin.in_mask,
    }
    if image_id:
        eddy_config["output_basename"] = wf.lzin.image_id
    if field:
        eddy_config["fieldmap_image"] = wf.lzin.field

    wf.add(fsl.Eddy(**eddy_config))

    wf.set_output(
        {
            "corrected_image": wf.eddy_fsl.lzout.corrected_image,
            "parameters_file": wf.eddy_fsl.lzout.parameters_file,
            "rotated_bvec_file": wf.eddy_fsl.lzout.rotated_bvec_file,
        }
    )

    return wf


def _build_eddy_fsl_input_specs(
    image_id: bool,
    field: bool,
    compute_mask: bool,
) -> List[str]:
    """Returns the input specs for eddy_fsl workflow."""
    fields = [
        "dwi_filename",
        "b_vectors_filename",
        "b_values_filename",
        "in_mask",
        "total_readout_time",
        "phase_encoding_direction",
    ]
    if image_id:
        fields.append("image_id")
    if field:
        fields.append("field")
    fields.append("reference_b0" if compute_mask else "in_mask")

    return fields
