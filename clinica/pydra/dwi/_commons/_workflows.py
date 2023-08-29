from pathlib import Path

from pydra import Workflow
from pydra.specs import SpecInfo


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
        "brain_mask": wf.mask_reference_b0.brain_mask
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
            "parameters_file": wf.eddy_fsl.lzout.parameters_file,
            "corrected_image": wf.eddy_fsl.lzout.corrected_image,
            "rotated_bvec_file": wf.eddy_fsl.lzout.rotated_bvec_file,
        }
    )

    return wf


def _build_eddy_fsl_input_specs(
    image_id: bool,
    field: bool,
    compute_mask: bool,
    name: str = "InputSpec",
) -> SpecInfo:
    fields = [
        ("dwi_filename", Path),
        ("b_vectors_filename", Path),
        ("b_values_filename", Path),
        ("in_mask", Path),
        ("total_readout_time", str),
        ("phase_encoding_direction", str),
    ]
    if image_id:
        fields.append(("image_id", str))
    if field:
        fields.append(("field", Path))
    if compute_mask:
        fields.append(("reference_b0", Path))

    return SpecInfo(name=name, fields=fields)
