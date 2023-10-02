from typing import List, Optional

from pydra import Workflow


def build_eddy_fsl_workflow(
    compute_mask: bool = True,
    image_id: bool = False,
    field: bool = False,
    output_dir: Optional[str] = None,
    name: str = "eddy_fsl",
) -> Workflow:
    """Build a FSL-based pipeline for head motion correction and eddy current distortion correction.

    This pipeline needs FSL as it relies on the Eddy executable.

    Parameters
    ----------
    compute_mask : bool, optional
        If True, the eddy_fsl pipeline computes the brain mask
        of the reference B0 volume using FSL BET. And provide the
        result to Eddy.
        If False, a mask should be provided to the eddy_fsl
        pipeline through the "in_mask" connection.
        Default=True.

    image_id : bool, optional
        Boolean to indicate whether the Eddy node should expect
        a value for its 'out_base' parameter. This is used for
        building the names of the output files.
        Default=False.

    field : bool, optional
        Boolean to indicate whether the Eddy node should expect
        a value for its 'field' input.
        Default=False.

    output_dir: str, optional
        Path to output directory.
        If provided, the pipeline will write its output in this folder.
        Default to None.
        TODO: Implement once we can actually write things.

    name : str, optional
        Name of the pipeline. Default='eddy_fsl'.

    Returns
    -------
    Workflow :
        The Pydra workflow.
        This workflow has the following inputs:
            - "dwi_filename": The path to the DWI image
            - "b_vectors_filename": The path to the associated B-vectors file
            - "b_values_filename": The path to the associated B-values file
            - "in_mask": The path to the mask image to be provided to Eddy
              This will be used only if compute_mask is False.
            - "image_id": Prefix to be used for output files
            - "field": The path to the field image to be used by Eddy
            - "reference_b0": The path to the reference B0 image
              This will be used to compute the mask only if compute_mask is True.
            - "total_readout_time": The total readout time extracted from JSON metadata
            - "phase_encoding_direction": The phase encoding direction extracted from JSON metadata

        And the following outputs:
            - "out_parameter": Path to the file storing the output parameters
            - "out_corrected": Path to the corrected image
            - "out_rotated_bvecs": Path to the file holding the rotated B-vectors
    """
    from pydra.tasks import fsl

    from ._tasks import generate_acq_file_task, generate_index_file_task

    wf = Workflow(
        name, input_spec=_build_eddy_fsl_input_specs(image_id, field, compute_mask)
    )

    if compute_mask:
        wf.add(
            fsl.BET(
                name="mask_reference_b0",
                fractional_intensity_threshold=0.3,
                save_brain_mask=True,
                with_robust_brain_center_estimation=True,
                input_image=wf.lzin.reference_b0,
            )
        )

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

    wf.add(
        fsl.Eddy(
            name="eddy_fsl",
            replace_outliers=True,
            bvec_file=wf.lzin.b_vectors_filename,
            bval_file=wf.lzin.b_values_filename,
            input_image=wf.lzin.dwi_filename,
            encoding_file=wf.generate_acq_file_task.lzout.out_file,
            index_file=wf.generate_index_file_task.lzout.out_file,
            brain_mask=wf.mask_reference_b0.lzout.brain_mask
            if compute_mask
            else wf.lzin.in_mask,
            output_basename=wf.lzin.image_id if image_id else "eddy",
        )
    )
    if field:
        wf.eddy_fsl.inputs.fieldmap_image = wf.lzin.field

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
        "total_readout_time",
        "phase_encoding_direction",
    ]
    if image_id:
        fields.append("image_id")
    if field:
        fields.append("field")
    fields.append("reference_b0" if compute_mask else "in_mask")

    return fields
