"""This module contains workflows used by multiple DWI preprocessing pipelines."""

from nipype.pipeline.engine import Workflow

__all__ = ["eddy_fsl_pipeline"]


def eddy_fsl_pipeline(
    base_dir: str,
    use_cuda: bool,
    initrand: bool,
    compute_mask: bool = True,
    image_id: bool = False,
    field: bool = False,
    output_dir=None,
    name: str = "eddy_fsl",
) -> Workflow:
    """Build a FSL-based pipeline for head motion correction and eddy current distortion correction.

    This pipeline needs FSL as it relies on the Eddy executable.

    Parameters
    ----------
    base_dir : str
        The base directory.

    use_cuda : bool
        Boolean to indicate whether cuda should be used or not.

    initrand : bool
        ????

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

    name : str, optional
        Name of the pipeline. Default='eddy_fsl'.

    Returns
    -------
    Workflow :
        The Nipype workflow.
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
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces.fsl import BET
    from nipype.interfaces.fsl.epi import Eddy

    from .tasks import generate_acq_file_task, generate_index_file_task

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "dwi_filename",
                "b_vectors_filename",
                "b_values_filename",
                "in_mask",
                "image_id",
                "field",
                "reference_b0",
                "total_readout_time",
                "phase_encoding_direction",
            ]
        ),
        name="inputnode",
    )

    mask_reference_b0 = pe.Node(
        BET(frac=0.3, mask=True, robust=True), name="mask_reference_b0"
    )

    generate_acq = pe.Node(
        niu.Function(
            input_names=[
                "dwi_filename",
                "fsl_phase_encoding_direction",
                "total_readout_time",
                "image_id",
                "output_dir",
            ],
            output_names=["out_file"],
            function=generate_acq_file_task,
        ),
        name="generate_acq",
    )
    generate_acq.inputs.output_dir = base_dir

    generate_index = pe.Node(
        niu.Function(
            input_names=["b_values_filename", "output_dir"],
            output_names=["out_file"],
            function=generate_index_file_task,
        ),
        name="generate_index",
    )
    generate_index.inputs.output_dir = base_dir

    eddy = pe.Node(interface=Eddy(), name="eddy_fsl")
    eddy.inputs.repol = True
    eddy.inputs.use_cuda = use_cuda
    eddy.inputs.initrand = initrand

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["out_parameter", "out_corrected", "out_rotated_bvecs"]
        ),
        name="outputnode",
    )

    if output_dir:
        write_results = pe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = output_dir
        write_results.inputs.parameterization = False

    wf = pe.Workflow(name=name)

    connections = [
        (inputnode, generate_acq, [("dwi_filename", "dwi_filename")]),
        (inputnode, generate_acq, [("total_readout_time", "total_readout_time")]),
        (
            inputnode,
            generate_acq,
            [("phase_encoding_direction", "fsl_phase_encoding_direction")],
        ),
        (inputnode, generate_index, [("b_values_filename", "b_values_filename")]),
        (inputnode, generate_index, [("image_id", "image_id")]),
        (inputnode, eddy, [("b_vectors_filename", "in_bvec")]),
        (inputnode, eddy, [("b_values_filename", "in_bval")]),
        (inputnode, eddy, [("dwi_filename", "in_file")]),
        (generate_acq, eddy, [("out_file", "in_acqp")]),
        (generate_index, eddy, [("out_file", "in_index")]),
        (eddy, outputnode, [("out_parameter", "out_parameter")]),
        (eddy, outputnode, [("out_corrected", "out_corrected")]),
        (eddy, outputnode, [("out_rotated_bvecs", "out_rotated_bvecs")]),
    ]

    if image_id:
        connections += [
            (inputnode, generate_acq, [("image_id", "image_id")]),
            (inputnode, eddy, [("image_id", "out_base")]),
        ]

    if field:
        connections += [(inputnode, eddy, [("field", "field")])]

    if compute_mask:
        connections += [
            (inputnode, mask_reference_b0, [("reference_b0", "in_file")]),
            (mask_reference_b0, eddy, [("mask_file", "in_mask")]),
        ]
    else:
        connections += [(inputnode, eddy, [("in_mask", "in_mask")])]

    if output_dir:
        connections += [
            (outputnode, write_results, [("out_parameter", "out_parameter")]),
            (outputnode, write_results, [("out_corrected", "out_corrected")]),
            (outputnode, write_results, [("out_rotated_bvecs", "out_rotated_bvecs")]),
        ]

    wf.connect(connections)

    return wf
