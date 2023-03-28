from typing import Optional

from nipype.pipeline.engine import Workflow


def prepare_phasediff_fmap(name="prepare_phasediff_fmap"):
    """This workflow adapts the fsl_prepare_fieldmap script from FSL for the FSL eddy command.

    Please note that the step 3 converts the fieldmap into Hz instead of rad/s
    in the initial script because FSL eddy --field expects the fieldmap to be in Hz.

    Input node:
        fmap_mask (str): Binary mask of the fieldmap.
        fmap_phasediff (str): Phase difference fieldmap.
        fmap_magnitude (str): Brain extracted magnitude fieldmap. Chose the fieldmap with the best contrast.
        delta_echo_time (float): DeltaEchoTime from BIDS specifications.

    Output node:
        calibrated_fmap (str): Calibrated fieldmap for eddy --field command.

    Warnings:
        This workflow can not be used for PRELUDE. You would need to change the step 3 (conversion into rad/s instead of Hz).
    """
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as nutil
    import nipype.pipeline.engine as npe
    from niflow.nipype1.workflows.dmri.fsl.utils import (
        cleanup_edge_pipeline,
        demean_image,
        siemens2rads,
    )

    from .dwi_preprocessing_using_phasediff_fmap_utils import rads2hz

    input_node = npe.Node(
        nutil.IdentityInterface(
            fields=["fmap_mask", "fmap_phasediff", "fmap_magnitude", "delta_echo_time"]
        ),
        name="input_node",
    )

    output_node = npe.Node(
        nutil.IdentityInterface(fields=["calibrated_fmap"]), name="output_node"
    )

    # Step 1 - Convert the fmap into radians
    pha2rads = npe.Node(
        nutil.Function(
            input_names=["in_file"], output_names=["out_file"], function=siemens2rads
        ),
        name="1-ConvertFMapToRads",
    )

    # Step 2 - Unwrap the fmap with PRELUDE
    prelude = npe.Node(fsl.PRELUDE(process3d=True), name="2-PhaseUnwrap")

    # Step 3 - Convert the fmap to Hz
    rads2hz = npe.Node(
        nutil.Function(
            input_names=["in_file", "delta_te"],
            output_names=["out_file"],
            function=rads2hz,
        ),
        name="3-ConvertFMapToHz",
    )

    # Step 4 - Call FUGUE to extrapolate from mask (fill holes, etc)
    pre_fugue = npe.Node(fsl.FUGUE(save_fmap=True), name="4-FugueExtrapolationFromMask")

    # Step 5 - Demean fmap to avoid gross shifting
    demean = npe.Node(
        nutil.Function(
            input_names=["in_file", "in_mask"],
            output_names=["out_file"],
            function=demean_image,
        ),
        name="5-DemeanFMap",
    )

    # Step 6 - Clean up edge voxels
    cleanup = cleanup_edge_pipeline(name="6-CleanUpEdgeVoxels")

    wf = npe.Workflow(name=name)
    # fmt: off
    wf.connect(
        [
            # Step 1 - Convert the fmap into radians
            (input_node, pha2rads, [("fmap_phasediff", "in_file")]),
            # Step 2 - Unwrap the fmap with PRELUDE
            (pha2rads, prelude, [("out_file", "phase_file")]),
            (input_node, prelude, [("fmap_magnitude", "magnitude_file")]),
            (input_node, prelude, [("fmap_mask", "mask_file")]),
            # Step 3 - Convert the fmap to Hz
            (prelude, rads2hz, [("unwrapped_phase_file", "in_file")]),
            (input_node, rads2hz, [("delta_echo_time", "delta_te")]),
            # Step 4 - Call FUGUE to extrapolate from mask (fill holes, etc)
            (rads2hz, pre_fugue, [("out_file", "fmap_in_file")]),
            (input_node, pre_fugue, [("fmap_mask", "mask_file")]),
            # Step 5 - Demean fmap to avoid gross shifting
            (pre_fugue, demean, [("fmap_out_file", "in_file")]),
            (input_node, demean, [("fmap_mask", "in_mask")]),
            # Step 6 - Clean up edge voxels
            (demean, cleanup, [("out_file", "inputnode.in_file")]),
            (input_node, cleanup, [("fmap_mask", "inputnode.in_mask")]),
            # Output node
            (cleanup, output_node, [("outputnode.out_file", "calibrated_fmap")]),
        ]
    )
    # fmt: on

    return wf


def compute_reference_b0(
    b_value_threshold: float,
    use_cuda: bool,
    initrand: bool,
    output_dir: Optional[str] = None,
    name: str = "compute_reference_b0",
) -> Workflow:
    """Step 1 of the DWI preprocessing using phase diff pipeline.

    Compute the reference b0 (i.e. average b0 with EPI distortions)

    This pipeline needs:
        - MRtrix3 to compute the whole brain mask
        - FSL to run Eddy and BET

    Parameters
    ----------
    b_value_threshold: float
        Threshold value to determine the B0 volumes in the DWI image

    use_cuda: bool
        Boolean to indicate whether cuda should be used or not

    initrand: bool
        ???

    output_dir: str, optional
        Path to output directory.
        If provided, the pipeline will write its output in this folder.
        Default to None.

    name: str, optional
        Name of the pipeline. Default='compute_reference_b0'.

    Returns
    -------
    Workflow :
        The Nipype workflow.
        This workflow has the following inputs:
            - "dwi_filename": The path to the DWI image
            - "b_vectors_filename": The path to the associated B-vectors file
            - "b_values_filename": The path to the associated B-values file
            - "total_readout_time": The total readout time extracted from JSON metadata
            - "phase_encoding_direction": The phase encoding direction extracted from JSON metadata
            - "image_id": Prefix to be used for output files
        And the following outputs:
            - "reference_b0": The path to the compute reference B0 volume
            - "brainmask": The path to the computed whole brain mask
    """
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.io as nio
    import nipype.interfaces.mrtrix3 as mrtrix3
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as npe

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        eddy_fsl_pipeline,
    )
    from clinica.utils.dwi import compute_average_b0

    from .dwi_preprocessing_using_phasediff_fmap_utils import get_grad_fsl

    inputnode = npe.Node(
        niu.IdentityInterface(
            fields=[
                "dwi_filename",
                "b_vectors_filename",
                "b_values_filename",
                "total_readout_time",
                "phase_encoding_direction",
                "image_id",
            ]
        ),
        name="inputnode",
    )

    fsl_gradient = npe.Node(
        niu.Function(
            input_names=["b_values_filename", "b_vectors_filename"],
            output_names=["grad_fsl"],
            function=get_grad_fsl,
        ),
        name="fsl_gradient",
    )

    # Compute whole brain mask
    brain_mask = npe.Node(mrtrix3.BrainMask(), name="brain_mask")
    brain_mask.inputs.out_file = "brainmask.nii.gz"

    # Run eddy without calibrated fmap
    pre_eddy = eddy_fsl_pipeline(
        use_cuda=use_cuda,
        initrand=initrand,
        image_id=True,
        name="pre_eddy",
    )

    # Compute the reference b0
    reference_b0 = npe.Node(
        niu.Function(
            input_names=["in_dwi", "in_bval"],
            output_names=["out_b0_average"],
            function=compute_average_b0,
        ),
        name="reference_b0",
    )
    reference_b0.inputs.low_bval = b_value_threshold

    # Compute brain mask from reference b0
    masked_reference_b0 = npe.Node(
        fsl.BET(mask=True, robust=True), name="masked_reference_b0"
    )

    outputnode = npe.Node(
        niu.IdentityInterface(fields=["reference_b0", "brainmask"]),
        name="outputnode",
    )

    if output_dir:
        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = output_dir
        write_results.inputs.parameterization = False

    wf = npe.Workflow(name=name)

    connections = [
        (
            inputnode,
            fsl_gradient,
            [
                ("b_values_filename", "b_values_filename"),
                ("b_vectors_filename", "b_vectors_filename"),
            ],
        ),
        (fsl_gradient, brain_mask, [("grad_fsl", "grad_fsl")]),
        (inputnode, brain_mask, [("dwi_filename", "in_file")]),
        (
            inputnode,
            pre_eddy,
            [
                ("total_readout_time", "inputnode.total_readout_time"),
                ("phase_encoding_direction", "inputnode.phase_encoding_direction"),
                ("dwi_filename", "inputnode.dwi_filename"),
                ("b_values_filename", "inputnode.b_values_filename"),
                ("b_vectors_filename", "inputnode.b_vectors_filename"),
                ("image_id", "inputnode.image_id"),
            ],
        ),
        (brain_mask, pre_eddy, [("out_file", "inputnode.in_mask")]),
        (inputnode, reference_b0, [("b_values_filename", "in_bval")]),
        (pre_eddy, reference_b0, [("outputnode.out_corrected", "in_dwi")]),
        (reference_b0, masked_reference_b0, [("out_b0_average", "in_file")]),
        (masked_reference_b0, outputnode, [("out_file", "reference_b0")]),
        (brain_mask, outputnode, [("out_file", "brainmask")]),
    ]
    if output_dir:
        connections += [
            (outputnode, write_results, [("reference_b0", "reference_b0")]),
            (outputnode, write_results, [("brainmask", "brainmask")]),
        ]

    wf.connect(connections)

    return wf
