from typing import Optional

from nipype.pipeline.engine import Workflow


def prepare_phasediff_fmap(
    output_dir: Optional[str] = None,
    name: Optional[str] = "prepare_phasediff_fmap",
) -> Workflow:
    """This workflow adapts the fsl_prepare_fieldmap script from FSL for the FSL eddy command.

    Please note that the step 3 converts the field map into Hz instead of rad/s
    in the initial script because FSL eddy --field expects the field map to be in Hz.

    This workflow needs FSL.

    The workflow performs the following steps:
        - Step 1 - Convert the field map into radians
        - Step 2 - Unwrap the field map with PRELUDE
        - Step 3 - Convert the field map to Hz
        - Step 4 - Call FUGUE to extrapolate from mask (fill holes, etc)
        - Step 5 - Demean the field map to avoid gross shifting
        - Step 6 - Clean up the edge voxels

    Parameters
    ----------
    output_dir: str, optional
         Path to output directory.
         If provided, the pipeline will write its output in this folder.
         Default to None.

    name : str, optional
        Name of the workflow. Default="prepare_phasediff_fmap".

    Returns
    -------
    Workflow :
        The Nipype workflow.
         This workflow has the following inputs:
            - "fmap_mask" : str, the path to the binary mask of the field map.
            - "fmap_phasediff" : str, the path to the phase difference field map.
            - "fmap_magnitude" : str, the path to the brain extracted magnitude field map.
              Chose the field map with the best contrast.
            - "delta_echo_time" : float, the DeltaEchoTime from BIDS specifications.
        And the following outputs:
            - "calibrated_fmap" : str, the path to the calibrated field map for eddy --field command.

    Warnings
    --------
    This workflow can not be used for PRELUDE.
    You would need to change the step 3 (conversion into rad/s instead of Hz).
    """
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.io as nio
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
            fields=[
                "fmap_mask",
                "fmap_phasediff",
                "fmap_magnitude",
                "delta_echo_time",
            ]
        ),
        name="input_node",
    )

    convert_fmap_to_rads = npe.Node(
        nutil.Function(
            input_names=["in_file"], output_names=["out_file"], function=siemens2rads
        ),
        name="convert_fmap_to_rads",
    )

    unwrap_fmap_with_prelude = npe.Node(
        fsl.PRELUDE(process3d=True), name="unwrap_fmap_with_prelude"
    )

    convert_fmap_to_hz = npe.Node(
        nutil.Function(
            input_names=["in_file", "delta_te"],
            output_names=["out_file"],
            function=rads2hz,
        ),
        name="convert_fmap_to_hz",
    )

    fugue_extrapolation_from_mask = npe.Node(
        fsl.FUGUE(save_fmap=True), name="fugue_extrapolation_from_mask"
    )

    demean_fmap = npe.Node(
        nutil.Function(
            input_names=["in_file", "in_mask"],
            output_names=["out_file"],
            function=demean_image,
        ),
        name="demean_fmap",
    )

    cleanup_edge_voxels = cleanup_edge_pipeline(name="cleanup_edge_voxels")

    output_node = npe.Node(
        nutil.IdentityInterface(fields=["calibrated_fmap"]),
        name="output_node",
    )

    if output_dir:
        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = output_dir
        write_results.inputs.parameterization = False

    wf = npe.Workflow(name=name)

    connections = [
        (input_node, convert_fmap_to_rads, [("fmap_phasediff", "in_file")]),
        (convert_fmap_to_rads, unwrap_fmap_with_prelude, [("out_file", "phase_file")]),
        (input_node, unwrap_fmap_with_prelude, [("fmap_magnitude", "magnitude_file")]),
        (input_node, unwrap_fmap_with_prelude, [("fmap_mask", "mask_file")]),
        (
            unwrap_fmap_with_prelude,
            convert_fmap_to_hz,
            [("unwrapped_phase_file", "in_file")],
        ),
        (input_node, convert_fmap_to_hz, [("delta_echo_time", "delta_te")]),
        (
            convert_fmap_to_hz,
            fugue_extrapolation_from_mask,
            [("out_file", "fmap_in_file")],
        ),
        (input_node, fugue_extrapolation_from_mask, [("fmap_mask", "mask_file")]),
        (fugue_extrapolation_from_mask, demean_fmap, [("fmap_out_file", "in_file")]),
        (input_node, demean_fmap, [("fmap_mask", "in_mask")]),
        (demean_fmap, cleanup_edge_voxels, [("out_file", "inputnode.in_file")]),
        (input_node, cleanup_edge_voxels, [("fmap_mask", "inputnode.in_mask")]),
        (
            cleanup_edge_voxels,
            output_node,
            [("outputnode.out_file", "calibrated_fmap")],
        ),
    ]

    if output_dir:
        connections += [
            (output_node, write_results, [("calibrated_fmap", "calibrated_fmap")]),
        ]
    wf.connect(connections)

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
        compute_mask=False,
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


def calibrate_and_register_fmap(
    output_dir: Optional[str] = None,
    name: str = "calibrate_and_register_fmap",
) -> Workflow:
    """Step 2 of the DWI preprocessing using phase diff pipeline.

    Calibrate and register the field map.

    This pipeline needs:
        - ANTS
        - FSL

    Parameters
    ----------
    output_dir: str, optional
        Path to output directory.
        If provided, the pipeline will write its output in this folder.
        Default to None.

    name: str, optional
        Name of the pipeline. Default='calibrate_and_register_fmap'.

    Returns
    -------
    Workflow :
        The Nipype workflow.
        This workflow has the following inputs:
            - "bias_magnitude_fmap": The path to the bias magnitude field map
            - "fmap_phasediff": ??
            - "reference_b0": The path to the computed reference B0 volume.
              This is one of the outputs from the step 1 compute_reference_b0 workflow.
            - "delta_echo_time": ???
        And the following outputs:
            - "smooth_calibrated_fmap": The path to the smoothed and calibrated field map
            - "bet_magnitude_fmap_registered_onto_b0": The path to the bet magnitude field
              map registered onto the B0 volume
    """
    import nipype.interfaces.ants as ants
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as npe

    inputnode = npe.Node(
        niu.IdentityInterface(
            fields=[
                "bias_magnitude_fmap",
                "fmap_phasediff",
                "reference_b0",
                "delta_echo_time",
            ]
        ),
        name="inputnode",
    )

    # Bias field correction of the magnitude image
    bias_magnitude_fmap = npe.Node(
        ants.N4BiasFieldCorrection(dimension=3),
        name="bias_magnitude_fmap",
    )

    # Brain extraction of the magnitude image
    bet_magnitude_fmap = npe.Node(
        fsl.BET(frac=0.4, mask=True),
        name="bet_magnitude_fmap",
    )

    # Calibrate FMap
    calibrate_fmap = prepare_phasediff_fmap(name="calibrate_fmap")

    # Register the BET magnitude fmap onto the BET b0
    register_bet_magnitude_fmap_onto_b0 = npe.Node(
        interface=fsl.FLIRT(),
        name="register_bet_magnitude_fmap_onto_b0",
    )
    register_bet_magnitude_fmap_onto_b0.inputs.dof = 6
    register_bet_magnitude_fmap_onto_b0.inputs.output_type = "NIFTI_GZ"

    # Apply the transformation on the calibrated fmap
    apply_xfm_on_calibrated_fmap = npe.Node(
        interface=fsl.ApplyXFM(),
        name="apply_xfm_on_calibrated_fmap",
    )
    apply_xfm_on_calibrated_fmap.inputs.output_type = "NIFTI_GZ"

    # Apply the transformation on the magnitude image
    apply_xfm_on_magnitude_fmap = apply_xfm_on_calibrated_fmap.clone("2e-2-MagFMapToB0")

    # Smooth the registered (calibrated) fmap
    smooth_calibrated_fmap = npe.Node(
        interface=fsl.maths.IsotropicSmooth(),
        name="smooth_calibrated_fmap",
    )
    smooth_calibrated_fmap.inputs.sigma = 4.0

    outputnode = npe.Node(
        niu.IdentityInterface(
            fields=[
                "smooth_calibrated_fmap",
                "bet_magnitude_fmap_registered_onto_b0",
                "registered_calibrated_fmap",
            ]
        ),
        name="outputnode",
    )

    if output_dir:
        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = output_dir
        write_results.inputs.parameterization = False

    wf = npe.Workflow(name=name)

    connections = [
        # Bias field correction of the magnitude image
        (inputnode, bias_magnitude_fmap, [("bias_magnitude_fmap", "input_image")]),
        # Brain extraction of the magnitude image
        (bias_magnitude_fmap, bet_magnitude_fmap, [("output_image", "in_file")]),
        # Calibration of the FMap
        (
            bet_magnitude_fmap,
            calibrate_fmap,
            [
                ("mask_file", "input_node.fmap_mask"),
                ("out_file", "input_node.fmap_magnitude"),
            ],
        ),
        (
            inputnode,
            calibrate_fmap,
            [
                ("fmap_phasediff", "input_node.fmap_phasediff"),
                ("delta_echo_time", "input_node.delta_echo_time"),
            ],
        ),
        # Register the BET magnitude fmap onto the BET b0
        (
            bet_magnitude_fmap,
            register_bet_magnitude_fmap_onto_b0,
            [("out_file", "in_file")],
        ),
        (
            inputnode,
            register_bet_magnitude_fmap_onto_b0,
            [("reference_b0", "reference")],
        ),
        # Apply the transformation on the magnitude image
        (
            register_bet_magnitude_fmap_onto_b0,
            apply_xfm_on_magnitude_fmap,
            [("out_matrix_file", "in_matrix_file")],
        ),
        (
            bias_magnitude_fmap,
            apply_xfm_on_magnitude_fmap,
            [("output_image", "in_file")],
        ),
        (inputnode, apply_xfm_on_magnitude_fmap, [("reference_b0", "reference")]),
        # Apply the transformation on the calibrated fmap
        (
            register_bet_magnitude_fmap_onto_b0,
            apply_xfm_on_calibrated_fmap,
            [("out_matrix_file", "in_matrix_file")],
        ),
        (
            calibrate_fmap,
            apply_xfm_on_calibrated_fmap,
            [("output_node.calibrated_fmap", "in_file")],
        ),
        (inputnode, apply_xfm_on_calibrated_fmap, [("reference_b0", "reference")]),
        # Smooth the registered (calibrated) fmap
        (
            apply_xfm_on_calibrated_fmap,
            smooth_calibrated_fmap,
            [("out_file", "in_file")],
        ),
        # Output connections
        (smooth_calibrated_fmap, outputnode, [("out_file", "smooth_calibrated_fmap")]),
        (
            register_bet_magnitude_fmap_onto_b0,
            outputnode,
            [("out_file", "bet_magnitude_fmap_registered_onto_b0")],
        ),
        (
            apply_xfm_on_calibrated_fmap,
            outputnode,
            [("out_file", "registered_calibrated_fmap")],
        ),
    ]

    if output_dir:
        connections += [
            (
                outputnode,
                write_results,
                [("smooth_calibrated_fmap", "smooth_calibrated_fmap")],
            ),
            (
                outputnode,
                write_results,
                [
                    (
                        "bet_magnitude_fmap_registered_onto_b0",
                        "bet_magnitude_fmap_registered_onto_b0",
                    )
                ],
            ),
            (
                outputnode,
                write_results,
                [("registered_calibrated_fmap", "registered_calibrated_fmap")],
            ),
        ]

    wf.connect(connections)

    return wf
