# coding: utf8


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
    from nipype.workflows.dmri.fsl.utils import (
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


def ants_bias_correction(low_bval=5.0, name="ants_bias_correction"):
    """This workflow reproduces dwibiascorrect script from MRtrix with ANTs algorithm."""
    import nipype.interfaces.ants as ants
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as npe
    from clinica.utils.dwi import compute_average_b0

    input_node = npe.Node(
        niu.IdentityInterface(fields=["dwi", "bval", "bvec", "mask"]), name="input_node"
    )

    output_node = npe.Node(
        niu.IdentityInterface(fields=["bias_corrected_dwi"]), name="output_node"
    )

    # Compute average b0
    compute_average_b0 = npe.Node(
        niu.Function(
            input_names=["in_dwi", "in_bval"],
            output_names=["out_b0_average"],
            function=compute_average_b0,
        ),
        name="ComputeB0Average",
    )
    compute_average_b0.inputs.low_bval = low_bval

    # Estimate bias field
    n4 = npe.Node(
        ants.N4BiasFieldCorrection(
            dimension=3,
            save_bias=True,
            shrink_factor=4,
            bspline_fitting_distance=100,
            bspline_order=3,
            n_iterations=[1000],
            convergence_threshold=0.0,
        ),
        name="BiasB0",
    )

    # Split DWI dataset
    split = npe.Node(fsl.Split(dimension="t"), name="SplitDWIs")

    # Remove bias to each DWI volume
    rm_bias = npe.MapNode(
        fsl.MultiImageMaths(op_string="-div %s"),
        iterfield=["in_file"],
        name="RemoveBiasOfDWIs",
    )

    # Remove negative values
    rm_negative = npe.MapNode(
        fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="RemoveNegative"
    )

    # Merge corrected DWIs
    merge = npe.Node(fsl.utils.Merge(dimension="t"), name="MergeDWIs")

    wf = npe.Workflow(name=name)
    # fmt: off
    wf.connect([
        # Compute average b0
        (input_node, compute_average_b0, [("dwi", "in_dwi"),
                                          ("bval", "in_bval")]),
        # Estimate bias field
        # Note from MRtrix developers:
        # Use the brain mask as a weights image rather than a mask; means that voxels at the edge of the mask
        # will have a smoothly-varying bias field correction applied, rather than multiplying by 1.0 outside the mask
        (input_node, n4, [("mask", "weight_image")]),
        (compute_average_b0, n4, [("out_b0_average", "input_image")]),
        # Split DWI dataset
        (input_node, split, [("dwi", "in_file")]),
        # Remove bias to each DWI volume
        (n4, rm_bias, [("bias_image", "operand_files")]),
        (split, rm_bias, [("out_files", "in_file")]),
        # Remove negative values
        (rm_bias, rm_negative, [("out_file", "in_file")]),
        # Merge corrected DWIs
        (rm_negative, merge, [("out_file", "in_files")]),
        # Output node
        (merge, output_node, [("merged_file", "bias_corrected_dwi")]),
    ])
    # fmt: on
    return wf
