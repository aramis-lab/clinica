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
