# coding: utf8


def prepare_phasediff_fmap(name='prepare_phasediff_fmap'):
    """
    This workflow adapts the fsl_prepare_fieldmap script from FSL for the FSL eddy command.

    Please note that the step 3 converts the fieldmap into Hz instead of rad/s
    in the initial script because FSL eddy --field expects the fieldmap to be in Hz.

    Input node:
        fmap_mask (str): Binary mask of the fieldmap.
        fmap_phasediff (str): Phase difference fieldmap.
        fmap_magnitude (str): Magnitude fieldmap. Chose the fieldmap with the best contrast.
        delta_echo_time (float): DeltaEchoTime from BIDS specifications.

    Output node:
        calibrated_fmap (str): Calibrated fieldmap for eddy --field command.

    Warnings:
        This workflow can not be used for PRELUDE. You would need
        to change the step 3 (conversion into rad/s instead of Hz)
    """
    import nipype.interfaces.utility as nutil
    import nipype.pipeline.engine as npe
    import nipype.interfaces.fsl as fsl
    from nipype.workflows.dmri.fsl.utils import (siemens2rads, cleanup_edge_pipeline, demean_image)
    from clinica.utils.fmap import rads2hz

    input_node = npe.Node(nutil.IdentityInterface(
        fields=['fmap_mask', 'fmap_phasediff', 'fmap_magnitude', 'delta_echo_time']),
        name='input_node')

    output_node = npe.Node(nutil.IdentityInterface(
        fields=['calibrated_fmap']),
        name='output_node')

    # Step 1 - Convert the fmap into radians
    pha2rads = npe.Node(nutil.Function(input_names=['in_file'],
                                       output_names=['out_file'],
                                       function=siemens2rads),
                        name='1-ConvertFMapToRads')

    # Step 2a - Dilate brain mask for PRELUDE
    dilate = npe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                                             # args='-kernel sphere 5 -dilM'),
                                             args='-ero'),
                      name='2a-DilateBrainMask')

    # Step 2b - Unwrap the fmap with PRELUDE
    prelude = npe.Node(fsl.PRELUDE(process3d=True),
                       name='2-PhaseUnwrap')

    # Step 3 - Convert the fmap to Hz
    rads2hz = npe.Node(nutil.Function(input_names=['in_file', 'delta_te'],
                                      output_names=['out_file'],
                                      function=rads2hz),
                       name='3-ConvertFMapToHz')

    # Step 4 - Call FUGUE to extrapolate from mask (fill holes, etc)
    pre_fugue = npe.Node(fsl.FUGUE(save_fmap=True),
                         name='4-FugueExtrapolationFromMask')

    # Step 5 - Demean fmap to avoid gross shifting
    demean = npe.Node(nutil.Function(input_names=['in_file', 'in_mask'],
                                     output_names=['out_file'],
                                     function=demean_image),
                      name='5-DemeanFMap')

    # Step 6 - Clean up edge voxels
    cleanup = cleanup_edge_pipeline(
        name='6-CleanUpEdgeVoxels')

    wf = npe.Workflow(name=name)
    wf.connect([
        # Step 1 - Convert the fmap into radians
        (input_node, pha2rads, [('fmap_phasediff', 'in_file')]),  # noqa
        # Step 2a - Dilate brain mask for prelude
        (input_node, dilate, [('fmap_mask', 'in_file')]),  # noqa
        # Step 2b - Unwrap the fmap with PRELUDE
        (pha2rads,   prelude, [('out_file', 'phase_file')]),  # noqa
        (input_node, prelude, [('fmap_magnitude', 'magnitude_file')]),  # noqa
        (dilate,     prelude, [('out_file', 'mask_file')]),  # noqa
        # Step 3 - Convert the fmap to Hz
        (prelude,    rads2hz, [('unwrapped_phase_file', 'in_file')]),  # noqa
        (input_node, rads2hz, [('delta_echo_time',      'delta_te')]),  # noqa
        # Step 4 - Call FUGUE to extrapolate from mask (fill holes, etc)
        (rads2hz,   pre_fugue, [('out_file', 'fmap_in_file')]),  # noqa
        (dilate,     pre_fugue, [('out_file', 'mask_file')]),  # noqa
#        (input_node, pre_fugue, [('fmap_mask',  'mask_file')]),  # noqa
        # Step 5 - Demean fmap to avoid gross shifting
        (pre_fugue, demean, [('fmap_out_file', 'in_file')]),  # noqa
        (dilate, demean, [('out_file', 'in_mask')]),  # noqa
#        (input_node, demean, [('fmap_mask', 'in_mask')]),  # noqa
        # Step 6 - Clean up edge voxels
        (demean,     cleanup, [('out_file', 'inputnode.in_file')]),  # noqa
        (dilate, cleanup, [('out_file', 'inputnode.in_mask')]),  # noqa
#        (input_node, cleanup, [('fmap_mask',  'inputnode.in_mask')]),  # noqa
        # Output node
        (cleanup, output_node, [('outputnode.out_file', 'calibrated_fmap')]),  # noqa
    ])

    return wf


def remove_bias(name='bias_correct'):
    """
    This workflow estimates a single multiplicative bias field from the
    averaged *b0* image, as suggested in [Jeurissen2014]_.
    .. admonition:: References
      .. [Jeurissen2014] Jeurissen B. et al., `Multi-tissue constrained
        spherical deconvolution for improved analysis of multi-shell diffusion
        MRI data <http://dx.doi.org/10.1016/j.neuroimage.2014.07.061>`_.squeue

        NeuroImage (2014). doi: 10.1016/j.neuroimage.2014.07.061
    """
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'b0_mask']),
                         name='outputnode')

    get_b0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')

    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                      name='mask_b0')

    n4 = pe.Node(ants.N4BiasFieldCorrection(
        dimension=3, save_bias=True, bspline_fitting_distance=600),
        name='Bias_b0')
    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    mult = pe.MapNode(fsl.MultiImageMaths(op_string='-div %s'),
                      iterfield=['in_file'], name='RemoveBiasOfDWIs')
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')
    merge = pe.Node(fsl.utils.Merge(dimension='t'), name='MergeDWIs')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, get_b0, [('in_file', 'in_file')]),
        (get_b0,   n4, [('roi_file', 'input_image')]),
        (get_b0, mask_b0, [('roi_file', 'in_file')]),
        (mask_b0, n4, [('mask_file', 'mask_image')]),
        (inputnode, split, [('in_file', 'in_file')]),
        (n4,    mult, [('bias_image', 'operand_files')]),
        (split, mult, [('out_files', 'in_file')]),
        (mult, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),
        (merge,   outputnode, [('merged_file', 'out_file')]),
        (mask_b0, outputnode, [('mask_file', 'b0_mask')])
    ])
    return wf
