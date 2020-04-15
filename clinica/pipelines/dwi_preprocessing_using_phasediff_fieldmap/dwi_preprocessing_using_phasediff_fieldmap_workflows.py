# coding: utf8


def susceptibility_distortion_correction_using_phasediff_fmap(
        fugue_params=dict(smooth3d=2.0),
        register_fmap_on_b0=True,
        name='susceptibility_distortion_correction_using_phasediff_fmap'):
    """
    The fieldmap based method implements susceptibility distortion correction
    by using a mapping of the B0 field as proposed by [Jezzard95]_.
    This workflow uses the implementation of FSL
    (`FUGUE <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE>`_). Phase
    unwrapping is performed using
    `PRELUDE <http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/fugue/prelude.html>`_
    [Jenkinson03]_. Preparation of the fieldmap is performed reproducing the
    script in FSL `fsl_prepare_fieldmap
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data>`_.

    Args:
        fugue_params: Regularisation to apply to the fieldmap with a
            3D Gaussian smoothing (default: 2.0)
        register_fmap_on_b0: Register the fmap on the b0 or not. If the fmap
            is already well aligned and
        name:

    Inputnode
        in_file (str): DWI dataset.
        in_bval (str): Bval file.
        in_mask (str): Mask file.
        in_fmap_magnitude (str): Magnitude fieldmap. Chose the fmap with the
            best contrast.
        in_fmap_phasediff (str): Phase difference fieldmap.
        delta_echo_time (float): Difference of echo time of the
        effective_echo_spacing (float):
        phase_encoding_direction (str):

    Outputnode:
        out_file (str): Output.
        out_vsm (str): The set of DWI volumes.
        out_warp (str): The bvalues corresponding to the out_dwi.


    .. warning:: Only SIEMENS format fieldmaps are supported.

    .. admonition:: References
<
      .. [Jezzard95] Jezzard P, and Balaban RS, `Correction for geometric
        distortion in echo planar images from B0 field variations
        <http://dx.doi.org/10.1002/mrm.1910340111>`_,
        MRM 34(1):65-73. (1995). doi: 10.1002/mrm.1910340111.

      .. [Jenkinson03] Jenkinson M., `Fast, automated, N-dimensional
        phase-unwrapping algorithm <http://dx.doi.org/10.1002/mrm.10354>`_,
        MRM 49(1):193-197, 2003, doi: 10.1002/mrm.10354.
    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants

    from clinica.utils.epi import bids_dir_to_fsl_dir
    from clinica.utils.fmap import resample_fmap_to_b0

    from nipype.workflows.dmri.fsl.utils import (rads2radsec, siemens2rads,
                                                 vsm2warp, add_empty_vol,
                                                 cleanup_edge_pipeline,
                                                 demean_image)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_dwi', 'in_mask', 'in_fmap_phasediff', 'in_fmap_magnitude',
                'delta_echo_time', 'effective_echo_spacing',
                'phase_encoding_direction']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_vsm', 'out_warp', 'out_prepared_fmap']),
        name='outputnode')

    fsl_dir = pe.Node(niu.Function(
        input_names=['bids_dir'],
        output_names=['fsl_dir'],
        function=bids_dir_to_fsl_dir), name='ConvertToFslDir')

    get_b0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='GetB0')

    first_mag = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='GetFirst')

    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='N4MagnitudeFmap')

    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='BetN4MagnitudeFmap')

    dilate = pe.Node(
            fsl.maths.MathsCommand(
                nan2zeros=True,
                args='-kernel sphere 5 -dilM'),
            name='DilateBet')

    pha2rads = pe.Node(niu.Function(
        input_names=['in_file'],
        output_names=['out_file'],
        function=siemens2rads), name='PreparePhase')

    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')

    rad2rsec = pe.Node(niu.Function(
        input_names=['in_file', 'delta_te'],
        output_names=['out_file'],
        function=rads2radsec), name='ToRadSec')

    # Node when we register the fmap onto the b0
    fmm2b0 = pe.Node(ants.Registration(
        output_warped_image=True), name="FmapMagnitudeToB0")
    fmm2b0.inputs.transforms = ['Rigid'] * 2
    fmm2b0.inputs.transform_parameters = [(1.0,)] * 2
    fmm2b0.inputs.number_of_iterations = [[50], [20]]
    fmm2b0.inputs.dimension = 3
    fmm2b0.inputs.metric = ['Mattes', 'Mattes']
    fmm2b0.inputs.metric_weight = [1.0] * 2
    fmm2b0.inputs.radius_or_number_of_bins = [64, 64]
    fmm2b0.inputs.sampling_strategy = ['Regular', 'Random']
    fmm2b0.inputs.sampling_percentage = [None, 0.2]
    fmm2b0.inputs.convergence_threshold = [1.e-5, 1.e-8]
    fmm2b0.inputs.convergence_window_size = [20, 10]
    fmm2b0.inputs.smoothing_sigmas = [[6.0], [2.0]]
    fmm2b0.inputs.sigma_units = ['vox'] * 2
    fmm2b0.inputs.shrink_factors = [[6], [1]]  # ,[1] ]
    fmm2b0.inputs.use_estimate_learning_rate_once = [True] * 2
    fmm2b0.inputs.use_histogram_matching = [True] * 2
    fmm2b0.inputs.initial_moving_transform_com = 0
    fmm2b0.inputs.collapse_output_transforms = True
    fmm2b0.inputs.winsorize_upper_quantile = 0.995

    # Node when we resample the fmap onto the b0
    res_fmap = pe.Node(niu.Function(
        input_names=['in_fmap', 'in_b0', 'out_file'],
        output_names=['out_resampled_fmap'],
        function=resample_fmap_to_b0), name='ResampleFmap')

    apply_xfm = pe.Node(ants.ApplyTransforms(
        dimension=3, interpolation='BSpline'), name='FmapPhaseToB0')

    pre_fugue = pe.Node(fsl.FUGUE(save_fmap=True), name='PreliminaryFugue')

    demean = pe.Node(
            niu.Function(
                input_names=['in_file', 'in_mask'],
                output_names=['out_file'],
                function=demean_image),
            name='DemeanFmap')

    cleanup = cleanup_edge_pipeline()

    add_vol = pe.Node(niu.Function(
        input_names=['in_file'],
        output_names=['out_file'],
        function=add_empty_vol), name='AddEmptyVol')

    vsm = pe.Node(fsl.FUGUE(save_shift=True, **fugue_params),
                  name="ComputeVSM")

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    unwarp = pe.MapNode(fsl.FUGUE(icorr=True, forward_warping=False),
                        iterfield=['in_file'], name='UnwarpDWIs')

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    vsm2dfm = vsm2warp()
    vsm2dfm.inputs.inputnode.scaling = 1.0

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, fsl_dir, [('phase_encoding_direction', 'bids_dir')]),  # noqa
        (inputnode, pha2rads, [('in_fmap_phasediff', 'in_file')]),  # noqa
        (inputnode, get_b0,   [('in_dwi', 'in_file')]),  # noqa
        (inputnode, first_mag, [('in_fmap_magnitude', 'in_file')]),  # noqa
        (first_mag, n4, [('roi_file', 'input_image')]),  # noqa
        (n4, bet, [('output_image', 'in_file')]),  # noqa
        (bet, dilate, [('mask_file', 'in_file')]),  # noqa
        (pha2rads, prelude, [('out_file',     'phase_file')]),  # noqa
        (n4,       prelude, [('output_image', 'magnitude_file')]),  # noqa
        (dilate,   prelude, [('out_file',     'mask_file')]),  # noqa
        (prelude,   rad2rsec, [('unwrapped_phase_file', 'in_file')]),  # noqa
        (inputnode, rad2rsec, [('delta_echo_time',      'delta_te')]),  # noqa
    ])
    if register_fmap_on_b0:
        wf.connect([
            (get_b0,    fmm2b0, [('roi_file',     'fixed_image')]),  # noqa
            (n4,        fmm2b0, [('output_image', 'moving_image')]),  # noqa
            (inputnode, fmm2b0, [('in_mask',      'fixed_image_mask')]),  # noqa
            (dilate,    fmm2b0, [('out_file',     'moving_image_mask')]),  # noqa
            (get_b0,   apply_xfm, [('roi_file',            'reference_image')]),  # noqa
            (rad2rsec, apply_xfm, [('out_file',            'input_image')]),  # noqa
            (fmm2b0,   apply_xfm, [('forward_transforms',  'transforms'),  # noqa
                                   ('forward_invert_flags', 'invert_transform_flags')]),  # noqa
            (apply_xfm, pre_fugue, [('output_image', 'fmap_in_file')]),  # noqa
            (inputnode, pre_fugue, [('in_mask', 'mask_file')]),  # noqa

            (apply_xfm, outputnode, [('output_image', 'out_prepared_fmap')])  # noqa
        ])
    else:
        wf.connect([
            (get_b0,   res_fmap, [('roi_file', 'in_b0')]),  # noqa
            (rad2rsec, res_fmap, [('out_file', 'in_fmap')]),  # noqa
            (res_fmap,  pre_fugue, [('out_resampled_fmap', 'fmap_in_file')]),  # noqa
            (inputnode, pre_fugue, [('in_mask',            'mask_file')]),  # noqa

            (res_fmap, outputnode, [('out_resampled_fmap', 'out_prepared_fmap')])  # noqa
        ])
    wf.connect([
        (pre_fugue, demean, [('fmap_out_file', 'in_file')]),  # noqa
        (inputnode, demean, [('in_mask', 'in_mask')]),  # noqa
        (demean,    cleanup, [('out_file', 'inputnode.in_file')]),  # noqa
        (inputnode, cleanup, [('in_mask',  'inputnode.in_mask')]),  # noqa
        (cleanup, add_vol, [('outputnode.out_file', 'in_file')]),  # noqa
        (inputnode, vsm, [('in_mask',                'mask_file')]),  # noqa
        (add_vol,   vsm, [('out_file',               'fmap_in_file')]),  # noqa
        (inputnode, vsm, [('delta_echo_time',        'asym_se_time')]),  # noqa
        (inputnode, vsm, [('effective_echo_spacing', 'dwell_time')]),  # noqa
        (inputnode, split, [('in_dwi', 'in_file')]),  # noqa
        (split,   unwarp, [('out_files',      'in_file')]),  # noqa
        (vsm,     unwarp, [('shift_out_file', 'shift_in_file')]),  # noqa
        (fsl_dir, unwarp, [('fsl_dir',        'unwarp_direction')]),  # noqa
        (unwarp, thres, [('unwarped_file', 'in_file')]),  # noqa
        (thres, merge, [('out_file', 'in_files')]),  # noqa
        (merge,   vsm2dfm, [('merged_file',    'inputnode.in_ref')]),  # noqa
        (vsm,     vsm2dfm, [('shift_out_file', 'inputnode.in_vsm')]),  # noqa
        (fsl_dir, vsm2dfm, [('fsl_dir',        'inputnode.enc_dir')]),  # noqa
        (rad2rsec, outputnode, [('out_file',            'out_native_fmap')]),  # noqa
        (merge,    outputnode, [('merged_file',         'out_file')]),  # noqa
        (vsm,      outputnode, [('shift_out_file',      'out_vsm')]),  # noqa
        (vsm2dfm,  outputnode, [('outputnode.out_warp', 'out_warp')]),  # noqa
    ])

    return wf
