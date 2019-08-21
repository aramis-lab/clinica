# coding: utf8

"""This module contains workflows for DWI preprocessing."""


def prepare_b0(num_b0s, low_bval=5, name='prepare_b0'):
    """
    Create a pipelines that prepare the data for further corrections. This
    pipelines coregister the B0 images and then
    average it in order to obtain only one average B0 images.
    The b-vectors and b-values are updated according to the modifications.

    Args:
        num_b0s (int): Number of b0 in the DWI dataset
        low_bval:
        name (Optional[str]): Name of the workflow

    Inputnode:
        in_dwi (str): Input DWI file.
        in_bvec (str): Vector file of the diffusion directions
            of the dwi dataset.
        in_bval (str): B-values file.

    Outputnode:
        dwi_b0_merge: Average of B0 images merged to the DWIs
        b0_reference: Average of the B0 images or the only B0 image
        out_bvec: Updated gradient vectors table
        out_bvals: Updated gradient values table
        mask_b0: Binary mask obtained from the average of the B0 images

    Returns:
        The workflow
    """
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.workflows.dmri.fsl.utils import b0_average

    from clinica.utils.dwi import b0_dwi_split, insert_b0_into_dwi

    inputnode = pe.Node(interface=niu.IdentityInterface(
        fields=["in_dwi", "in_bvec", "in_bval"]),
        name="inputnode")

    b0_dwi_split = pe.Node(niu.Function(
        input_names=['in_file', 'in_bvals', 'in_bvecs', 'low_bval'],
        output_names=['out_b0', 'out_dwi', 'out_bvals', 'out_bvecs'],
        function=b0_dwi_split), name='b0_dwi_split')
    b0_dwi_split.inputs.lowbval = low_bval

    b0_flirt = b0_flirt_pipeline(num_b0s=num_b0s, name='b0_co_registration')

    b0_avg = pe.Node(niu.Function(
        input_names=['in_file'], output_names=['out_file'],
        function=b0_average), name='b0_average')

    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                      name='mask_b0')

    insert_b0_into_dwi = pe.Node(niu.Function(
        input_names=['in_b0', 'in_dwi', 'in_bvals', 'in_bvecs'],
        output_names=['out_dwi', 'out_bvals', 'out_bvecs'],
        function=insert_b0_into_dwi), name='insert_b0avg_into_dwi')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['mask_b0', 'b0_reference', 'dwi_b0_merge',
                'out_bvecs', 'out_bvals']),
        name='outputnode')

    wf = pe.Workflow(name=name)

    if num_b0s == 1:
        wf.connect([
            # Split dataset into two datasets (b0, b>low_bval)
            (inputnode,          b0_dwi_split, [('in_bval',       'in_bvals'),
                                                ('in_bvec',       'in_bvecs'),
                                                ('in_dwi',       'in_file')]),
            # Merge datasets such that bval(DWI) = (0 b1 ... bn)
            (b0_dwi_split, insert_b0_into_dwi, [('out_b0',           'in_b0'),
                                                ('out_dwi',         'in_dwi'),
                                                ('out_bvals',     'in_bvals'),
                                                ('out_bvecs',   'in_bvecs')]),
            # Compute b0 mask
            (b0_dwi_split,            mask_b0, [('out_b0',       'in_file')]),
            # Outputnode
            (insert_b0_into_dwi,   outputnode, [('out_dwi',   'dwi_b0_merge'),
                                                ('out_bvals',    'out_bvals'),
                                                ('out_bvecs',  'out_bvecs')]),
            (mask_b0,              outputnode, [('mask_file',    'mask_b0')]),
            (b0_dwi_split,         outputnode, [('out_b0',  'b0_reference')])
        ])
    elif num_b0s > 1:
        wf.connect([
            # Split dataset into two datasets (b0s, b>low_bval)
            (inputnode,             b0_dwi_split, [('in_bval',              'in_bvals'),  # noqa
                                                   ('in_bvec',              'in_bvecs'),  # noqa
                                                   ('in_dwi',              'in_file')]),  # noqa
            # Register the b0 onto the first b0
            (b0_dwi_split,              b0_flirt, [('out_b0',    'inputnode.in_file')]),  # noqa
            # Average the b0s
            (b0_flirt,                    b0_avg, [('outputnode.out_file', 'in_file')]),  # noqa
            # Compute b0 mask from b0avg
            (b0_avg,                     mask_b0, [('out_file',            'in_file')]),  # noqa
            # Merge datasets such that bval(DWI) = (0 b1 ... bn)
            (b0_avg,          insert_b0_into_dwi, [('out_file',              'in_b0')]),  # noqa
            (b0_dwi_split,    insert_b0_into_dwi, [('out_dwi',                'in_dwi'),  # noqa
                                                   ('out_bvals',            'in_bvals'),  # noqa
                                                   ('out_bvecs',          'in_bvecs')]),  # noqa
            # Outputnode
            (insert_b0_into_dwi,      outputnode, [('out_dwi',          'dwi_b0_merge'),  # noqa
                                                   ('out_bvals',           'out_bvals'),  # noqa
                                                   ('out_bvecs',         'out_bvecs')]),  # noqa
            (mask_b0,                 outputnode, [('mask_file',           'mask_b0')]),  # noqa
            (b0_avg,                  outputnode, [('out_file',       'b0_reference')])   # noqa
        ])
    else:
        raise ValueError('The number of b0s should be strictly positive.')

    return wf


def b0_flirt_pipeline(num_b0s, name='b0_coregistration'):
    """
    Rigid registration of the B0 dataset onto the first volume. Rigid
    registration is achieved using FLIRT and the normalized
    correlation.

    Args:
        num_b0s (int): Number of the B0 volumes in the dataset.
        name (str): Name of the workflow.

    Inputnode:
        in_file(str): B0 dataset.

    Outputnode
        out_b0_reg(str): The set of B0 volumes registered to the first volume.

    Returns:
        The workflow
    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces import fsl
    import nipype.interfaces.utility as niu

    from clinica.utils.dwi import merge_volumes_tdim

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']),
                        name='inputnode')
    fslroi_ref = pe.Node(fsl.ExtractROI(args='0 1'), name='b0_reference')
    tsize = num_b0s - 1
    fslroi_moving = pe.Node(fsl.ExtractROI(args='1 '+str(tsize)),
                            name='b0_moving')
    split_moving = pe.Node(fsl.Split(dimension='t'), name='split_b0_moving')

    bet_ref = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                      name='bet_ref')

    dilate = pe.Node(
            fsl.maths.MathsCommand(
                nan2zeros=True,
                args='-kernel sphere 5 -dilM'),
            name='mask_dilate')

    flirt = pe.MapNode(fsl.FLIRT(
        interp='spline', dof=6, bins=50, save_log=True,
        cost='corratio', cost_func='corratio', padding_size=10,
        searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4],
        fine_search=1, coarse_search=10),
        name='b0_co_registration', iterfield=['in_file'])

    merge = pe.Node(fsl.Merge(dimension='t'), name='merge_registered_b0s')
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='remove_negative')
    insert_ref = pe.Node(niu.Function(input_names=['in_file1', 'in_file2'],
                                      output_names=['out_file'],
                                      function=merge_volumes_tdim),
                         name='concat_ref_moving')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_xfms']),
        name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  fslroi_ref,   [('in_file', 'in_file')]),
        (inputnode,  fslroi_moving,   [('in_file', 'in_file')]),
        (fslroi_moving, split_moving,   [('roi_file', 'in_file')]),
        (fslroi_ref, bet_ref, [('roi_file', 'in_file')]),
        (bet_ref, dilate, [('mask_file', 'in_file')]),
        (dilate, flirt, [('out_file', 'ref_weight'),
                         ('out_file', 'in_weight')]),
        (fslroi_ref, flirt, [('roi_file', 'reference')]),
        (split_moving, flirt, [('out_files', 'in_file')]),
        (flirt, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),
        (merge, insert_ref, [('merged_file', 'in_file2')]),
        (fslroi_ref, insert_ref, [('roi_file', 'in_file1')]),
        (insert_ref, outputnode, [('out_file', 'out_file')]),
        (flirt, outputnode, [('out_matrix_file', 'out_xfms')])
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
    Example
    -------
    >>> from nipype.workflows.dmri.fsl.artifacts import remove_bias
    >>> bias = remove_bias()
    >>> bias.inputs.inputnode.in_file = 'epi.nii'
    >>> bias.inputs.inputnode.in_bval = 'diffusion.bval'
    >>> bias.inputs.inputnode.in_mask = 'mask.nii'
    >>> bias.run() # doctest: +SKIP
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


def eddy_fsl_pipeline(epi_param, name='eddy_fsl'):
    """
    Using eddy from fsl for head motion correction and eddy current distortion correction.

    """
    from nipype.interfaces.fsl import Eddy
    import nipype.interfaces.utility as niu     # utility
    import nipype.pipeline.engine as pe          # pypeline engine
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import eddy_fsl, generate_acq, generate_index, b0_indices

    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['in_file',
                        'in_bvec',
                        'in_bval',
                        'in_mask',
                        'ref_b0']),
            name='inputnode')

    generate_acq = pe.Node(niu.Function(function=generate_acq,
                                        input_names=['in_b0', 'epi_param'],
                                        output_names=['out_file']),
                           name='generate_acq')
    generate_acq.inputs.epi_param = epi_param

    list_b0 = pe.Node(niu.Function(input_names=['in_bval'],
                                   output_names=['out_idx'],
                                   function=b0_indices),
                      name='find_b0_indices')

    generate_index = pe.Node(niu.Function(function=generate_index,
                                          input_names=['in_bval', 'b0_index'],
                                          output_names=['eddy_index']),
                             name='generate_index')

    eddy = pe.Node(niu.Function(input_names=['in_bvec', 'in_bval', 'in_file', 'in_mask', 'in_acqp', 'in_index'],
                                output_names=['out_parameter', 'out_corrected', 'out_rotated_bvecs'],
                                function=eddy_fsl),
                   name='eddy_fsl')

    eddy = pe.Node(interface=Eddy(), name='eddy_fsl')
    eddy.inputs.flm = 'linear'

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_parameter',
                                                       'out_corrected',
                                                       'out_rotated_bvecs']),
                         name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,     generate_acq,   [('ref_b0', 'in_b0')]),
        (inputnode,  generate_index,      [('in_bval', 'in_bval')]),
        (list_b0,  generate_index,      [('out_idx', 'b0_index')]),
        (inputnode,      list_b0,      [('in_bval', 'in_bval')]),

        (inputnode,      eddy,     [('in_bvec', 'in_bvec')]),
        (inputnode,      eddy,     [('in_bval', 'in_bval')]),
        (inputnode,  eddy,   [('in_file', 'in_file')]),
        (inputnode,     eddy,   [('in_mask', 'in_mask')]),
        (generate_acq,      eddy, [('out_file', 'in_acqp')]),
        (generate_index,      eddy, [('eddy_index', 'in_index')]),

        (eddy,   outputnode, [('out_parameter', 'out_parameter')]),
        (eddy,      outputnode, [('out_corrected', 'out_corrected')]),
        (eddy,      outputnode, [('out_rotated_bvecs', 'out_rotated_bvecs')])
    ])
    return wf
