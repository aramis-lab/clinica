# coding: utf8


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
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from .dwi_preprocessing_using_t1_utils import eddy_fsl, generate_acq, generate_index, b0_indices

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


def epi_pipeline(name='susceptibility_distortion_correction_using_t1'):
    """
    This workflow allows to correct for echo-planar induced susceptibility artifacts without fieldmap
    (e.g. ADNI Database) by elastically register DWIs to their respective baseline T1-weighted
    structural scans using an inverse consistent registration algorithm with a mutual information cost
    function (SyN algorithm). This workflow allows also a coregistration of DWIs with their respective
    baseline T1-weighted structural scans in order to latter combine tracks and cortex parcelation.
    ..  warning:: This workflow rotates the `b`-vectors'
    .. References
      .. Nir et al. (Neurobiology of Aging 2015)- Connectivity network measures predict volumetric atrophy in mild cognitive impairment

        Leow et al. (IEEE Trans Med Imaging 2007)- Statistical Properties of Jacobian Maps and the Realization of Unbiased Large Deformation Nonlinear Image Registration
    Example
    -------
    >>> epi = epi_pipeline()
    >>> epi.inputs.inputnode.DWI = 'DWI.nii'
    >>> epi.inputs.inputnode.bvec = 'bvec.txt'
    >>> epi.inputs.inputnode.T1 = 'T1.nii'
    >>> epi.run() # doctest: +SKIP
    """

    from .dwi_preprocessing_using_t1_utils import (create_jacobian_determinant_image,
                                                   change_itk_transform_type,
                                                   expend_matrix_list,
                                                   rotate_bvecs,
                                                   ants_registration_syn_quick,
                                                   ants_warp_image_multi_transform,
                                                   ants_combin_transform)
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.c3 as c3

    inputnode = pe.Node(niu.IdentityInterface(fields=['T1', 'DWI', 'bvec']), name='inputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_2_T1 = pe.Node(interface=fsl.FLIRT(dof=6), name='flirt_B0_2_T1')
    flirt_b0_2_T1.inputs.interp = "spline"
    flirt_b0_2_T1.inputs.cost = 'normmi'
    flirt_b0_2_T1.inputs.cost_func = 'normmi'

    apply_xfm = pe.Node(interface=fsl.preprocess.ApplyXFM(), name='apply_xfm')
    apply_xfm.inputs.apply_xfm = True

    expend_matrix = pe.Node(
            interface=niu.Function(
                input_names=['in_matrix', 'in_bvec'],
                output_names=['out_matrix_list'],
                function=expend_matrix_list),
            name='expend_matrix')

    rot_bvec = pe.Node(
            niu.Function(
                input_names=['in_matrix', 'in_bvec'],
                output_names=['out_file'],
                function=rotate_bvecs),
            name='Rotate_Bvec')

    antsRegistrationSyNQuick = pe.Node(
            interface=niu.Function(
                input_names=['fix_image', 'moving_image'],
                output_names=['image_warped',
                              'affine_matrix',
                              'warp',
                              'inverse_warped',
                              'inverse_warp'],
                function=ants_registration_syn_quick),
            name='antsRegistrationSyNQuick')

    c3d_flirt2ants = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk')
    c3d_flirt2ants.inputs.itk_transform = True
    c3d_flirt2ants.inputs.fsl2ras = True

    change_transform = pe.Node(niu.Function(
            input_names=['input_affine_file'],
            output_names=['updated_affine_file'],
            function=change_itk_transform_type),
            name='change_transform_type')

    merge_transform = pe.Node(niu.Merge(3), name='MergeTransforms')

    apply_transform = pe.MapNode(interface=niu.Function(input_names=['fix_image', 'moving_image', 'ants_warp_affine'],
                                                        output_names=['out_warp_field', 'out_warped'],
                                                        function=ants_combin_transform),
                                 iterfield=['moving_image'],
                                 name='warp_filed')

    jacobian = pe.MapNode(interface=niu.Function(input_names=['imageDimension', 'deformationField', 'outputImage'],
                                                 output_names=['outputImage'],
                                                 function=create_jacobian_determinant_image),
                          iterfield=['deformationField'],
                          name='jacobian')

    jacobian.inputs.imageDimension = 3
    jacobian.inputs.outputImage = 'Jacobian_image.nii.gz'

    jacmult = pe.MapNode(fsl.MultiImageMaths(op_string='-mul %s'),
                         iterfield=['in_file', 'operand_files'],
                         name='ModulateDWIs')

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    outputnode = pe.Node(niu.IdentityInterface(fields=['DWI_2_T1_Coregistration_matrix',
                                                       'epi_correction_deformation_field',
                                                       'epi_correction_affine_transform',
                                                       'epi_correction_image_warped',
                                                       'DWIs_epicorrected',
                                                       'warp_epi',
                                                       'out_bvec'
                                                       ]), name='outputnode')

    wf = pe.Workflow(name='epi_pipeline')

    wf.connect([(inputnode, split, [('DWI', 'in_file')])])
    wf.connect([(split, pick_ref, [('out_files', 'inlist')])])
    wf.connect([(pick_ref, flirt_b0_2_T1, [('out', 'in_file')])])
    wf.connect([(inputnode, flirt_b0_2_T1, [('T1', 'reference')])])
    wf.connect([(inputnode, rot_bvec, [('bvec', 'in_bvec')])])
    wf.connect([(flirt_b0_2_T1, expend_matrix, [('out_matrix_file', 'in_matrix')])])
    wf.connect([(inputnode, expend_matrix, [('bvec', 'in_bvec')])])
    wf.connect([(expend_matrix, rot_bvec, [('out_matrix_list', 'in_matrix')])])
    wf.connect([(inputnode, antsRegistrationSyNQuick, [('T1', 'fix_image')])])
    wf.connect([(flirt_b0_2_T1, antsRegistrationSyNQuick, [('out_file', 'moving_image')])])

    wf.connect([(inputnode, c3d_flirt2ants, [('T1', 'reference_file')])])
    wf.connect([(pick_ref, c3d_flirt2ants, [('out', 'source_file')])])
    wf.connect([(flirt_b0_2_T1, c3d_flirt2ants, [('out_matrix_file', 'transform_file')])])
    wf.connect([(c3d_flirt2ants, change_transform, [('itk_transform', 'input_affine_file')])])

    wf.connect([(antsRegistrationSyNQuick, merge_transform, [('warp', 'in1')])])
    wf.connect([(antsRegistrationSyNQuick, merge_transform, [('affine_matrix', 'in2')])])
    wf.connect([(change_transform, merge_transform, [('updated_affine_file', 'in3')])])
    wf.connect([(inputnode, apply_transform, [('T1', 'fix_image')])])
    wf.connect([(split, apply_transform, [('out_files', 'moving_image')])])

    wf.connect([(merge_transform, apply_transform, [('out', 'ants_warp_affine')])])
    wf.connect([(apply_transform, jacobian, [('out_warp_field', 'deformationField')])])
    wf.connect([(apply_transform, jacmult, [('out_warped', 'operand_files')])])
    wf.connect([(jacobian, jacmult, [('outputImage', 'in_file')])])
    wf.connect([(jacmult, thres, [('out_file', 'in_file')])])
    wf.connect([(thres, merge, [('out_file', 'in_files')])])

    wf.connect([(merge, outputnode, [('merged_file', 'DWIs_epicorrected')])])
    wf.connect([(flirt_b0_2_T1, outputnode, [('out_matrix_file', 'DWI_2_T1_Coregistration_matrix')])])
    wf.connect([(antsRegistrationSyNQuick, outputnode, [('warp', 'epi_correction_deformation_field'),
                                                        ('affine_matrix', 'epi_correction_affine_transform'),
                                                        ('image_warped', 'epi_correction_image_warped')])])
    wf.connect([(merge_transform, outputnode, [('out', 'warp_epi')])])
    wf.connect([(rot_bvec, outputnode, [('out_file', 'out_bvec')])])

    return wf
