# coding: utf8


def susceptibility_distortion_correction_using_t1(
        name='susceptibility_distortion_correction_using_t1'):
    """

    Args:
        name (Optional[str]): Name of the workflow.

    Inputnode:

    Outputnode:
        outputnode.out_dwi - corrected dwi file
        outputnode.out_bvec - rotated gradient vectors table
        outputnode.out_b0_to_t1_rigid_body_matrix - B0 to T1 image FLIRT rigid body fsl coregistration matrix
        outputnode.out_t1_coregistered_to_b0 - T1 image rigid body coregistered to the B0 image
        outputnode.out_b0_to_t1_affine_matrix - B0 to T1 image ANTs affine itk coregistration matrix
        outputnode.out_b0_to_t1_syn_deformation_field - B0 to T1 image ANTs SyN itk warp
        outputnode.out_warp - Out warp allowing DWI to T1 registration and susceptibilty induced artifacts correction



    Returns:

    Example:
        >>> epi = susceptibility_distortion_correction_using_t1()
        >>> epi.inputs.inputnode.in_dwi = 'dwi.nii'
        >>> epi.inputs.inputnode.in_t1 = 'T1w.nii'
        >>> epi.run() # doctest: +SKIP
    """

    """
    SDC stands for susceptibility distortion correction and SYB stand for SyN based. This workflow
    allows to correct for echo-planar induced susceptibility artifacts without fieldmap
    (e.g. ADNI Database) by elastically register DWIs to their respective baseline T1-weighted
    structural scans using an inverse consistent registration algorithm with a mutual information cost
    function (SyN algorithm).
    .. References
      .. Nir et al. (Neurobiology of Aging 2015)- Connectivity network measures predict volumetric atrophy in mild cognitive impairment

         Leow et al. (IEEE Trans Med Imaging 2007)- Statistical Properties of Jacobian Maps and the Realization of Unbiased Large Deformation Nonlinear Image Registration

    Inputnode
    ---------
    DWI : FILE
      Mandatory input. Input dwi file.
    T1 : FILE
      Mandatory input. Input T1 file.

    Outputnode
    ----------


    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    import clinica.pipeline.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils as utils

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_t1', 'in_dwi']),
                        name='inputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_to_t1 = pe.Node(interface=fsl.FLIRT(dof=6), name='flirt_b0_to_t1')
    flirt_b0_to_t1.inputs.interp = "spline"
    flirt_b0_to_t1.inputs.cost = 'normmi'
    flirt_b0_to_t1.inputs.cost_func = 'normmi'

    invert_xfm = pe.Node(interface=fsl.ConvertXFM(), name='invert_xfm')
    invert_xfm.inputs.invert_xfm = True

    apply_xfm = pe.Node(interface=fsl.ApplyXfm(), name='apply_xfm')
    apply_xfm.inputs.apply_xfm = True
    apply_xfm.inputs.interp = 'spline'
    apply_xfm.inputs.cost = 'normmi'
    apply_xfm.inputs.cost_func = 'normmi'

    ants_registration_syn_quick = pe.Node(interface=niu.Function(
        input_names=['fix_image', 'moving_image'],
        output_names=['image_warped', 'affine_matrix',
                      'warp', 'inverse_warped', 'inverse_warp'],
        function=utils.ants_registration_syn_quick),
        name='ants_registration_syn_quick')

    merge_transform = pe.Node(niu.Merge(2), name='MergeTransforms')

    combine_warp = pe.Node(interface=niu.Function(
        input_names=['in_file', 'transforms_list', 'reference'],
        output_names=['out_warp'],
        function=utils.ants_combine_transform), name='combine_warp')

    coeffs = pe.Node(fsl.WarpUtils(out_format='spline'), name='CoeffComp')

    fsl_transf = pe.Node(fsl.WarpUtils(out_format='field'), name='fsl_transf')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(), iterfield=['in_file'],
                            name='apply_warp')
    apply_warp.inputs.interp = 'spline'

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_b0_to_t1_rigid_body_matrix',
                'out_t1_to_b0_rigid_body_matrix',
                'out_t1_coregistered_to_b0',
                'out_b0_to_t1_syn_deformation_field',
                'out_b0_to_t1_affine_matrix',
                'out_dwi',
                'out_warp']),
        name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, split, [('in_dwi', 'in_file')]),
        (split, pick_ref, [('out_files', 'inlist')]),
        (pick_ref, flirt_b0_to_t1, [('out', 'in_file')]),
        (inputnode, flirt_b0_to_t1, [('in_t1', 'reference')]),
        (flirt_b0_to_t1, invert_xfm, [('out_matrix_file', 'in_file')]),
        (invert_xfm, apply_xfm, [('out_file', 'in_matrix_file')]),
        (inputnode, apply_xfm, [('in_t1', 'in_file')]),
        (pick_ref, apply_xfm, [('out', 'reference')]),
        (apply_xfm, ants_registration_syn_quick, [('out_file', 'fix_image')]),
        (pick_ref, ants_registration_syn_quick, [('out', 'moving_image')]),
        (ants_registration_syn_quick, merge_transform, [('affine_matrix', 'in2'),
                                                        ('warp', 'in1')]),
        (pick_ref, combine_warp, [('out', 'in_file')]),
        (merge_transform, combine_warp, [('out', 'transforms_list')]),
        (apply_xfm, combine_warp, [('out_file', 'reference')]),
        (apply_xfm, coeffs, [('out_file', 'reference')]),
        (combine_warp, coeffs, [('out_warp', 'in_file')]),
        (coeffs, fsl_transf, [('out_file', 'in_file')]),
        (apply_xfm, fsl_transf, [('out_file', 'reference')]),
        (fsl_transf, apply_warp, [('out_file', 'field_file')]),
        (split, apply_warp, [('out_files', 'in_file')]),
        (apply_xfm, apply_warp, [('out_file', 'ref_file')]),
        (apply_warp, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),
        (merge,                       outputnode, [('merged_file',                            'out_dwi')]),  # noqa
        (flirt_b0_to_t1,              outputnode, [('out_matrix_file', 'out_b0_to_t1_rigid_body_matrix')]),  # noqa
        (invert_xfm,                  outputnode, [('out_file',        'out_t1_to_b0_rigid_body_matrix')]),  # noqa
        (apply_xfm,                   outputnode, [('out_file',             'out_t1_coregistered_to_b0')]),  # noqa
        (ants_registration_syn_quick, outputnode, [('warp',          'out_b0_to_t1_syn_deformation_field'),  # noqa
                                                   ('affine_matrix',       'out_b0_to_t1_affine_matrix')]),  # noqa
        (fsl_transf,                  outputnode, [('out_file',                              'out_warp')])   # noqa
    ])

    return wf


def apply_all_corrections_using_ants(name='apply_all_corrections_using_ants'):
    """
    Combines two lists of linear transforms with the deformation field
    map obtained epi_correction by Ants.
    Additionally, computes the corresponding bspline coefficients and
    the map of determinants of the jacobian.
    """
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_sdc_syb', 'in_hmc', 'in_ecc', 'in_dwi', 'in_t1']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_warp', 'out_coeff', 'out_jacobian']),
        name='outputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')

    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_to_t1 = pe.Node(interface=fsl.FLIRT(dof=6), name='flirt_b0_to_t1')
    flirt_b0_to_t1.inputs.interp = "spline"
    flirt_b0_to_t1.inputs.cost = 'normmi'
    flirt_b0_to_t1.inputs.cost_func = 'normmi'

    invert_xfm = pe.Node(interface=fsl.ConvertXFM(), name='invert_xfm')
    invert_xfm.inputs.invert_xfm = True

    apply_xfm = pe.Node(interface=fsl.ApplyXfm(), name='apply_xfm')
    apply_xfm.inputs.apply_xfm = True
    apply_xfm.inputs.interp = "spline"
    apply_xfm.inputs.cost = 'normmi'
    apply_xfm.inputs.cost_func = 'normmi'

    concat_hmc_ecc = pe.MapNode(fsl.ConvertXFM(),
                                iterfield=['in_file', 'in_file2'],
                                name="concat_hmc_ecc")
    concat_hmc_ecc.inputs.concat_xfm = True

    warps = pe.MapNode(fsl.ConvertWarp(),
                       iterfield=['premat'],
                       name='ConvertWarp')

    unwarp = pe.MapNode(interface=fsl.ApplyWarp(),
                        iterfield=['in_file', 'field_file'],
                        name='unwarp_warp')
    unwarp.inputs.interp = 'spline'

    coeffs = pe.MapNode(fsl.WarpUtils(out_format='spline'),
                        iterfield=['in_file'], name='CoeffComp')
    jacobian = pe.MapNode(fsl.WarpUtils(write_jacobian=True),
                          iterfield=['in_file'], name='JacobianComp')
    jacmult = pe.MapNode(fsl.MultiImageMaths(op_string='-mul %s'),
                         iterfield=['in_file', 'operand_files'],
                         name='ModulateDWIs')

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')
    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, concat_hmc_ecc, [('in_ecc', 'in_file2')]),
        (inputnode, concat_hmc_ecc, [('in_hmc', 'in_file')]),
        (concat_hmc_ecc, warps,     [('out_file', 'premat')]),
        (inputnode, warps,        [('in_sdc_syb', 'warp1')]),
        (inputnode, split, [('in_dwi', 'in_file')]),
        (split, pick_ref, [('out_files', 'inlist')]),
        (pick_ref, flirt_b0_to_t1, [('out', 'in_file')]),
        (inputnode, flirt_b0_to_t1, [('in_t1', 'reference')]),
        (flirt_b0_to_t1, invert_xfm, [('out_matrix_file', 'in_file')]),
        (invert_xfm, apply_xfm, [('out_file', 'in_matrix_file')]),
        (inputnode, apply_xfm, [('in_t1', 'in_file')]),
        (pick_ref, apply_xfm, [('out', 'reference')]),
        (apply_xfm, warps, [('out_file', 'reference')]),
        (warps, unwarp, [('out_file', 'field_file')]),
        (split, unwarp, [('out_files', 'in_file')]),
        (apply_xfm, unwarp, [('out_file', 'ref_file')]),
        (apply_xfm, coeffs, [('out_file', 'reference')]),
        (warps, coeffs, [('out_file', 'in_file')]),
        (apply_xfm, jacobian, [('out_file', 'reference')]),
        (coeffs, jacobian, [('out_file', 'in_file')]),
        (unwarp, jacmult, [('out_file', 'in_file')]),
        (jacobian, jacmult, [('out_jacobian', 'operand_files')]),
        (jacmult, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),
        (warps, outputnode, [('out_file', 'out_warp')]),
        (coeffs, outputnode, [('out_file', 'out_coeff')]),
        (jacobian, outputnode, [('out_jacobian', 'out_jacobian')]),
        (merge, outputnode, [('merged_file', 'out_file')])
    ])
    return wf
