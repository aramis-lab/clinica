# coding: utf8


def susceptibility_distortion_correction_using_t1(
        name='susceptibility_distortion_correction_using_t1'):
    """
    Susceptibility distortion correction using the T1w image.

    This workflow allows to correct for echo-planar induced susceptibility
    artifacts without fieldmap (e.g. ADNI Database) by elastically register
    DWIs to their respective baseline T1-weighted structural scans using an
    inverse consistent registration algorithm with a mutual information cost
    function (SyN algorithm).

    Args:
        name (Optional[str]): Name of the workflow.

    Inputnode:
        in_t1 (str): T1w image.
        in_dwi (str): DWI dataset

    Outputnode:
        out_dwi (str): Corrected DWI dataset
        out_warp (str): Out warp allowing DWI to T1 registration and
            susceptibilty induced artifacts correction
        out_b0_to_t1_rigid_body_matrix (str): B0 to T1 image FLIRT rigid body
            FSL coregistration matrix
        out_t1_to_b0_rigid_body_matrix (str): T1 to B0 image FLIRT rigid body
            FSL coregistration matrix
        out_t1_coregistered_to_b0 (str): T1 image rigid body coregistered to
            the B0 image
        out_b0_to_t1_syn_deformation_field (str): B0 to T1 image ANTs SyN
            ITK warp
        out_b0_to_t1_affine_matrix (str): B0 to T1 image ANTs affine ITK
            coregistration matrix

    References:
      .. Nir et al. (Neurobiology of Aging 2015): Connectivity network measures
        predict volumetric atrophy in mild cognitive impairment

      .. Leow et al. (IEEE Trans Med Imaging 2007): Statistical Properties of
        Jacobian Maps and the Realization of Unbiased Large Deformation
        Nonlinear Image Registration


    Returns:
        The workflow

    Example:
        >>> epi = susceptibility_distortion_correction_using_t1()
        >>> epi.inputs.inputnode.in_dwi = 'dwi.nii'
        >>> epi.inputs.inputnode.in_t1 = 'T1w.nii'
        >>> epi.run() # doctest: +SKIP
    """
    import nipype
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    import \
        clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils as utils

    def expend_matrix_list(in_matrix, in_bvec):
        import numpy as np

        bvecs = np.loadtxt(in_bvec).T
        out_matrix_list = [in_matrix]

        out_matrix_list = out_matrix_list * len(bvecs)

        return out_matrix_list

    def rotate_bvecs(in_bvec, in_matrix):
        """
        Rotates the input bvec file accordingly with a list of matrices.
        .. note:: the input affine matrix transforms points in the destination
          image to their corresponding coordinates in the original image.
          Therefore, this matrix should be inverted first, as we want to know
          the target position of :math:`\\vec{r}`.
        """
        import os
        import numpy as np

        name, fext = os.path.splitext(os.path.basename(in_bvec))
        if fext == '.gz':
            name, _ = os.path.splitext(name)
        out_file = os.path.abspath('%s_rotated.bvec' % name)
        bvecs = np.loadtxt(
            in_bvec).T  # Warning, bvecs.txt are not in the good configuration, need to put '.T'
        new_bvecs = []

        if len(bvecs) != len(in_matrix):
            raise RuntimeError(('Number of b-vectors (%d) and rotation '
                                'matrices (%d) should match.') %
                               (len(bvecs), len(in_matrix)))

        for bvec, mat in zip(bvecs, in_matrix):
            if np.all(bvec == 0.0):
                new_bvecs.append(bvec)
            else:
                invrot = np.linalg.inv(np.loadtxt(mat))[:3, :3]
                newbvec = invrot.dot(bvec)
                new_bvecs.append((newbvec / np.linalg.norm(newbvec)))

        np.savetxt(out_file, np.array(new_bvecs).T, fmt='%0.15f')
        return out_file

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['in_t1', 'in_dwi', 'in_bvec']),
        name='inputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_to_t1 = pe.Node(interface=fsl.FLIRT(dof=6),
                             name='flirt_b0_to_t1')
    flirt_b0_to_t1.inputs.interp = "spline"
    flirt_b0_to_t1.inputs.cost = 'normmi'
    flirt_b0_to_t1.inputs.cost_func = 'normmi'

    if nipype.__version__.split('.') < ['0', '13', '0']:
        apply_xfm = pe.Node(interface=fsl.ApplyXfm(),
                            name='apply_xfm')
    else:
        apply_xfm = pe.Node(interface=fsl.ApplyXFM(),
                            name='apply_xfm')
    apply_xfm.inputs.apply_xfm = True

    expend_matrix = pe.Node(
        interface=niu.Function(input_names=['in_matrix', 'in_bvec'],
                               output_names=['out_matrix_list'],
                               function=expend_matrix_list),
        name='expend_matrix')

    rot_bvec = pe.Node(niu.Function(input_names=['in_matrix', 'in_bvec'],
                                    output_names=['out_file'],
                                    function=rotate_bvecs),
                       name='Rotate_Bvec')

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

    fsl_transf = pe.Node(fsl.WarpUtils(out_format='field'),
                         name='fsl_transf')

    warp_epi = pe.Node(fsl.ConvertWarp(), name='warp_epi')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                            iterfield=['in_file'], name='apply_warp')
    apply_warp.inputs.interp = 'spline'

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['dwi_to_t1_coregistration_matrix',
                'itk_dwi_t1_coregistration_matrix',
                'epi_correction_deformation_field',
                'epi_correction_affine_transform',
                'merge_epi_transform', 'out_dwi', 'out_warp',
                'out_bvec']), name='outputnode')

    wf = pe.Workflow(name=name)

    wf.connect([
        (inputnode, split, [('in_dwi', 'in_file')]),  # noqa

        (split, pick_ref, [('out_files', 'inlist')]),  # noqa

        (pick_ref, flirt_b0_to_t1, [('out', 'in_file')]),  # noqa
        (inputnode, flirt_b0_to_t1, [('in_t1', 'reference')]),  # noqa

        (flirt_b0_to_t1, expend_matrix, [('out_matrix_file', 'in_matrix')]),
        # noqa
        (inputnode, expend_matrix, [('in_bvec', 'in_bvec')]),  # noqa

        (inputnode, rot_bvec, [('in_bvec', 'in_bvec')]),  # noqa
        (expend_matrix, rot_bvec, [('out_matrix_list', 'in_matrix')]),  # noqa

        (inputnode, ants_registration_syn_quick, [('in_t1', 'fix_image')]),
        # noqa
        (flirt_b0_to_t1, ants_registration_syn_quick,
         [('out_file', 'moving_image')]),  # noqa

        (ants_registration_syn_quick, merge_transform,
         [('affine_matrix', 'in2'),  # noqa
          ('warp', 'in1')]),  # noqa

        (flirt_b0_to_t1, combine_warp, [('out_file', 'in_file')]),  # noqa
        (merge_transform, combine_warp, [('out', 'transforms_list')]),  # noqa
        (inputnode, combine_warp, [('in_t1', 'reference')]),  # noqa

        (inputnode, coeffs, [('in_t1', 'reference')]),  # noqa
        (combine_warp, coeffs, [('out_warp', 'in_file')]),  # noqa

        (coeffs, fsl_transf, [('out_file', 'in_file')]),  # noqa
        (inputnode, fsl_transf, [('in_t1', 'reference')]),  # noqa

        (inputnode, warp_epi, [('in_t1', 'reference')]),  # noqa
        (flirt_b0_to_t1, warp_epi, [('out_matrix_file', 'premat')]),  # noqa
        (fsl_transf, warp_epi, [('out_file', 'warp1')]),  # noqa

        (warp_epi, apply_warp, [('out_file', 'field_file')]),  # noqa
        (split, apply_warp, [('out_files', 'in_file')]),  # noqa
        (inputnode, apply_warp, [('in_t1', 'ref_file')]),  # noqa

        (apply_warp, thres, [('out_file', 'in_file')]),  # noqa

        (thres, merge, [('out_file', 'in_files')]),  # noqa
        # Outputnode
        (merge, outputnode, [('merged_file', 'out_dwi')]),  # noqa
        (flirt_b0_to_t1, outputnode,
         [('out_matrix_file', 'dwi_to_t1_coregistration_matrix')]),  # noqa
        (ants_registration_syn_quick, outputnode,
         [('warp', 'epi_correction_deformation_field'),  # noqa
          ('affine_matrix', 'epi_correction_affine_transform')]),  # noqa
        (warp_epi, outputnode, [('out_file', 'out_warp')]),  # noqa
        (rot_bvec, outputnode, [('out_file', 'out_bvec')]),  # noqa
    ])
    return wf


def apply_all_corrections_using_ants(name='UnwarpArtifacts'):
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

    concat_hmc_ecc = pe.MapNode(fsl.ConvertXFM(), name="concat_hmc_ecc",
                                iterfield=['in_file', 'in_file2'])
    concat_hmc_ecc.inputs.concat_xfm = True

    warps = pe.MapNode(fsl.ConvertWarp(), iterfield=['premat'],
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
        (inputnode, concat_hmc_ecc, [('in_ecc', 'in_file2')]),  # noqa
        (inputnode, concat_hmc_ecc, [('in_hmc', 'in_file')]),  # noqa

        (concat_hmc_ecc, warps, [('out_file', 'premat')]),  # noqa
        (inputnode, warps, [('in_sdc_syb', 'warp1')]),  # noqa
        (inputnode, warps, [('in_t1', 'reference')]),  # noqa
        (inputnode, split, [('in_dwi', 'in_file')]),  # noqa
        (warps, unwarp, [('out_file', 'field_file')]),  # noqa
        (split, unwarp, [('out_files', 'in_file')]),  # noqa
        (inputnode, unwarp, [('in_t1', 'ref_file')]),  # noqa
        (inputnode, coeffs, [('in_t1', 'reference')]),  # noqa
        (warps, coeffs, [('out_file', 'in_file')]),  # noqa
        (inputnode, jacobian, [('in_t1', 'reference')]),  # noqa
        (coeffs, jacobian, [('out_file', 'in_file')]),  # noqa
        (unwarp, jacmult, [('out_file', 'in_file')]),  # noqa
        (jacobian, jacmult, [('out_jacobian', 'operand_files')]),  # noqa
        (jacmult, thres, [('out_file', 'in_file')]),  # noqa
        (thres, merge, [('out_file', 'in_files')]),  # noqa
        (warps, outputnode, [('out_file', 'out_warp')]),  # noqa
        (coeffs, outputnode, [('out_file', 'out_coeff')]),  # noqa
        (jacobian, outputnode, [('out_jacobian', 'out_jacobian')]),  # noqa
        (merge, outputnode, [('merged_file', 'out_file')])  # noqa
    ])

    return wf
