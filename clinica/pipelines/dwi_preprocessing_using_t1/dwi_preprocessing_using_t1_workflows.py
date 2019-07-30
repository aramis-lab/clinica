# coding: utf8


def eddy_fsl_pipeline(low_bval, name='eddy_fsl'):
    """
    Using eddy from FSL for head motion correction and eddy current distortion correction.
    """
    from nipype.interfaces.fsl import Eddy
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from clinica.utils.dwi import generate_acq_file, generate_index_file

    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['in_file',
                        'in_bvec',
                        'in_bval',
                        'in_mask',
                        'ref_b0',
                        'total_readout_time',
                        'phase_encoding_direction']),
            name='inputnode')

    generate_acq = pe.Node(niu.Function(input_names=['in_dwi', 'fsl_phase_encoding_direction', 'total_readout_time'],
                                        output_names=['out_file'],
                                        function=generate_acq_file),
                           name='generate_acq')

    generate_index = pe.Node(niu.Function(input_names=['in_bval', 'low_bval'],
                                          output_names=['out_file'],
                                          function=generate_index_file),
                             name='generate_index')
    generate_index.inputs.low_bval = low_bval

    eddy = pe.Node(interface=Eddy(), name='eddy_fsl')
    eddy.inputs.flm = 'linear'
    eddy.inputs.use_cuda = True

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_parameter',
                                                       'out_corrected',
                                                       'out_rotated_bvecs']),
                         name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, generate_acq, [('in_file', 'in_dwi')]),
        (inputnode, generate_acq, [('total_readout_time', 'total_readout_time')]),
        (inputnode, generate_acq, [('phase_encoding_direction', 'fsl_phase_encoding_direction')]),

        (inputnode, generate_index, [('in_bval', 'in_bval')]),

        (inputnode,  eddy,     [('in_bvec', 'in_bvec')]),
        (inputnode,  eddy,     [('in_bval', 'in_bval')]),
        (inputnode,  eddy,   [('in_file', 'in_file')]),
        (inputnode,  eddy,   [('in_mask', 'in_mask')]),
        (generate_acq, eddy, [('out_file', 'in_acqp')]),
        (generate_index, eddy, [('out_file', 'in_index')]),
        (eddy, outputnode, [('out_parameter', 'out_parameter')]),
        (eddy, outputnode, [('out_corrected', 'out_corrected')]),
        (eddy, outputnode, [('out_rotated_bvecs', 'out_rotated_bvecs')])
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
                                                   ants_combin_transform)
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants
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

    ants_registration = pe.Node(interface=ants.registration.RegistrationSynQuick(
        transform_type='br',
        dimension=3),
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

    wf.connect([
        (inputnode, split, [('DWI', 'in_file')]),
        (split, pick_ref, [('out_files', 'inlist')]),
        (pick_ref, flirt_b0_2_T1, [('out', 'in_file')]),
        (inputnode, flirt_b0_2_T1, [('T1', 'reference')]),
        (inputnode, rot_bvec, [('bvec', 'in_bvec')]),
        (flirt_b0_2_T1, expend_matrix, [('out_matrix_file', 'in_matrix')]),
        (inputnode, expend_matrix, [('bvec', 'in_bvec')]),
        (expend_matrix, rot_bvec, [('out_matrix_list', 'in_matrix')]),
        (inputnode, ants_registration, [('T1', 'fixed_image')]),
        (flirt_b0_2_T1, ants_registration, [('out_file', 'moving_image')]),

        (inputnode, c3d_flirt2ants, [('T1', 'reference_file')]),
        (pick_ref, c3d_flirt2ants, [('out', 'source_file')]),
        (flirt_b0_2_T1, c3d_flirt2ants, [('out_matrix_file', 'transform_file')]),
        (c3d_flirt2ants, change_transform, [('itk_transform', 'input_affine_file')]),

        (ants_registration, merge_transform, [('forward_warp_field', 'in1')]),
        (ants_registration, merge_transform, [('out_matrix', 'in2')]),
        (change_transform, merge_transform, [('updated_affine_file', 'in3')]),
        (inputnode, apply_transform, [('T1', 'fix_image')]),
        (split, apply_transform, [('out_files', 'moving_image')]),

        (merge_transform, apply_transform, [('out', 'ants_warp_affine')]),
        (apply_transform, jacobian, [('out_warp_field', 'deformationField')]),
        (apply_transform, jacmult, [('out_warped', 'operand_files')]),
        (jacobian, jacmult, [('outputImage', 'in_file')]),
        (jacmult, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),

        (merge, outputnode, [('merged_file', 'DWIs_epicorrected')]),
        (flirt_b0_2_T1, outputnode, [('out_matrix_file', 'DWI_2_T1_Coregistration_matrix')]),
        (ants_registration, outputnode, [('forward_warp_field', 'epi_correction_deformation_field'),
                                         ('out_matrix', 'epi_correction_affine_transform'),
                                         ('warped_image', 'epi_correction_image_warped')]),
        (merge_transform, outputnode, [('out', 'warp_epi')]),
        (rot_bvec, outputnode, [('out_file', 'out_bvec')]),
    ])

    return wf
