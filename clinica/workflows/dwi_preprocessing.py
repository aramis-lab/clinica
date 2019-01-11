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


def dwi_flirt(name='DWICoregistration', excl_nodiff=False, flirt_param={}):
    """
    Generates a workflow for linear registration of dwi volumes using flirt.

    Inputnode
    ---------
    reference : FILE
      Mandatory input. Reference data set.
    in_file : FILE
      Mandatory input. Moving data set.
    ref_mask : FILE
      Mandatory input. Binary mask of the reference volume.
    in_xfms : FILE
      Mandatory input. Intialisation matrices for flirt.
    in_bval : FILE
      Mandatory input. B values file.

    """
    import nipype.interfaces.ants as ants
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from nipype.workflows.dmri.fsl.utils import _checkinitxfm

    from nipype.workflows.dmri.fsl.utils import enhance

    inputnode = pe.Node(
            niu.IdentityInterface(
                fields=['reference',
                        'in_file',
                        'ref_mask',
                        'in_xfms',
                        'in_bval']),
            name='inputnode')

    initmat = pe.Node(
            niu.Function(
                input_names=['in_bval',
                             'in_xfms',
                             'excl_nodiff'],
                output_names=['init_xfms'],
                function=_checkinitxfm),
            name='InitXforms')
    initmat.inputs.excl_nodiff = excl_nodiff
    dilate = pe.Node(
            fsl.maths.MathsCommand(
                nan2zeros=True,
                args='-kernel sphere 5 -dilM'),
            name='MskDilate')
    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='Bias')
    flirt = pe.MapNode(fsl.FLIRT(**flirt_param), name='CoRegistration',
                       iterfield=['in_file', 'in_matrix_file'])
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')
    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
    outputnode = pe.Node(
            niu.IdentityInterface(
                fields=['out_file',
                        'out_xfms',
                        'out_ref']),
            name='outputnode')
    enhb0 = pe.Node(niu.Function(
        input_names=['in_file', 'in_mask', 'clip_limit'],
        output_names=['out_file'], function=enhance), name='B0Equalize')
    enhb0.inputs.clip_limit = 0.015
    enhdw = pe.MapNode(niu.Function(
        input_names=['in_file', 'in_mask'], output_names=['out_file'],
        function=enhance), name='DWEqualize', iterfield=['in_file'])
    # enhb0.inputs.clip_limit = clip_limit

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  split,      [('in_file', 'in_file')]),
        (inputnode,  dilate,     [('ref_mask', 'in_file')]),
        (inputnode,   n4,        [('reference', 'input_image'),
                                  ('ref_mask', 'mask_image')]),
        #        (inputnode,  flirt,      [('ref_mask', 'reference')]),
        (n4, enhb0, [('output_image', 'in_file')]),
        (enhb0, flirt, [('out_file', 'reference')]),
        (inputnode, initmat, [('in_xfms', 'in_xfms'),
                              ('in_bval', 'in_bval')]),
        (split, enhdw, [('out_files', 'in_file')]),
        (dilate, enhdw, [('out_file', 'in_mask')]),
        (dilate, flirt, [('out_file', 'ref_weight'),
                         ('out_file', 'in_weight')]),
        (enhdw, flirt, [('out_file', 'in_file')]),
        (initmat, flirt, [('init_xfms', 'in_matrix_file')]),
        (flirt,      thres,      [('out_file', 'in_file')]),
        (thres,      merge,      [('out_file', 'in_files')]),
        (merge,     outputnode, [('merged_file', 'out_file')]),
        (enhb0, outputnode, [('out_file', 'out_ref')]),
        (flirt,     outputnode, [('out_matrix_file', 'out_xfms')])
    ])
    return wf


def hmc_pipeline(name='motion_correct'):
    """
    HMC stands for head-motion correction.

    Creates a pipelines that corrects for head motion artifacts in dMRI
    sequences. It takes a series of diffusion weighted images and
    rigidly co-registers them to one reference image (FLIRT normalised
    mutual information). Finally, the `b`-matrix is rotated
    accordingly [Leemans09]_ making use of the rotation matrix
    obtained by FLIRT.

    A list of rigid transformation matrices is provided, so that transforms
    can be chained.
    This is useful to correct for artifacts with only one interpolation process
    and also to compute nuisance regressors as proposed by [Yendiki13]_.

    .. warning:: This workflow rotates the `b`-vectors, so please be advised
      that not all the dicom converters ensure the consistency between the
      resulting nifti orientation and the gradients table (e.g. dcm2nii
      checks it).

    .. admonition:: References

      .. [Leemans09] Leemans A, and Jones DK, `The B-matrix must be rotated
        when correcting for subject motion in DTI data
        <http://dx.doi.org/10.1002/mrm.21890>`_,
        Magn Reson Med. 61(6):1336-49. 2009. doi: 10.1002/mrm.21890.

      .. [Yendiki13] Yendiki A et al., `Spurious group differences due to head
        motion in a diffusion MRI study
        <http://dx.doi.org/10.1016/j.neuroimage.2013.11.027>`_.
        Neuroimage. 21(88C):79-90. 2013. doi: 10.1016/j.neuroimage.2013.11.027

    Inputnode
    ---------
    in_file : FILE
      Mandatory input. Input dwi file.
    in_bvec : FILE
      Mandatory input. Vector file of the diffusion directions of the dwi dataset.
    in_bval : FILE
      Mandatory input. B values file.
    in_mask : FILE
      Mandatory input. Weights mask of reference image (a file with data
      range in [0.0, 1.0], indicating the weight of each voxel when computing the metric
    ref_num : INT
      Optional input. Default=0. Index of the b0 volume that should be taken as reference.

    Outputnode
    ----------
        outputnode.out_file - corrected dwi file
        outputnode.out_bvec - rotated gradient vectors table
        outputnode.out_xfms - list of transformation matrices

    """
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from nipype.workflows.data import get_flirt_schedule
    from nipype.workflows.dmri.fsl.utils import insert_mat
    from nipype.workflows.dmri.fsl.utils import rotate_bvecs

    from clinica.utils.dwi import merge_volumes_tdim
    from clinica.utils.dwi import hmc_split
    from clinica.workflows.dwi_preprocessing import dwi_flirt

    params = dict(dof=6, bgvalue=0, save_log=True, no_search=True,
                  # cost='mutualinfo', cost_func='mutualinfo', bins=64,
                  schedule=get_flirt_schedule('hmc'))

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_bvec', 'in_bval', 'in_mask', 'ref_num']),
        name='inputnode')

    split = pe.Node(
            niu.Function(
                function=hmc_split,
                input_names=['in_file', 'in_bval', 'ref_num'],
                output_names=['out_ref', 'out_mov', 'out_bval', 'volid']),
            name='split_ref_moving')

    flirt = dwi_flirt(flirt_param=params)

    insmat = pe.Node(niu.Function(
        input_names=['inlist', 'volid'],
        output_names=['out'],
        function=insert_mat), name='insert_ref_matrix')

    rot_bvec = pe.Node(niu.Function(
        input_names=['in_bvec', 'in_matrix'],
        output_names=['out_file'],
        function=rotate_bvecs), name='Rotate_Bvec')

    merged_volumes = pe.Node(niu.Function(
        input_names=['in_file1', 'in_file2'],
        output_names=['out_file'],
        function=merge_volumes_tdim), name='merge_reference_moving')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_bvec', 'out_xfms', 'mask_B0']),
        name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,          split,  [('in_file',                'in_file'),
                                      ('in_bval',                'in_bval'),
                                      ('ref_num',              'ref_num')]),
        (inputnode,          flirt,  [('in_mask',   'inputnode.ref_mask')]),
        (split,              flirt,  [('out_ref',    'inputnode.reference'),
                                      ('out_mov',      'inputnode.in_file'),
                                      ('out_bval',   'inputnode.in_bval')]),
        (flirt,              insmat, [('outputnode.out_xfms',   'inlist')]),
        (split,              insmat, [('volid',                  'volid')]),
        (inputnode,        rot_bvec, [('in_bvec',              'in_bvec')]),
        (insmat,           rot_bvec, [('out',                'in_matrix')]),
        (rot_bvec,       outputnode, [('out_file',            'out_bvec')]),
        (flirt,      merged_volumes, [('outputnode.out_ref',    'in_file1'),
                                      ('outputnode.out_file', 'in_file2')]),
        (merged_volumes, outputnode, [('out_file',            'out_file')]),
        (insmat,         outputnode, [('out',                 'out_xfms')]),
        (flirt,          outputnode, [('outputnode.out_ref',   'mask_B0')])
    ])
    return wf


def ecc_pipeline(name='eddy_correct'):
    """
    ECC stands for Eddy currents correction.
    Creates a pipelines that corrects for artifacts induced by Eddy currents in
    dMRI sequences.
    It takes a series of diffusion weighted images and linearly co-registers
    them to one reference image (the average of all b0s in the dataset).
    DWIs are also modulated by the determinant of the Jacobian as indicated by
    [Jones10]_ and [Rohde04]_.
    A list of rigid transformation matrices can be provided, sourcing from a
    :func:`.hmc_pipeline` workflow, to initialize registrations in a *motion
    free* framework.
    A list of affine transformation matrices is available as output, so that
    transforms can be chained (discussion
    `here <https://github.com/nipy/nipype/pull/530#issuecomment-14505042>`_).
    .. admonition:: References
      .. [Jones10] Jones DK, `The signal intensity must be modulated by the
        determinant of the Jacobian when correcting for eddy currents in
        diffusion MRI
        <http://cds.ismrm.org/protected/10MProceedings/files/1644_129.pdf>`_,
        Proc. ISMRM 18th Annual Meeting, (2010).
      .. [Rohde04] Rohde et al., `Comprehensive Approach for Correction of
        Motion and Distortion in Diffusion-Weighted MRI
        <http://stbb.nichd.nih.gov/pdf/com_app_cor_mri04.pdf>`_, MRM
        51:103-114 (2004).
    Example
    -------
    from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline
    ecc = ecc_pipeline()
    ecc.inputs.inputnode.in_file = 'diffusion.nii'
    ecc.inputs.inputnode.in_bval = 'diffusion.bval'
    ecc.inputs.inputnode.in_mask = 'mask.nii'
    ecc.run() # doctest: +SKIP
    Inputs::
        inputnode.in_file - input dwi file
        inputnode.in_mask - weights mask of reference image (a file with data \
range sin [0.0, 1.0], indicating the weight of each voxel when computing the \
metric.
        inputnode.in_bval - b-values table
        inputnode.in_xfms - list of matrices to initialize registration (from \
head-motion correction)
    Outputs::
        outputnode.out_file - corrected dwi file
        outputnode.out_xfms - list of transformation matrices
    """

    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl

    from nipype.workflows.data import get_flirt_schedule
    from nipype.workflows.dmri.fsl.utils import extract_bval
    from nipype.workflows.dmri.fsl.utils import recompose_xfm
    from nipype.workflows.dmri.fsl.utils import recompose_dwi
    from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian

    from clinica.workflows.dwi_preprocessing import dwi_flirt
    from clinica.utils.dwi import merge_volumes_tdim

    params = dict(dof=12, no_search=True, interp='spline', bgvalue=0,
                  schedule=get_flirt_schedule('ecc'))

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_bval', 'in_mask', 'in_xfms']), name='inputnode')

    getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')

    pick_dws = pe.Node(niu.Function(
        input_names=['in_dwi', 'in_bval', 'b'], output_names=['out_file'],
        function=extract_bval), name='extract_dwi')
    pick_dws.inputs.b = 'diff'

    flirt = dwi_flirt(flirt_param=params, excl_nodiff=True)

    mult = pe.MapNode(fsl.BinaryMaths(operation='mul'), name='ModulateDWIs',
                      iterfield=['in_file', 'operand_value'])
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    get_mat = pe.Node(niu.Function(
        input_names=['in_bval', 'in_xfms'], output_names=['out_files'],
        function=recompose_xfm), name='GatherMatrices')
    merge = pe.Node(niu.Function(
        input_names=['in_dwi', 'in_bval', 'in_corrected'],
        output_names=['out_file'], function=recompose_dwi), name='MergeDWIs')

    merged_volumes = pe.Node(niu.Function(
        input_names=['in_file1', 'in_file2'],
        output_names=['out_file'],
        function=merge_volumes_tdim), name='merge_enhanced_ref_dwis')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_xfms']), name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  getb0,        [('in_file', 'in_file')]),
        (inputnode,  pick_dws,     [('in_file', 'in_dwi'),
                                    ('in_bval', 'in_bval')]),
        (flirt, merged_volumes, [('outputnode.out_ref', 'in_file1'),
                                 ('outputnode.out_file', 'in_file2')]),

        (merged_volumes,  merge,        [('out_file', 'in_dwi')]),
        (inputnode, merge, [('in_bval', 'in_bval')]),
        (inputnode,  flirt,        [('in_mask', 'inputnode.ref_mask'),
                                    ('in_xfms', 'inputnode.in_xfms'),
                                    ('in_bval', 'inputnode.in_bval')]),
        (inputnode,  get_mat,      [('in_bval', 'in_bval')]),
        (getb0,      flirt,        [('roi_file', 'inputnode.reference')]),
        (pick_dws,   flirt,        [('out_file', 'inputnode.in_file')]),
        (flirt,      get_mat,      [('outputnode.out_xfms', 'in_xfms')]),
        (flirt,      mult,         [(('outputnode.out_xfms', _xfm_jacobian),
                                     'operand_value')]),
        (flirt,      split,        [('outputnode.out_file', 'in_file')]),
        (split,      mult,         [('out_files', 'in_file')]),
        (mult,       thres,        [('out_file', 'in_file')]),
        (thres,      merge,        [('out_file', 'in_corrected')]),
        (get_mat,    outputnode,   [('out_files', 'out_xfms')]),
        (merge,      outputnode,   [('out_file', 'out_file')])
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


def epi_pipeline(name='susceptibility_distortion_correction_using_t1'):
    """
    This workflow allows to correct for echo-planareinduced susceptibility artifacts without fieldmap
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

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import create_jacobian_determinant_image, change_itk_transform_type, expend_matrix_list, rotate_bvecs, ants_registration_syn_quick, ants_warp_image_multi_transform, ants_combin_transform
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
