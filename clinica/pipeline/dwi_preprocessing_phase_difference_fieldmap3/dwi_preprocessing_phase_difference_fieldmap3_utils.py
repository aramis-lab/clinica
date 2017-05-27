#!/usr/bin/python

import nipype.interfaces.ants as ants
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.utils import add_empty_vol
from nipype.workflows.dmri.fsl.utils import cleanup_edge_pipeline
from nipype.workflows.dmri.fsl.utils import demean_image
from nipype.workflows.dmri.fsl.utils import insert_mat
from nipype.workflows.dmri.fsl.utils import rads2radsec
from nipype.workflows.dmri.fsl.utils import rotate_bvecs
from nipype.workflows.dmri.fsl.utils import siemens2rads
from nipype.workflows.dmri.fsl.utils import vsm2warp
from nipype.workflows.dmri.fsl.utils import enhance


def diffusion_preprocessing_phasediff_fieldmap(in_dwi, in_bvals, in_bvecs, in_fmap_magnitude, in_fmap_phasediff,
                                               delta_te, echospacing, enc_dir, caps_dir,
                                               register_fmap_on_b0, register_b0_on_b0, work_dir, n_procs,
                                               name='Diffusion_preprocessing_phasediff_fmb'):
    import dwi_preprocessing_phase_difference_fieldmap3_utils as utils
    import nipype.interfaces.utility as nutil
    import nipype.pipeline.engine as npe
    from nipype.workflows.dmri.fsl.utils import apply_all_corrections
    import nipype.interfaces.io as nio
    from os.path import join

    pre = utils.prepare_data(register_b0_on_b0=register_b0_on_b0)
    pre.inputs.inputnode.in_file = in_dwi
    pre.inputs.inputnode.in_bvals = in_bvals
    pre.inputs.inputnode.in_bvecs = in_bvecs

    hmc = utils.hmc_pipeline(name='HeadMotionCorrection')
    hmc.inputs.inputnode.ref_num = 0

    sdc = utils.sdc_fmb(name='EPICorrectionWithPhaseDiffFmap',
                        fugue_params=dict(smooth3d=2.0),
                        register_fmap_on_b0=register_fmap_on_b0)
    sdc.inputs.inputnode.in_fmap_magnitude = in_fmap_magnitude
    sdc.inputs.inputnode.in_fmap_phasediff = in_fmap_phasediff
    sdc.inputs.inputnode.delta_te = delta_te
    sdc.inputs.inputnode.echospacing = echospacing
    sdc.inputs.inputnode.enc_dir = enc_dir

    ecc = utils.ecc_pipeline(name='EddyCurrentCorrection')

    unwarp = apply_all_corrections(name='ApplyAllCorrections')

    bias = utils.remove_bias(name='RemoveBias')

    outputnode = npe.Node(nutil.IdentityInterface(
        fields=['out_preprocessed_dwi', 'out_bvecs', 'out_bvals', 'out_b0_mask']),
        name='outputnode')

    datasink = npe.Node(nio.DataSink(), name='datasink')
    # get the participant_id and session_id from the input in_dwi
    participant_id = (in_dwi.split('/')[-1]).split('_')[0]
    session_id = (in_dwi.split('/')[-1]).split('_')[1]
    caps_identifier = participant_id + '_' + session_id
    datasink.inputs.base_directory = join(caps_dir, 'subjects', participant_id, session_id, 'dwi')
    datasink.inputs.substitutions = [
        ('vol0000_warp_maths_thresh_merged_roi_brain_mask.nii.gz', caps_identifier + '_b0Mask.nii.gz'),
        ('vol0000_maths_thresh_merged.nii.gz', caps_identifier + '_dwi.nii.gz'),
        ('bvecs_rotated.bvec', caps_identifier + '_dwi.bvec'),
        ('bvals', caps_identifier + '_dwi.bval'),
        ('merged_files.nii.gz', caps_identifier + '_dwi.nii.gz'),
        ('merged_files.nii.gz', caps_identifier + '_dwi.nii.gz'),
        ('merged_files.nii.gz', caps_identifier + '_dwi.nii.gz'),
        ('vol0000_unwarped_thresh_merged.nii.gz', caps_identifier + '_unwarped_thresh_merged.nii.gz'),
        ('vol0000_unwarped_thresh_merged_concatwarp.nii.gz',
         caps_identifier + '_unwarped_thresh_merged_concatwarp.nii.gz'),
        ('vol0', caps_identifier + '_hmc-dwi-0'),
        ('_flirt.mat', '.mat')
    ]

    wf = npe.Workflow(name=name)
    wf.base_dir = work_dir
    wf.connect([

        # Head-motion correction:
        (pre, hmc, [('outputnode.dwi_b0_merge', 'inputnode.in_file'),
                    ('outputnode.out_bvals', 'inputnode.in_bval'),
                    ('outputnode.out_bvecs', 'inputnode.in_bvec')]),
        (pre, hmc, [('outputnode.mask_b0', 'inputnode.in_mask')]),
        # Susceptibility distortion correction:
        (hmc, sdc, [('outputnode.out_file', 'inputnode.in_file')]),
        (hmc, sdc, [('outputnode.mask_B0', 'inputnode.in_mask')]),
        # Eddy-currents correction:
        (hmc, ecc, [('outputnode.out_xfms', 'inputnode.in_xfms')]),
        (pre, ecc, [('outputnode.out_bvals', 'inputnode.in_bval')]),
        (pre, ecc, [('outputnode.dwi_b0_merge', 'inputnode.in_file')]),
        (pre, ecc, [('outputnode.mask_b0', 'inputnode.in_mask')]),
        # Apply all corrections:
        (pre, unwarp, [('outputnode.dwi_b0_merge', 'inputnode.in_dwi')]),
        (hmc, unwarp, [('outputnode.out_xfms', 'inputnode.in_hmc')]),
        (ecc, unwarp, [('outputnode.out_xfms', 'inputnode.in_ecc')]),
        (sdc, unwarp, [('outputnode.out_warp', 'inputnode.in_sdc')]),
        # Remove bias:
        (unwarp, bias, [('outputnode.out_file', 'inputnode.in_file')]),
        #        (mask_b0,   bias, [('mask_file', 'inputnode.in_mask')]),
        # # Outputnode:
        (hmc, outputnode, [('outputnode.out_bvec', 'out_bvecs')]),
        (pre, outputnode, [('outputnode.out_bvals', 'out_bvals')]),
        (bias, outputnode, [('outputnode.out_file', 'out_preprocessed_dwi')]),
        (bias, outputnode, [('outputnode.b0_mask', 'out_b0_mask')]),
        # Datasink:
        (hmc, datasink, [('outputnode.out_xfms', 'preprocessing.head-motion-correction.@out_matrices')]),
        (pre, datasink, [('outputnode.out_bvals', 'preprocessing.head-motion-correction.@out_bval')]),
        (hmc, datasink, [('outputnode.out_bvec', 'preprocessing.head-motion-correction.@out_bvec')]),
        (hmc, datasink, [('outputnode.out_file', 'preprocessing.head-motion-correction.@out_file')]),
        (ecc, datasink, [('outputnode.out_file', 'preprocessing.eddy-currents-correction.@out_file')]),
        (sdc, datasink,
         [('outputnode.out_file', 'preprocessing.susceptibility-distortion-correction.@out_file')]),
        (
        sdc, datasink, [('outputnode.out_vsm', 'preprocessing.susceptibility-distortion-correction.@out_vsm')]),
        (sdc, datasink,
         [('outputnode.out_warp', 'preprocessing.susceptibility-distortion-correction.@out_warp')]),
        (sdc, datasink, [('outputnode.out_registered_fmap',
                          'preprocessing.susceptibility-distortion-correction.@out_registered_fmap')]),
        (pre, datasink, [('outputnode.out_bvals', 'preprocessing.@out_bvals')]),
        (hmc, datasink, [('outputnode.out_bvec', 'preprocessing.@out_bvecs')]),
        (bias, datasink, [('outputnode.b0_mask', 'preprocessing.@out_b0_mask')]),
        (bias, datasink, [('outputnode.out_file', 'preprocessing.@out_preprocessed_dwi')])
    ])

    return wf.run(plugin='MultiProc', plugin_args={'n_procs': n_procs})
    # return wf

def prepare_data(register_b0_on_b0, name='prepare_data', low_bval=5.0):
    """
    Create a pipeline that prepare the data for further corrections. This pipeline coregister the B0 images and then average it in order
    to obtain only one average B0 images. The bvectors and bvales are update according to the modifications.

    Parameters
    ----------
    num_b0s : INT
      Mandatory input. Number of the B0 volumes in the dataset.

    Inputnode
    ---------
    in_file : FILE
      Mandatory input. Input dwi file.
    in_bvecs : FILE
      Mandatory input. Vector file of the diffusion directions of the dwi dataset.
    in_bvals : FILE
      Mandatory input. B values file.

    Outputnode
    ----------

        outputnode.dwi_b0_merge - average of B0 images merged to the DWIs
        outputnode.b0_reference - average of the B0 images or the only B0 image
        outputnode.out_bvec - updated gradient vectors table
        outputnode.out_bvals - updated gradient values table
        outputnode.mask_b0 - Binary mask obtained from the average of the B0 images

    """
    import clinica.pipeline.dwi.dwi_preprocessing_utils as dwi_utils

    inputnode = pe.Node(interface=niu.IdentityInterface(fields=["in_file", "in_bvecs", "in_bvals"]), name="inputnode")

    b0_dwi_split = pe.Node(niu.Function(input_names=['in_file', 'in_bvals', 'in_bvecs'], output_names=['out_b0', 'out_dwi', 'out_bvals', 'out_bvecs'], function=dwi_utils.b0_dwi_split), name='b0_dwi_split')

    b0_flirt = dwi_utils.b0_flirt_pipeline(name='b0_co_registration')

    b0_avg = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'], function=dwi_utils.b0_average), name='b0_average')

    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True), name='mask_b0')

    insert_b0_into_dwi = pe.Node(niu.Function(input_names=['in_b0', 'in_dwi', 'in_bvals', 'in_bvecs'], output_names=['out_dwi', 'out_bvals', 'out_bvecs'],
                                              function=dwi_utils.insert_b0_into_dwi), name='insert_b0avg_into_dwi')

    outputnode = pe.Node(niu.IdentityInterface(fields=['mask_b0', 'b0_reference','out_bvecs', 'dwi_b0_merge', 'out_bvals'  ]), name='outputnode')

    wf = pe.Workflow(name=name)

    if  register_b0_on_b0:
        wf.connect([
            (inputnode,    b0_dwi_split, [('in_bvals', 'in_bvals'),
                                          ('in_bvecs', 'in_bvecs'),
                                          ('in_file', 'in_file')]),
            (b0_dwi_split, b0_flirt,           [('out_b0', 'inputnode.in_file')]),
            (b0_flirt,     b0_avg,             [('outputnode.out_file', 'in_file')]),
            (b0_avg,       insert_b0_into_dwi, [('out_file', 'in_b0')]),
            (b0_avg,       mask_b0,            [('out_file', 'in_file')]),
            (b0_dwi_split, insert_b0_into_dwi, [('out_dwi','in_dwi'),
                                                ('out_bvals','in_bvals'),
                                                ('out_bvecs','in_bvecs')]),
            (insert_b0_into_dwi, outputnode, [('out_dwi','dwi_b0_merge'),
                                              ('out_bvals','out_bvals'),
                                              ('out_bvecs','out_bvecs')]),
            (mask_b0,            outputnode, [('mask_file','mask_b0')]),
            (b0_avg,             outputnode, [('out_file','b0_reference')])
            ])
    elif register_b0_on_b0 == False:
        wf.connect([
            (inputnode,             b0_dwi_split,       [('in_bvals', 'in_bvals'),
                                                         ('in_bvecs', 'in_bvecs'),
                                                         ('in_file', 'in_file')]),
            (b0_dwi_split,          insert_b0_into_dwi, [('out_b0', 'in_b0'),
                                                         ('out_dwi','in_dwi'),
                                                         ('out_bvals','in_bvals'),
                                                         ('out_bvecs','in_bvecs')]),
            (b0_dwi_split,          mask_b0,            [('out_b0', 'in_file')]),
            (insert_b0_into_dwi,    outputnode,         [('out_dwi','dwi_b0_merge'),
                                                         ('out_bvals','out_bvals'),
                                                         ('out_bvecs','out_bvecs')]),
            (mask_b0,               outputnode,         [('mask_file','mask_b0')]),
            (insert_b0_into_dwi,    outputnode,         [('out_dwi','b0_reference')])
        ])
    else:
        raise()

    return wf



def hmc_pipeline(name='motion_correct'):
    """
    HMC stands for head-motion correction.

    Creates a pipeline that corrects for head motion artifacts in dMRI
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
    from nipype.workflows.data import get_flirt_schedule
    from clinica.pipeline.dwi.dwi_preprocessing_utils import merge_volumes_tdim
    from clinica.pipeline.dwi.dwi_preprocessing_utils import hmc_split
    from clinica.pipeline.dwi.dwi_preprocessing_workflows import dwi_flirt

#    params = dict(dof=6, interp='spline', cost='normmi', cost_func='normmi', bins=50, save_log=True, padding_size=10,
#                  schedule=get_flirt_schedule('hmc'),
#                  searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4], fine_search=1, coarse_search=10 )

    # params = dict(dof=6, interp='spline', cost='normmi', cost_func='normmi', save_log=True,
    #               no_search=True, bgvalue=0, padding_size=10,
    #               schedule=get_flirt_schedule('hmc'),
    #               searchr_x=[-5, 5], searchr_y=[-5, 5], searchr_z=[-25, 25])
    params = dict(dof=6, bgvalue=0, save_log=True, no_search=True,
                  # cost='mutualinfo', cost_func='mutualinfo', bins=64,
                  schedule=get_flirt_schedule('hmc'))

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_bvec', 'in_bval', 'in_mask', 'ref_num']),
        name='inputnode')

    split = pe.Node(niu.Function(function=hmc_split,
                    input_names=['in_file', 'in_bval', 'ref_num'],
                    output_names=['out_ref', 'out_mov', 'out_bval', 'volid']),
                    name='split_ref_moving')

    flirt = dwi_flirt(flirt_param=params)

    insmat = pe.Node(niu.Function(input_names=['inlist', 'volid'],
                     output_names=['out'], function=insert_mat), name='insert_ref_matrix')

    rot_bvec = pe.Node(niu.Function(input_names=['in_bvec', 'in_matrix'],
                       output_names=['out_file'], function=rotate_bvecs),
                       name='Rotate_Bvec')

    merged_volumes = pe.Node(niu.Function(input_names=['in_file1', 'in_file2'], output_names=['out_file'], function=merge_volumes_tdim), name='merge_reference_moving')

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file',
                         'out_bvec', 'out_xfms', 'mask_B0']), name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, split, [('in_file', 'in_file'),
                            ('in_bval', 'in_bval'),
                            ('ref_num', 'ref_num')]),
        (inputnode, flirt, [('in_mask', 'inputnode.ref_mask')]),
        (split,     flirt, [('out_ref', 'inputnode.reference'),
                            ('out_mov', 'inputnode.in_file'),
                            ('out_bval', 'inputnode.in_bval')]),
        (flirt,     insmat, [('outputnode.out_xfms', 'inlist')]),
        (split,     insmat, [('volid', 'volid')]),
        (inputnode, rot_bvec, [('in_bvec', 'in_bvec')]),
        (insmat,    rot_bvec, [('out', 'in_matrix')]),
        (rot_bvec,         outputnode,       [('out_file', 'out_bvec')]),
        (flirt,            merged_volumes,   [('outputnode.out_ref', 'in_file1'),
                                              ('outputnode.out_file', 'in_file2')]),
        (merged_volumes,   outputnode,       [('out_file', 'out_file')]),
        (insmat,           outputnode,       [('out', 'out_xfms')]),
        (flirt,           outputnode,       [('outputnode.out_ref', 'mask_B0')])
    ])
    return wf




def ecc_pipeline(name='eddy_correct'):
    """
    ECC stands for Eddy currents correction.
    Creates a pipeline that corrects for artifacts induced by Eddy currents in
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
    from clinica.pipeline.dwi.dwi_preprocessing_workflows import dwi_flirt
    from nipype.workflows.data import get_flirt_schedule
    from nipype.workflows.dmri.fsl.utils import extract_bval
    from nipype.workflows.dmri.fsl.utils import recompose_xfm
    from nipype.workflows.dmri.fsl.utils import recompose_dwi
    from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian
    from clinica.pipeline.dwi.dwi_preprocessing_utils import merge_volumes_tdim


#    params = dict(dof=12, no_search=True, interp='spline', bgvalue=0,
#                  schedule=get_flirt_schedule('ecc'))

    # params = dict(dof=12, interp='spline', cost='normmi', cost_func='normmi', save_log=True,
    #               no_search=True, bgvalue=0, padding_size=10,
    #               schedule=get_flirt_schedule('ecc'),
    #               searchr_x=[-5, 5], searchr_y=[-5, 5], searchr_z=[-25, 25])

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

    merged_volumes = pe.Node(niu.Function(input_names=['in_file1', 'in_file2'], output_names=['out_file'], function=merge_volumes_tdim), name='merge_enhanced_ref_dwis')

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



def  sdc_fmb(fugue_params=dict(smooth3d=2.0),
            register_fmap_on_b0=True,
            name='fmb_correction'):
    """
    SDC stands for susceptibility distortion correction. FMB stands for
    fieldmap-based.

    The fieldmap based method (FMB) implements SDC by using a mapping of the
    B0 field as proposed by [Jezzard95]_. This workflow uses the implementation
    of FSL (`FUGUE <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE>`_). Phase
    unwrapping is performed using `PRELUDE
    <http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/fugue/prelude.html>`_
    [Jenkinson03]_. Preparation of the fieldmap is performed reproducing the
    script in FSL `fsl_prepare_fieldmap
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data>`_.

    Parameters
    ----------
    in_file : FILE
      Mandatory input. Dwi dataset.
    in_bval : FILE
      Mandatory input. Bval file.
    in_mask : FILE
      Mandatory input. Mask file.
    bmap_mag : FILE
      Mandatory input. Grefield map. Magnitude.
    bmap_pha : FILE
      Mandatory input. Grefield map. Phase.

    Outputs
    ------
    out_file : FILE
      Output.
    out_vsm : FILE
      Output. The set of dwi volumes.
    out_warp : FILE
      Output. The bvalues corresponding to the out_dwi.

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

    echo_spacing = 1/(BandwidthPerPixelPhaseEncode x (AcquisitionMatrixText component #1))
    """

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file',
                                                      'in_mask', 'in_fmap_phasediff', 'in_fmap_magnitude', 'delta_te',
                                                      'echospacing', 'enc_dir']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_vsm',
                         'out_warp', 'out_registered_fmap']),
                         name='outputnode')

    getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')
    firstmag = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='GetFirst')
    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='n4_magnitude')
    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='bet_n4_magnitude')
    dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='dilate_bet')
    pha2rads = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
                       function=siemens2rads), name='PreparePhase')
    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')
    rad2rsec = pe.Node(niu.Function(input_names=['in_file', 'delta_te'],
                       output_names=['out_file'], function=rads2radsec), name='ToRadSec')

    ############ Solution 3: use ants
    fmm2b0 = pe.Node(ants.Registration(output_warped_image=True),
                     name="FMm_to_B0")
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

    applyxfm = pe.Node(ants.ApplyTransforms(
        dimension=3, interpolation='BSpline'), name='FMp_to_B0')


    pre_fugue = pe.Node(fsl.FUGUE(save_fmap=True), name='PreliminaryFugue')
    demean = pe.Node(niu.Function(input_names=['in_file', 'in_mask'],
                     output_names=['out_file'], function=demean_image),
                     name='DemeanFmap')

    cleanup = cleanup_edge_pipeline()

    addvol = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
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
        (inputnode, pha2rads, [('in_fmap_phasediff', 'in_file')]),
        (inputnode, getb0, [('in_file', 'in_file')]),
        (inputnode, firstmag, [('in_fmap_magnitude', 'in_file')]),
        (firstmag, n4, [('roi_file', 'input_image')]),
        (n4, bet, [('output_image', 'in_file')]),
        (bet, dilate, [('mask_file', 'in_file')]),
        (pha2rads, prelude, [('out_file', 'phase_file')]),
        (n4, prelude, [('output_image', 'magnitude_file')]),
        (dilate, prelude, [('out_file', 'mask_file')]),
        (prelude, rad2rsec, [('unwrapped_phase_file', 'in_file')]),
        (inputnode, rad2rsec, [('delta_te', 'delta_te')]),
    ])
    if register_fmap_on_b0:
        wf.connect([
            (getb0, applyxfm, [('roi_file', 'reference_image')]),
            (rad2rsec, applyxfm, [('out_file', 'input_image')]),
            (fmm2b0, applyxfm, [
                ('forward_transforms', 'transforms'),
                ('forward_invert_flags', 'invert_transform_flags')]),
            (applyxfm, pre_fugue, [('output_image', 'fmap_in_file')]),
            ###### Solution using ants
            (getb0, fmm2b0, [('roi_file', 'fixed_image')]),
            (n4, fmm2b0, [('output_image', 'moving_image')]),
            (inputnode, fmm2b0, [('in_mask', 'fixed_image_mask')]),
            (dilate, fmm2b0, [('out_file', 'moving_image_mask')]),
            (inputnode, pre_fugue, [('in_mask', 'mask_file')]),
            (pre_fugue, demean, [('fmap_out_file', 'in_file')]),
            (inputnode, demean, [('in_mask', 'in_mask')]),
            (demean, cleanup, [('out_file', 'inputnode.in_file')]),
            (inputnode, cleanup, [('in_mask', 'inputnode.in_mask')]),
            (cleanup, addvol, [('outputnode.out_file', 'in_file')]),
            (inputnode, vsm, [('in_mask', 'mask_file')]),
            (inputnode, vsm, [('echospacing', 'dwell_time')]),
            (inputnode, vsm, [('delta_te', 'asym_se_time')]),
            (addvol, vsm, [('out_file', 'fmap_in_file')]),
            (inputnode, split, [('in_file', 'in_file')]),
            (split, unwarp, [('out_files', 'in_file')]),
            (vsm, unwarp, [('shift_out_file', 'shift_in_file')]),
            (inputnode, unwarp, [('enc_dir', 'unwarp_direction')]),
            (unwarp, thres, [('unwarped_file', 'in_file')]),
            (thres, merge, [('out_file', 'in_files')]),
            (inputnode, vsm2dfm, [('enc_dir', 'inputnode.enc_dir')]),
            (merge, vsm2dfm, [('merged_file', 'inputnode.in_ref')]),
            (vsm, vsm2dfm, [('shift_out_file', 'inputnode.in_vsm')]),
            (applyxfm, outputnode, [('output_image', 'out_registered_fmap')]),
            (rad2rsec, outputnode, [('out_file', 'out_native_fmap')]),
            (merge, outputnode, [('merged_file', 'out_file')]),
            (vsm, outputnode, [('shift_out_file', 'out_vsm')]),
            (vsm2dfm, outputnode, [('outputnode.out_warp', 'out_warp')])
        ])
    else:
        wf.connect([
            (rad2rsec, pre_fugue, [('out_file', 'fmap_in_file')]),
            (inputnode, pre_fugue, [('in_mask', 'mask_file')]),
            (pre_fugue, demean, [('fmap_out_file', 'in_file')]),
            (inputnode, demean, [('in_mask', 'in_mask')]),
            (demean, cleanup, [('out_file', 'inputnode.in_file')]),
            (inputnode, cleanup, [('in_mask', 'inputnode.in_mask')]),
            (cleanup, addvol, [('outputnode.out_file', 'in_file')]),
            (inputnode, vsm, [('in_mask', 'mask_file')]),
            (inputnode, vsm, [('echospacing', 'echospacing')]),
            (inputnode, vsm, [('delta_te', 'delta_te')]),
            (addvol, vsm, [('out_file', 'fmap_in_file')]),
            (inputnode, split, [('in_file', 'in_file')]),
            (split, unwarp, [('out_files', 'in_file')]),
            (vsm, unwarp, [('shift_out_file', 'shift_in_file')]),
            (inputnode, unwarp, [('enc_dir', 'enc_dir')]),
            (unwarp, thres, [('unwarped_file', 'in_file')]),
            (thres, merge, [('out_file', 'in_files')]),
            (merge, vsm2dfm, [('merged_file', 'inputnode.in_ref')]),
            (vsm, vsm2dfm, [('shift_out_file', 'inputnode.in_vsm')]),
            (inputnode, vsm2dfm, [('enc_dir', 'inputnode.enc_dir')]),
            (rad2rsec, outputnode, [('out_file', 'out_native_fmap')]),
            (merge, outputnode, [('merged_file', 'out_file')]),
            (vsm, outputnode, [('shift_out_file', 'out_vsm')]),
            (vsm2dfm, outputnode, [('outputnode.out_warp', 'out_warp')])
        ])

    return wf


def sdc_fmb_twophase(name='fmb_correction',
            fugue_params=dict(smooth3d=2.0),
            fmap_params=dict(delta_te=2.46e-3),
            epi_params=dict(echospacing=0.39e-3,
                            enc_dir='y')):
    """
    SDC stands for susceptibility distortion correction. FMB stands for
    fieldmap-based.

    The fieldmap based method (FMB) implements SDC by using a mapping of the
    B0 field as proposed by [Jezzard95]_. This workflow uses the implementation
    of FSL (`FUGUE <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE>`_). Phase
    unwrapping is performed using `PRELUDE
    <http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/fugue/prelude.html>`_
    [Jenkinson03]_. Preparation of the fieldmap is performed reproducing the
    script in FSL `fsl_prepare_fieldmap
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#SIEMENS_data>`_.

    Parameters
    ----------
    in_file : FILE
      Mandatory input. Dwi dataset.
    in_bval : FILE
      Mandatory input. Bval file.
    in_mask : FILE
      Mandatory input. Mask file.
    bmap_mag : FILE
      Mandatory input. Grefield map. Magnitude.
    bmap_pha : FILE
      Mandatory input. Grefield map. Phase.

    Outputs
    ------
    out_file : FILE
      Output.
    out_vsm : FILE
      Output. The set of dwi volumes.
    out_warp : FILE
      Output. The bvalues corresponding to the out_dwi.

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

    echo_spacing = 1/(BandwidthPerPixelPhaseEncode x (AcquisitionMatrixText component #1))
    """
    from clinica.pipeline.dwi.dwi_preprocessing_utils import convert_phase_in_radians
    from clinica.pipeline.dwi.dwi_preprocessing_utils import create_phase_in_radsec

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_mask', 'in_fmap_phase1', 'in_fmap_phase2', 'in_fmap_magnitude1', 'in_fmap_magnitude2']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_vsm', 'out_warp', 'out_native_fmap']),
        name='outputnode')

    getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')
#    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='n4_magnitude')
    n4_1 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='n4_magnitude1')
    n4_2 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='n4_magnitude2')

#    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='bet_n4_magnitude')
    bet_1 = pe.Node(fsl.BET(frac=0.4, mask=True), name='bet_n4_magnitude1')
    bet_2 = pe.Node(fsl.BET(frac=0.4, mask=True), name='bet_n4_magnitude2')
#    dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
#                     args='-kernel sphere 5 -dilM'), name='dilate_bet')
    dilate_1 = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='dilate_bet_1')
    dilate_2 = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='dilate_bet_2')

#    phase1_in_rad = pe.Node(niu.Function(input_names=['in_file', 'name_output_file'], output_names=['out_file'],
#                       function=convert_phase_in_radians), name='Phase1InRad')
#    phase2_in_rad = pe.Node(niu.Function(input_names=['in_file', 'name_output_file'], output_names=['out_file'],
#                       function=convert_phase_in_radians), name='Phase2InRad')
    phase1_in_rad = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
                       function=siemens2rads), name='Phase1InRad')
    phase2_in_rad = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
                       function=siemens2rads), name='Phase2InRad')

#    phase1_unwarp = pe.Node(fsl.PRELUDE(process3d=True), name='Phase1Unwarp')
#    phase2_unwarp = pe.Node(fsl.PRELUDE(process3d=True), name='Phase2Unwarp')

    phase_in_rsec =pe.Node(niu.Function(input_names=['in_phase1', 'in_phase2', 'delta_te', 'out_file'], output_names=['out_file'],
                       function=create_phase_in_radsec), name='PhaseInRadSec')
    phase_in_rsec.inputs.delta_te = fmap_params['delta_te']

#    pha2rads = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
#                       function=siemens2rads), name='PreparePhase')
#    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')
#    rad2rsec = pe.Node(niu.Function(input_names=['in_file', 'delta_te'],
#                       output_names=['out_file'], function=rads2radsec), name='ToRadSec')
#    rad2rsec.inputs.delta_te = fmap_params['delta_te']

    flirt = pe.Node(fsl.FLIRT(interp='spline', cost='normmi', cost_func='normmi',
                    dof=6, bins=64, save_log=True, padding_size=10,
                    searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4],
                    fine_search=1, coarse_search=10),
                    name='BmapMag2B0')
    applyxfm = pe.Node(fsl.ApplyXfm(interp='spline', padding_size=10, apply_xfm=True),
                       name='BmapPha2B0')

    pre_fugue = pe.Node(fsl.FUGUE(save_fmap=True), name='PreliminaryFugue')
    demean = pe.Node(niu.Function(input_names=['in_file', 'in_mask'],
                     output_names=['out_file'], function=demean_image),
                     name='DemeanFmap')

    cleanup = cleanup_edge_pipeline()

    addvol = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
                     function=add_empty_vol), name='AddEmptyVol')

    vsm = pe.Node(fsl.FUGUE(save_shift=True, **fugue_params),
                  name="ComputeVSM")
    vsm.inputs.asym_se_time = fmap_params['delta_te']
    vsm.inputs.dwell_time = epi_params['echospacing']

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
    unwarp = pe.MapNode(fsl.FUGUE(icorr=True, forward_warping=False),
                        iterfield=['in_file'], name='UnwarpDWIs')
    unwarp.inputs.unwarp_direction = epi_params['enc_dir']
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')
    vsm2dfm = vsm2warp()
    vsm2dfm.inputs.inputnode.scaling = 1.0
    vsm2dfm.inputs.inputnode.enc_dir = epi_params['enc_dir']

    wf = pe.Workflow(name=name)
    wf.connect([
#        (inputnode, pha2rads, [('bmap_pha', 'in_file')]),
        (inputnode, phase1_in_rad, [('in_fmap_phase1', 'in_file')]),
        (inputnode, phase2_in_rad, [('in_fmap_phase2', 'in_file')]),
        (inputnode, getb0, [('in_file', 'in_file')]),
#        (inputnode, n4, [('bmap_mag', 'input_image')]),
        (inputnode, n4_1, [('in_fmap_magnitude1', 'input_image')]),
        (inputnode, n4_2, [('in_fmap_magnitude2', 'input_image')]),
#        (n4, bet, [('output_image', 'in_file')]),
        (n4_1, bet_1, [('output_image', 'in_file')]),
        (n4_2, bet_2, [('output_image', 'in_file')]),
#        (bet, dilate, [('mask_file', 'in_file')]),
        (bet_1, dilate_1, [('mask_file', 'in_file')]),
        (bet_2, dilate_2, [('mask_file', 'in_file')]),
        #        (pha2rads, prelude, [('out_file', 'phase_file')]),
#        (n4,       prelude, [('output_image', 'magnitude_file')]),
#        (dilate,   prelude, [('out_file', 'mask_file')]),
#        (phase1_in_rad, phase1_unwarp, [('out_file', 'phase_file')]),
#        (n4_1,          phase1_unwarp, [('output_image', 'magnitude_file')]),
#        (dilate_1,      phase1_unwarp, [('out_file', 'mask_file')]),
#        (phase2_in_rad, phase2_unwarp, [('out_file', 'phase_file')]),
#        (n4_2,          phase2_unwarp, [('output_image', 'magnitude_file')]),
#        (dilate_2,      phase2_unwarp, [('out_file', 'mask_file')]),
#        (prelude, rad2rsec, [('unwrapped_phase_file', 'in_file')]),
        (phase1_in_rad, phase_in_rsec, [('out_file', 'in_phase1')]),
        (phase2_in_rad, phase_in_rsec, [('out_file', 'in_phase2')]),
#        (phase1_unwarp, phase_in_rsec, [('unwrapped_phase_file', 'in_phase1')]),
#        (phase2_unwarp, phase_in_rsec, [('unwrapped_phase_file', 'in_phase2')]),

        (getb0,     flirt, [('roi_file', 'reference')]),
        (inputnode, flirt, [('in_mask', 'ref_weight')]),
#        (n4,        flirt, [('output_image', 'in_file')]),
#        (dilate,    flirt, [('out_file', 'in_weight')]),
        (n4_1,      flirt, [('output_image', 'in_file')]),
        (dilate_1,  flirt, [('out_file', 'in_weight')]),
        (getb0,    applyxfm, [('roi_file', 'reference')]),
#        (rad2rsec, applyxfm, [('out_file', 'in_file')]),
        (phase_in_rsec, applyxfm, [('out_file', 'in_file')]),
        (flirt,    applyxfm, [('out_matrix_file', 'in_matrix_file')]),
        (applyxfm,  pre_fugue, [('out_file', 'fmap_in_file')]),
        (inputnode, pre_fugue, [('in_mask', 'mask_file')]),
        (pre_fugue, demean, [('fmap_out_file', 'in_file')]),
        (inputnode, demean, [('in_mask', 'in_mask')]),
        (demean,    cleanup, [('out_file', 'inputnode.in_file')]),
        (inputnode, cleanup, [('in_mask', 'inputnode.in_mask')]),
        (cleanup, addvol, [('outputnode.out_file', 'in_file')]),
        (inputnode, vsm, [('in_mask', 'mask_file')]),
        (addvol,    vsm, [('out_file', 'fmap_in_file')]),
        (inputnode, split, [('in_file', 'in_file')]),
        (split, unwarp, [('out_files', 'in_file')]),
        (vsm,   unwarp, [('shift_out_file', 'shift_in_file')]),
        (unwarp, thres, [('unwarped_file', 'in_file')]),
        (thres,  merge, [('out_file', 'in_files')]),
        (merge, vsm2dfm, [('merged_file', 'inputnode.in_ref')]),
        (vsm,   vsm2dfm, [('shift_out_file', 'inputnode.in_vsm')]),
        (applyxfm, outputnode, [('out_file', 'out_registered_fmap')]),
        (merge,   outputnode, [('merged_file', 'out_file')]),
        (vsm,     outputnode, [('shift_out_file', 'out_vsm')]),
        (vsm2dfm, outputnode, [('outputnode.out_warp', 'out_warp')])
    ])
    return wf


def sdc_syb_pipeline(name='sdc_syb_correct'):

    """
    SDC stands for susceptibility distortion correction and SYB stand for SyN based. This workflow
    allows to correct for echo-planare induced susceptibility artifacts without fieldmap
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

    outputnode.out_dwi - corrected dwi file
    outputnode.out_bvec - rotated gradient vectors table
    outputnode.out_b0_to_t1_rigid_body_matrix - B0 to T1 image FLIRT rigid body fsl coregistration matrix
    outputnode.out_t1_coregistered_to_b0 - T1 image rigid body coregistered to the B0 image
    outputnode.out_b0_to_t1_affine_matrix - B0 to T1 image ANTs affine itk coregistration matrix
    outputnode.out_b0_to_t1_syn_defomation_field - B0 to T1 image ANTs SyN itk warp
    outputnode.out_warp - Out warp allowing DWI to T1 registration and susceptibilty induced artifacts correction

    Example
    -------
    >>> epi = epi_pipeline()
    >>> epi.inputs.inputnode.in_dwi = 'DWI.nii'
    >>> epi.inputs.inputnode.in_t1 = 'T1.nii'
    >>> epi.run() # doctest: +SKIP
    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl
    from clinica.pipeline.dwi.dwi_registration import ants_registration_syn_quick, antscombintransform

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_t1', 'in_dwi']), name='inputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_to_t1 = pe.Node(interface=fsl.FLIRT(dof=6), name = 'flirt_b0_to_t1')
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

    ants_registration_syn_quick = pe.Node(interface=niu.Function(
        input_names=['fixe_image', 'moving_image'],
        output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
        function=ants_registration_syn_quick), name='ants_registration_syn_quick')

    merge_transform = pe.Node(niu.Merge(2), name='MergeTransforms')

    combin_warp = pe.Node(interface=niu.Function(
        input_names=['in_file', 'transforms_list', 'reference'],
        output_names=['out_warp'],
        function=antscombintransform), name='combin_warp')

    coeffs = pe.Node(fsl.WarpUtils(out_format='spline'), name='CoeffComp')

    fsl_transf = pe.Node(fsl.WarpUtils(out_format='field'), name='fsl_transf')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(), iterfield=['in_file'],name='apply_warp')
    apply_warp.inputs.interp = 'spline'

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_b0_to_t1_rigid_body_matrix', 'out_t1_to_b0_rigid_body_matrix', 'out_t1_coregistered_to_b0',
                'out_b0_to_t1_syn_defomation_field', 'out_b0_to_t1_affine_matrix', 'out_dwi', 'out_warp']),
        name='outputnode')

    wf = pe.Workflow(name='sdc_syb_pipeline')
    wf.connect([
        (inputnode, split, [('in_dwi', 'in_file')]),
        (split, pick_ref, [('out_files', 'inlist')]),
        (pick_ref, flirt_b0_to_t1, [('out', 'in_file')]),
        (inputnode, flirt_b0_to_t1, [('in_t1', 'reference')]),
        (flirt_b0_to_t1, invert_xfm, [('out_matrix_file', 'in_file')]),
        (invert_xfm, apply_xfm, [('out_file', 'in_matrix_file')]),
        (inputnode, apply_xfm, [('in_t1', 'in_file')]),
        (pick_ref, apply_xfm, [('out', 'reference')]),
        (apply_xfm, ants_registration_syn_quick, [('out_file', 'fixe_image')]),
        (pick_ref, ants_registration_syn_quick, [('out', 'moving_image')]),
        (ants_registration_syn_quick, merge_transform, [('affine_matrix', 'in2'),
                                                     ('warp', 'in1')]),
        (pick_ref, combin_warp, [('out', 'in_file')]),
        (merge_transform, combin_warp, [('out', 'transforms_list')]),
        (apply_xfm, combin_warp, [('out_file', 'reference')]),
        (apply_xfm, coeffs, [('out_file', 'reference')]),
        (combin_warp, coeffs, [('out_warp', 'in_file')]),
        (coeffs, fsl_transf, [('out_file', 'in_file')]),
        (apply_xfm, fsl_transf, [('out_file', 'reference')]),
        (fsl_transf, apply_warp, [('out_file', 'field_file')]),
        (split, apply_warp, [('out_files', 'in_file')]),
        (apply_xfm, apply_warp, [('out_file', 'ref_file')]),
        (apply_warp, thres, [('out_file', 'in_file')]),
        (thres, merge, [('out_file', 'in_files')]),
        (merge,                    outputnode, [('merged_file', 'out_dwi')]),
        (flirt_b0_to_t1,            outputnode, [('out_matrix_file', 'out_b0_to_t1_rigid_body_matrix')]),
        (invert_xfm,               outputnode, [('out_file', 'out_t1_to_b0_rigid_body_matrix')]),
        (apply_xfm,                outputnode, [('out_file', 'out_t1_coregistered_to_b0')]),
        (ants_registration_syn_quick, outputnode, [('warp', 'out_b0_to_t1_syn_defomation_field'),
                                                ('affine_matrix', 'out_b0_to_t1_affine_matrix')]),
        (fsl_transf,               outputnode, [('out_file', 'out_warp')])
    ])

    return wf



def apply_all_corrections_syb(name='UnwarpArtifacts'):
    """
    Combines two lists of linear transforms with the deformation field
    map obtained epi_correction by Ants.
    Additionally, computes the corresponding bspline coefficients and
    the map of determinants of the jacobian.
    """

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_sdc_syb', 'in_hmc','in_ecc', 'in_dwi', 'in_t1']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_warp', 'out_coeff', 'out_jacobian']),
        name='outputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')

    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_to_t1 = pe.Node(interface=fsl.FLIRT(dof=6), name = 'flirt_b0_to_t1')
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

    concat_hmc_ecc = pe.MapNode(fsl.ConvertXFM(), name="concat_hmc_ecc", iterfield=['in_file', 'in_file2'])
    concat_hmc_ecc.inputs.concat_xfm = True

    warps = pe.MapNode(fsl.ConvertWarp(), iterfield=['premat'], name='ConvertWarp')

    unwarp = pe.MapNode(interface=fsl.ApplyWarp(), iterfield=['in_file', 'field_file'],name='unwarp_warp')
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

#    datasink = pe.Node(nio.DataSink(), name='datasink')
#    datasink.inputs.base_directory = op.join(datasink_directory, 'apply_all_correction/')

    wf = pe.Workflow(name=name)
    wf.connect([(inputnode, concat_hmc_ecc, [('in_ecc','in_file2')]),
                (inputnode, concat_hmc_ecc, [('in_hmc','in_file')]),
                (concat_hmc_ecc, warps, [('out_file','premat')]),
                (inputnode, warps, [('in_sdc_syb','warp1')]),
                (inputnode, split, [('in_dwi', 'in_file')]),
                (split, pick_ref, [('out_files','inlist')]),
                (pick_ref, flirt_b0_to_t1, [('out','in_file')]),
                (inputnode, flirt_b0_to_t1, [('in_t1','reference')]),
                (flirt_b0_to_t1, invert_xfm, [('out_matrix_file','in_file')]),
                (invert_xfm, apply_xfm, [('out_file','in_matrix_file')]),
                (inputnode, apply_xfm, [('in_t1','in_file')]),
                (pick_ref, apply_xfm, [('out','reference')]),
                (apply_xfm, warps, [('out_file','reference')]),
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

#    wf.connect([(warps, datasink, [('out_file', 'out_warp')]),
#                (coeffs, datasink, [('out_file', 'out_coeff')]),
#                (jacobian, datasink, [('out_jacobian', 'out_jacobian')]),
#                (merge, datasink, [('merged_file', 'out_file')])])

    return wf




def dwi_flirt(name='DWICoregistration', excl_nodiff=False,
              flirt_param={}):
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
    from nipype.workflows.dmri.fsl.utils import _checkinitxfm

    inputnode = pe.Node(niu.IdentityInterface(fields=['reference',
                        'in_file', 'ref_mask', 'in_xfms', 'in_bval']),
                        name='inputnode')

    initmat = pe.Node(niu.Function(input_names=['in_bval', 'in_xfms',
                      'excl_nodiff'], output_names=['init_xfms'],
                                   function=_checkinitxfm), name='InitXforms')
    initmat.inputs.excl_nodiff = excl_nodiff
    dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='MskDilate')
    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='Bias')
    flirt = pe.MapNode(fsl.FLIRT(**flirt_param), name='CoRegistration',
                       iterfield=['in_file', 'in_matrix_file'])
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')
    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file',
                         'out_xfms', 'out_ref']), name='outputnode')
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
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file','b0_mask']),
                         name='outputnode')

    getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')

    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True), name='mask_b0')

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
        (inputnode, getb0, [('in_file', 'in_file')]),
        (getb0,   n4, [('roi_file', 'input_image')]),
        (getb0, mask_b0, [('roi_file', 'in_file')]),
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

####################################################
####################### Function
####################################################


def merge_volumes_tdim(in_file1, in_file2):
    """
    Merge 'in_file1' and 'in_file2' in the t dimension.

    Parameters
    ----------
    in_file1 : FILE
      Mandatory input. First volume.
    in_file2 : FILE
      Mandatory input. Second volume.

    Output
    ------
    out_file : FILE
      Output. The two volumes merged.

    """
    import os.path as op
    import os

    out_file = op.abspath('merged_files.nii.gz')
    cmd = 'fslmerge -t '+ out_file + ' ' + in_file1 + ' ' + in_file2
    os.system(cmd)
    return out_file

def b0_dwi_split(in_file, in_bvals, in_bvecs, lowbval=5.0):
    """
    Split the volumes into two datasets :
     - the first dataset contains the set of B0 volumes.
     - the second dataset contains the set of dwi volumes.

    Parameters
    ----------
    in_file : FILE
      Mandatory input. Dwi dataset.
    in_bvals : FILE
      Mandatory input. Bval file.
    in_bvecs : FILE
      Mandatory input. Bvecs file.
    lowbval : FLOAT
      Optional input. Define the B0 volumes as all volume bval <= lowbval. Default lowbval=5.0

    Outputs
    ------
    out_b0 : FILE
      Output. The set of b0 volumes.
    out_dwi : FILE
      Output. The set of dwi volumes.
    out_bvals : FILE
      Output. The bvalues corresponding to the out_dwi.
    out_bvecs : FILE
      Output : The bvecs corresponding to the out_dwi.

    """
    import numpy as np
    import nibabel as nib
    import os.path as op
    import warnings

    assert(op.isfile(in_file))
    assert(op.isfile(in_bvals))
    assert(op.isfile(in_bvecs))
    assert(lowbval >= 0)

    im = nib.load(in_file)
    data = im.get_data()
    hdr = im.get_header().copy()
    bvals = np.loadtxt(in_bvals)
    bvecs = np.loadtxt(in_bvecs)

    if bvals.shape[0] == bvecs.shape[0]:
        warnings.warn("Warning - The b-vectors file should be column-wise: the b-vectors will be transposed", UserWarning)
        bvecs = bvecs.T

    lowbs = np.where(bvals <= lowbval)[0]
    out_b0 = op.abspath('b0.nii.gz')
    b0data = data[..., lowbs]
    hdr.set_data_shape(b0data.shape)
    nib.Nifti1Image(b0data, im.get_affine(), hdr).to_filename(out_b0)

    dwi_bvals = np.where(bvals > lowbval)[0]
    out_dwi = op.abspath('dwi.nii.gz')
    dwidata = data[..., dwi_bvals]
    hdr.set_data_shape(dwidata.shape)
    nib.Nifti1Image(dwidata, im.get_affine(), hdr).to_filename(out_dwi)

    bvals_dwi = bvals[dwi_bvals]
    out_bvals = op.abspath('bvals')
    np.savetxt(out_bvals, bvals_dwi, fmt='%d', delimiter=' ')

    bvecs_dwi = np.array([bvecs[0][dwi_bvals].tolist(), bvecs[1][dwi_bvals].tolist(), bvecs[2][dwi_bvals].tolist()])
    out_bvecs = op.abspath('bvecs')
    np.savetxt(out_bvecs, bvecs_dwi, fmt='%10.5f', delimiter=' ')

    return out_b0, out_dwi, out_bvals, out_bvecs



def b0_flirt_pipeline(name='b0coregistration', nb_b0= 10, excl_nodiff=False):
    """
    Rigid registration of the B0 dataset onto the first volume. Rigid
    registration is achieved using FLIRT and the normalized
    correlation.

    Parameters
    ----------
    nb_b0 : INT
      Mandatory input. Number of the B0 volumes in the dataset.

    Inputnode
    ---------
    in_file : FILE
      Mandatory input. B0 dataset.

    Outputnode
    ----------
    out_b0_reg : FILE
      Output. The set of B0 volumes registered to the first volume.
    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces import fsl
    import nipype.interfaces.utility as niu

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']),
                        name='inputnode')
    fslroi_ref = pe.Node(fsl.ExtractROI(args='0 1'), name='b0_reference')
    tsize = nb_b0 - 1
    fslroi_moving = pe.Node(fsl.ExtractROI(args='1 '+str(tsize)), name='b0_moving')
    split_moving = pe.Node(fsl.Split(dimension='t'), name='split_b0_moving')

    bet_ref = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                       name='bet_ref')

    dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='mask_dilate')

    flirt = pe.MapNode(fsl.FLIRT(interp='spline', dof=6, bins=50, save_log=True, cost='corratio', cost_func='corratio', padding_size=10, searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4], fine_search=1, coarse_search=10), name='b0_co_registration', iterfield=['in_file'])

    merge = pe.Node(fsl.Merge(dimension='t'), name='merge_registered_b0s')
    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='remove_negative')
    insert_ref = pe.Node(niu.Function(input_names=['in_file1', 'in_file2'], output_names=['out_file'], function=merge_volumes_tdim), name='concat_ref_moving')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_xfms']), name='outputnode')

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


def count_b0s(in_bvals, low_bval=5.0):
    """
    Count the number of B0 volumes.

    Parameters
    ----------
    in_file : FILE
      Mandatory input. Bvals file.

    Outputs
    -------
    num_b0s : Number of the B0 volumes.
    """
    import numpy as np

    # Counting the number of b0s
    bvals = np.loadtxt(in_bvals)
    num_b0s = len(np.where(bvals <= low_bval)[0])

    return num_b0s


def b0_average(in_file, out_file=None):
    """
    Average the B0 volumes.
    Warning: the B0 volumes must be registered.

    Parameters
    ----------
    in_file : FILE
      Mandatory input. B0 dataset registered.

    Outputs
    -------
    out_file : FILE
      Output. The mean of the B0 volumes.
    """
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("%s_avg_b0%s" % (fname, ext))

    imgs = np.array(nb.four_to_three(nb.load(in_file)))
    b0s = [im.get_data().astype(np.float32)
           for im in imgs]
    b0 = np.average(np.array(b0s), axis=0)

    hdr = imgs[0].get_header().copy()
    hdr.set_data_shape(b0.shape)
    hdr.set_xyzt_units('mm')
    hdr.set_data_dtype(np.float32)
    nb.Nifti1Image(b0, imgs[0].get_affine(), hdr).to_filename(out_file)

    return out_file



def insert_b0_into_dwi(in_b0, in_dwi, in_bvals, in_bvecs):
    """
    This function inserts a b0 volume into the dwi dataset as the
    first volume and updates the bvals and bvecs files.

    Parameters
    ----------
    in_b0 : FILE
      Mandatory input. One B0 volume (could be the average of a b0 dataset).
    in_dwi : FILE
      Mandatory input. Dwi dataset.
    in_bvals : FILE
      Mandatory input. File describing the bvalues of the dwi dataset.
    in_bvecs : FILE
      Mandatory input. File describing the directions of the dwi dataset.

    Outputs
    -------
    out_dwi : FILE
      Output. Diffusion dataset : b0 volume + dwi volumes.
    out_bvals : FILE
      Output. B values update.
    out_bvecs. Directions of diffusion update.
    """
    from clinica.pipeline.dwi.dwi_preprocessing_utils import merge_volumes_tdim
    import os.path as op
    import numpy as np

    assert(op.isfile(in_b0))
    assert(op.isfile(in_dwi))
    assert(op.isfile(in_bvals))
    assert(op.isfile(in_bvecs))

    out_dwi = merge_volumes_tdim(in_b0, in_dwi)

    lst = np.loadtxt(in_bvals).tolist()
    lst.insert(0,0)
    out_bvals = op.abspath('bvals')
    np.savetxt(out_bvals, np.matrix(lst), fmt='%d', delimiter=' ')

    bvecs = np.loadtxt(in_bvecs)
    bvecs_0 = bvecs[0].tolist()
    bvecs_0.insert(0,0.0)
    bvecs_1 = bvecs[1].tolist()
    bvecs_1.insert(0,0.0)
    bvecs_2 = bvecs[2].tolist()
    bvecs_2.insert(0,0.0)
    bvecs_dwi = np.array([bvecs_0, bvecs_1, bvecs_2])
    out_bvecs = op.abspath('bvecs')
    np.savetxt(out_bvecs, bvecs_dwi, fmt='%10.5f', delimiter=' ')

    return out_dwi, out_bvals, out_bvecs



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
    bvecs = np.loadtxt(in_bvec).T
    new_bvecs = []

    if len(bvecs) != len(in_matrix):
        raise RuntimeError(('Number of b-vectors (%d) and rotation '
                            'matrices (%d) should match.') % (len(bvecs),
                                                              len(in_matrix)))

    for bvec, mat in zip(bvecs, in_matrix):
        if np.all(bvec == 0.0):
            new_bvecs.append(bvec)
        else:
            invrot = np.linalg.inv(np.loadtxt(mat))[:3, :3]
            newbvec = invrot.dot(bvec)
            new_bvecs.append((newbvec / np.linalg.norm(newbvec)))

    np.savetxt(out_file, np.array(new_bvecs).T, fmt='%0.15f')
    return out_file





def hmc_split(in_file, in_bval, ref_num=0, lowbval=5.0):
    """
    Selects the reference ('out_ref') and moving ('out_mov') volumes
    from a dwi dataset for the purpose of head motion correction (HMC).

    Parameters
    ----------
    in_file : FILE
      Mandatory input. Dwi dataset.
    in_bval : FILE
      Mandatory input. Bval file.
    ref_num : INT
      Optional input. The reference volume in the dwi dataset. Default ref_num= 0.
    lowbval : FLOAT
      Optional input. Define the volumes with low bval. All volume bval <= lowbval. Default lowbval=5.0

    Outputs
    ------
    out_ref : FILE
      Output. The reference volume.
    out_mov : FILE
     Output. The moving volume to align to the reference volume.
    out_bval : FILE
     Output. The bvalues corresonding to the moving volume.
    volid : INT
      Index of the reference volume.
    """
    import numpy as np
    import nibabel as nib
    import os.path as op

    im = nib.load(in_file)
    data = im.get_data()
    hdr = im.get_header().copy()
    bval = np.loadtxt(in_bval)

    lowbs = np.where(bval <= lowbval)[0]
    assert(ref_num in lowbs)
    volid = ref_num

    out_ref = op.abspath('hmc_ref.nii.gz')
    refdata = data[..., volid]
    hdr.set_data_shape(refdata.shape)
    nib.Nifti1Image(refdata, im.get_affine(), hdr).to_filename(out_ref)

    if volid == 0:
        data = data[..., 1:]
        bval = bval[1:]
    elif volid == (data.shape[-1] - 1):
        data = data[..., :-1]
        bval = bval[:-1]
    else:
        data = np.concatenate((data[..., :volid], data[..., (volid + 1):]),
                              axis=3)
        bval = np.hstack((bval[:volid], bval[(volid + 1):]))

    out_mov = op.abspath('hmc_mov.nii.gz')
    out_bval = op.abspath('bval_split.txt')

    hdr.set_data_shape(data.shape)
    nib.Nifti1Image(data, im.get_affine(), hdr).to_filename(out_mov)
    np.savetxt(out_bval, bval)
    return out_ref, out_mov, out_bval, volid


def convert_phase_in_radians(in_file, out_file=None):
    """
    Convert phase image in radians.
    """
    import math
    import nibabel as nb
    import numpy as np
    import os.path as op
    import os

    assert(op.isfile(in_file))

    img = nb.load(in_file)
    imin = np.amin(img.get_data())
    imax = np.amax(img.get_data())

    if out_file is None:
        out_file = op.abspath('phase_in_rad.nii.gz')


    cmd = 'fslmaths '+in_file+ ' -mul '+str(2.0*math.pi)+" -div "+str(imax)+" "+out_file+' -odt float'
    os.system(cmd)

    data = (img2.get_data().astype(np.float32)-img1.get_data().astype(np.float32)) * (1.0/delta_te)
    nb.Nifti1Image(data, img.get_affine(), img.get_header()).to_filename(out_file)

    return out_file



def create_phase_in_radsec(in_phase1, in_phase2, delta_te, out_file=None):
    """
    Converts input (unwarpped) phase1 and phase2 map to into a fieldmap in rads. (delta_te should be in seconds)
    """
    import numpy as np
    import nibabel as nb
    import os.path as op

    if out_file is None:
        out_file = op.abspath('fmap_radsec.nii.gz')

    img1 = nb.load(in_phase1)
    img2 = nb.load(in_phase2)
    data = (img2.get_data().astype(np.float32)-img1.get_data().astype(np.float32)) * (1.0/delta_te)
    nb.Nifti1Image(data, img1.get_affine(), img1.get_header()).to_filename(out_file)
    return out_file

def distribute_subjects(in_dwi, in_bval, in_bvec, in_fmap_magnitude, in_fmap_phasediff):
    return in_dwi, in_bval, in_bvec, in_fmap_magnitude, in_fmap_phasediff