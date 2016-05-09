#!/usr/bin/python

import os
import os.path as op
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.ants as ants
from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian
from nipype.workflows.dmri.fsl.utils import insert_mat
from nipype.workflows.dmri.fsl.utils import recompose_xfm
from nipype.workflows.dmri.fsl.utils import recompose_dwi
from nipype.workflows.dmri.fsl.utils import extract_bval
from nipype.workflows.dmri.fsl.utils import rotate_bvecs
from nipype.workflows.dmri.fsl.utils import siemens2rads
from nipype.workflows.dmri.fsl.utils import rads2radsec
from nipype.workflows.dmri.fsl.utils import demean_image
from nipype.workflows.dmri.fsl.utils import cleanup_edge_pipeline
from nipype.workflows.dmri.fsl.utils import add_empty_vol
from nipype.workflows.dmri.fsl.utils import vsm2warp
import clinica.pipeline.preprocessing.dwi_utils as predifutils

def prepare_data(datasink_directory, name='prepare_data'):
    """
    Create a pipeline that prepare the data for further corrections. This pipeline coregister the B0 images and then average it in order
    to obtain only one average B0 images. The bvectors and bvales are update according to the modifications.

    Inputnode
    ---------
    dwi_image : FILE
      Mandatory input. Input dwi file.
    bvectors_directions : FILE
      Mandatory input. Vector file of the diffusion directions of the dwi dataset.
    bvalues : FILE
      Mandatory input. B values file.

    Outputnode
    ----------

        outputnode.dwi_b0_merge - average of B0 images merged to the DWIs
        outputnode.b0_average - average of the B0 images
        outputnode.out_bvec - updated gradient vectors table
        outputnode.out_bvals - updated gradient values table
        outputnode.mask_b0 - Binary mask obtained from the average of the B0 images

    """
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=["dwi_image", "bvectors_directions", "bvalues"]), name="inputnode")

    b0_dwi_split = pe.Node(niu.Function(input_names=['in_file', 'in_bvals', 'in_bvecs'], output_names=['out_b0', 'out_dwi', 'out_bvals', 'out_bvecs'],
                                        function=predifutils.b0_dwi_split), name='b0_dwi_split')

    b0_flirt = predifutils.b0_flirt_pipeline(name='b0_co_registration')

    b0_avg = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'], function=predifutils.b0_average), name='b0_average')

    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True), name='mask_b0')

    insert_b0_into_dwi = pe.Node(niu.Function(input_names=['in_b0', 'in_dwi', 'in_bvals', 'in_bvecs'], output_names=['out_dwi', 'out_bvals', 'out_bvecs'],
                                              function=predifutils.insert_b0_into_dwi), name='insert_b0avg_into_dwi')

    datasink = pe.Node(nio.DataSink(), name='datasink_prep')
    datasink.inputs.base_directory = op.join(datasink_directory, 'pre_preprocess/')

    outputnode = pe.Node(niu.IdentityInterface(fields=['mask_b0', 'b0_average','out_bvecs', 'dwi_b0_merge', 'out_bvals'  ]), name='outputnode')


    wf = pe.Workflow(name=name)

    wf.connect([(inputnode, b0_dwi_split, [('bvalues', 'in_bvals'),('bvectors_directions', 'in_bvecs'), ('dwi_image', 'in_file')])])
    wf.connect([(b0_dwi_split, b0_flirt, [('out_b0', 'inputnode.in_file')])])
    wf.connect([(b0_flirt, b0_avg, [('outputnode.out_file', 'in_file')])])
    wf.connect([(b0_avg, insert_b0_into_dwi, [('out_file', 'in_b0')])])
    wf.connect([(b0_avg, mask_b0,[('out_file', 'in_file')])])
    wf.connect([(b0_dwi_split, insert_b0_into_dwi, [('out_dwi','in_dwi'), ('out_bvals','in_bvals'), ('out_bvecs','in_bvecs')])])
    wf.connect([(insert_b0_into_dwi, outputnode, [('out_dwi','dwi_b0_merge'), ('out_bvals','out_bvals'),('out_bvecs','out_bvecs')])])
    wf.connect([(mask_b0, outputnode, [('mask_file','mask_b0')])])
    wf.connect([(b0_avg, outputnode, [('out_file','b0_average')])])
    wf.connect([(insert_b0_into_dwi, datasink, [('out_dwi','dwi_b0_merge'), ('out_bvals','out_bvals'),('out_bvecs','out_bvecs')])])
    wf.connect([(mask_b0, datasink, [('mask_file','mask_b0')])])
    wf.connect([(b0_avg, datasink, [('out_file','b0_average')])])

    return wf



def hmc_pipeline(datasink_directory, name='motion_correct'):
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
    from clinica.pipeline.preprocessing.dwi_utils import merge_volumes_tdim
    from clinica.pipeline.preprocessing.dwi_utils import hmc_split
    from clinica.pipeline.preprocessing.dwi_corrections import dwi_flirt
    import nipype.interfaces.io as nio

    params = dict(dof=6, interp='spline', cost='normmi', cost_func='normmi', bins=50, save_log=True, padding_size=10,
                  schedule=get_flirt_schedule('hmc'),
                  searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4], fine_search=1, coarse_search=10 )

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file',
                        'in_bvec', 'in_bval', 'in_mask', 'ref_num']), name='inputnode')

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

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'hmc_correction/')

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file',
                         'out_bvec', 'out_xfms']), name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,        split,            [('in_file', 'in_file'),
                                              ('in_bval', 'in_bval'),
                                              ('ref_num', 'ref_num')]),
        (inputnode,        flirt,            [('in_mask', 'inputnode.ref_mask')]),
        (split,            flirt,            [('out_ref', 'inputnode.reference'),
                                              ('out_mov', 'inputnode.in_file'),
                                              ('out_bval', 'inputnode.in_bval')]),
        (flirt,            insmat,           [('outputnode.out_xfms', 'inlist')]),
        (split,            insmat,           [('volid', 'volid')]),
        (inputnode,        rot_bvec,         [('in_bvec', 'in_bvec')]),
        (insmat,           rot_bvec,         [('out', 'in_matrix')]),
        (rot_bvec,         outputnode,       [('out_file', 'out_bvec')]),
        (flirt,            merged_volumes,   [('outputnode.out_ref', 'in_file1'),
                                              ('outputnode.out_file', 'in_file2')]),
        (merged_volumes,   outputnode,       [('out_file', 'out_file')]),
        (insmat,           outputnode,       [('out', 'out_xfms')]),
        (merged_volumes,   datasink,         [('out_file', 'out_file')]),
        (insmat,           datasink,         [('out', 'out_xfms')]),
        (rot_bvec,         datasink,         [('out_file', 'out_bvec')])
    ])
    return wf




def ecc_pipeline(datasink_directory, name='eddy_correct'):
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
    >>> from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline
    >>> ecc = ecc_pipeline()
    >>> ecc.inputs.inputnode.in_file = 'diffusion.nii'
    >>> ecc.inputs.inputnode.in_bval = 'diffusion.bval'
    >>> ecc.inputs.inputnode.in_mask = 'mask.nii'
    >>> ecc.run() # doctest: +SKIP
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
    from clinica.pipeline.preprocessing.dwi_corrections import dwi_flirt
    from nipype.workflows.data import get_flirt_schedule
    from nipype.workflows.dmri.fsl.utils import extract_bval
    from nipype.workflows.dmri.fsl.utils import recompose_xfm
    from nipype.workflows.dmri.fsl.utils import recompose_dwi
    from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian
    import nipype.interfaces.io as nio

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

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_xfms']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'ecc_correction/')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  getb0,        [('in_file', 'in_file')]),
        (inputnode,  pick_dws,     [('in_file', 'in_dwi'),
                                    ('in_bval', 'in_bval')]),
        (inputnode,  merge,        [('in_file', 'in_dwi'),
                                    ('in_bval', 'in_bval')]),
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
        (merge,      outputnode,   [('out_file', 'out_file')]),
        (get_mat,    datasink,     [('out_files', 'out_xfms')]),
        (merge,      datasink,     [('out_file', 'out_file')])
    ])
    return wf



def sdc_fmb(datasink_directory, name='fmb_correction',
            fugue_params=dict(smooth3d=2.0),
            bmap_params=dict(delta_te=2.46e-3),
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
    import nipype.interfaces.io as nio

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file',
                        'in_mask', 'bmap_pha', 'bmap_mag']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file', 'out_vsm',
                         'out_warp']),
                         name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'sdc_correction/')

    getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')
    n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='n4_magnitude')
    bet = pe.Node(fsl.BET(frac=0.4, mask=True), name='bet_n4_magnitude')
    dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                     args='-kernel sphere 5 -dilM'), name='dilate_bet')
    pha2rads = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'],
                       function=siemens2rads), name='PreparePhase')
    prelude = pe.Node(fsl.PRELUDE(process3d=True), name='PhaseUnwrap')
    rad2rsec = pe.Node(niu.Function(input_names=['in_file', 'delta_te'],
                       output_names=['out_file'], function=rads2radsec), name='ToRadSec')
    rad2rsec.inputs.delta_te = bmap_params['delta_te']

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
    vsm.inputs.asym_se_time = bmap_params['delta_te']
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
        (inputnode,   pha2rads,     [('bmap_pha', 'in_file')]),
        (inputnode,   getb0,        [('in_file', 'in_file')]),
        (inputnode,   n4,           [('bmap_mag', 'input_image')]),
        (n4,          bet,          [('output_image', 'in_file')]),
        (bet,         dilate,       [('mask_file', 'in_file')]),
        (pha2rads,    prelude,      [('out_file', 'phase_file')]),
        (n4,          prelude,      [('output_image', 'magnitude_file')]),
        (dilate,      prelude,      [('out_file', 'mask_file')]),
        (prelude,     rad2rsec,     [('unwrapped_phase_file', 'in_file')]),
        (getb0,       flirt,        [('roi_file', 'reference')]),
        (inputnode,   flirt,        [('in_mask', 'ref_weight')]),
        (n4,          flirt,        [('output_image', 'in_file')]),
        (dilate,      flirt,        [('out_file', 'in_weight')]),
        (getb0,       applyxfm,     [('roi_file', 'reference')]),
        (rad2rsec,    applyxfm,     [('out_file', 'in_file')]),
        (flirt,       applyxfm,     [('out_matrix_file', 'in_matrix_file')]),
        (applyxfm,    pre_fugue,    [('out_file', 'fmap_in_file')]),
        (inputnode,   pre_fugue,    [('in_mask', 'mask_file')]),
        (pre_fugue,   demean,       [('fmap_out_file', 'in_file')]),
        (inputnode,   demean,       [('in_mask', 'in_mask')]),
        (demean,      cleanup,      [('out_file', 'inputnode.in_file')]),
        (inputnode,   cleanup,      [('in_mask', 'inputnode.in_mask')]),
        (cleanup,     addvol,       [('outputnode.out_file', 'in_file')]),
        (inputnode,   vsm,          [('in_mask', 'mask_file')]),
        (addvol,      vsm,          [('out_file', 'fmap_in_file')]),
        (inputnode,   split,        [('in_file', 'in_file')]),
        (split,       unwarp,       [('out_files', 'in_file')]),
        (vsm,         unwarp,       [('shift_out_file', 'shift_in_file')]),
        (unwarp,      thres,        [('unwarped_file', 'in_file')]),
        (thres,       merge,        [('out_file', 'in_files')]),
        (merge,       vsm2dfm,      [('merged_file', 'inputnode.in_ref')]),
        (vsm,         vsm2dfm,      [('shift_out_file', 'inputnode.in_vsm')]),
        (merge,       outputnode,   [('merged_file', 'out_file')]),
        (vsm,         outputnode,   [('shift_out_file', 'out_vsm')]),
        (vsm2dfm,     outputnode,   [('outputnode.out_warp', 'out_warp')]),
        (merge,       datasink,     [('merged_file', 'out_file')]),
        (vsm,         datasink,     [('shift_out_file', 'out_vsm')]),
        (vsm2dfm,     datasink,     [('outputnode.out_warp', 'out_warp')])
    ])
    return wf



def sdc_syb_pipeline(datasink_directory, name='sdc_syb_correct'):

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
    outputnode.B0_2_T1_rigid_body_matrix - B0 to T1 image FLIRT rigid body fsl coregistration matrix
    outputnode.T1_coregistered_2_B0 - T1 image rigid body coregistered to the B0 image
    outputnode.B0_2_T1_affine_matrix - B0 to T1 image ANTs affine itk coregistration matrix
    outputnode.B0_2_T1_SyN_defomation_field - B0 to T1 image ANTs SyN itk warp
    outputnode.out_warp - Out warp allowing DWI to T1 registration and susceptibilty induced artifacts correction

    Example
    -------
    >>> epi = epi_pipeline()
    >>> epi.inputs.inputnode.DWI = 'DWI.nii'
    >>> epi.inputs.inputnode.T1 = 'T1.nii'
    >>> epi.run() # doctest: +SKIP
    """

    import nipype.interfaces.io as nio
    import nipype.interfaces.ants as ants
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl as fsl

    def antsRegistrationSyNQuick(fixe_image, moving_image):

        import subprocess
        import os.path as op

        image_warped = op.abspath('SyN_QuickWarped.nii.gz')
        affine_matrix = op.abspath('SyN_Quick0GenericAffine.mat')
        warp = op.abspath('SyN_Quick1Warp.nii.gz')
        inverse_warped = op.abspath('SyN_QuickInverseWarped.nii.gz')
        inverse_warp = op.abspath('SyN_Quick1InverseWarp.nii.gz')

        cmd = 'antsRegistrationSyNQuick.sh -t br -d 3 -f ' + fixe_image + ' -m ' + moving_image + ' -o SyN_Quick'
        subprocess.call([cmd], shell=True)

        return image_warped, affine_matrix, warp, inverse_warped, inverse_warp

    def antscombintransform(in_file, transforms_list, reference):

        import os
        import os.path as op

        out_warp = op.abspath('out_warp.nii.gz')

        transforms = ""
        for trans in transforms_list:
            transforms += " " + trans
        cmd = 'antsApplyTransforms -o [out_warp.nii.gz,1] -i ' + in_file + ' -r ' + reference + ' -t' + transforms
        os.system(cmd)

        return out_warp

    inputnode = pe.Node(niu.IdentityInterface(fields=['T1', 'DWI']), name='inputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_2_T1 = pe.Node(interface=fsl.FLIRT(dof=6), name = 'flirt_B0_2_T1')
    flirt_b0_2_T1.inputs.interp = "spline"
    flirt_b0_2_T1.inputs.cost = 'normmi'
    flirt_b0_2_T1.inputs.cost_func = 'normmi'

    invertxfm = pe.Node(interface=fsl.ConvertXFM(), name='invert_xfm')
    invertxfm.inputs.invert_xfm = True

    apply_xfm = pe.Node(interface=fsl.ApplyXfm(), name='apply_xfm')
    apply_xfm.inputs.apply_xfm = True
    apply_xfm.inputs.interp = "spline"
    apply_xfm.inputs.cost = 'normmi'
    apply_xfm.inputs.cost_func = 'normmi'

    antsRegistrationSyNQuick = pe.Node(interface=niu.Function(input_names=['fixe_image', 'moving_image'], output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
                                                              function=antsRegistrationSyNQuick), name='antsRegistrationSyNQuick')

    merge_transform = pe.Node(niu.Merge(2), name='MergeTransforms')

    combin_warp = pe.Node(interface=niu.Function(input_names=['in_file', 'transforms_list', 'reference'], output_names=['out_warp'],
                                                    function=antscombintransform), name='combin_warp')

    coeffs = pe.Node(fsl.WarpUtils(out_format='spline'), name='CoeffComp')

    fsl_transf = pe.Node(fsl.WarpUtils(out_format='field'), name='fsl_transf')

    apply_warp = pe.MapNode(interface=fsl.ApplyWarp(), iterfield=['in_file'],name='apply_warp')
    apply_warp.inputs.interp = 'spline'

    thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                       name='RemoveNegative')

    merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')

    outputnode = pe.Node(niu.IdentityInterface(fields=['B0_2_T1_rigid_body_matrix', 'T1_2_B0_rigid_body_matrix',
                                                       'T1_coregistered_2_B0', 'B0_2_T1_SyN_defomation_field',
                                                       'B0_2_T1_affine_matrix', 'out_dwi', 'out_warp']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'epi_correction/')

    wf = pe.Workflow(name='sdc_syb_pipeline')

    wf.connect([(inputnode, split,[('DWI','in_file')])])
    wf.connect([(split, pick_ref, [('out_files','inlist')])])
    wf.connect([(pick_ref, flirt_b0_2_T1, [('out','in_file')])])
    wf.connect([(inputnode, flirt_b0_2_T1, [('T1','reference')])])
    wf.connect([(flirt_b0_2_T1, invertxfm, [('out_matrix_file','in_file')])])
    wf.connect([(invertxfm, apply_xfm, [('out_file','in_matrix_file')])])
    wf.connect([(inputnode, apply_xfm, [('T1','in_file')])])
    wf.connect([(pick_ref, apply_xfm, [('out','reference')])])
    wf.connect([(apply_xfm, antsRegistrationSyNQuick, [('out_file','fixe_image')])])
    wf.connect([(pick_ref, antsRegistrationSyNQuick,[('out','moving_image')])])
    wf.connect([(antsRegistrationSyNQuick, merge_transform, [('affine_matrix','in2'), ('warp','in1')])])
    wf.connect([(pick_ref, combin_warp, [('out','in_file')])])
    wf.connect([(merge_transform, combin_warp, [('out','transforms_list')])])
    wf.connect([(apply_xfm, combin_warp, [('out_file','reference')])])
    wf.connect([(apply_xfm, coeffs, [('out_file', 'reference')])])
    wf.connect([(combin_warp, coeffs, [('out_warp', 'in_file')])])
    wf.connect([(coeffs, fsl_transf, [('out_file', 'in_file')])])
    wf.connect([(apply_xfm, fsl_transf, [('out_file', 'reference')])])
    wf.connect([(fsl_transf, apply_warp, [('out_file','field_file')])])
    wf.connect([(split, apply_warp, [('out_files','in_file')])])
    wf.connect([(apply_xfm, apply_warp, [('out_file','ref_file')])])
    wf.connect([(apply_warp, thres, [('out_file','in_file')])])
    wf.connect([(thres, merge, [('out_file','in_files')])])
    wf.connect([(merge, outputnode, [('merged_file','out_dwi')])])
    wf.connect([(flirt_b0_2_T1, outputnode, [('out_matrix_file','B0_2_T1_rigid_body_matrix')])])
    wf.connect([(invertxfm, outputnode, [('out_file', 'T1_2_B0_rigid_body_matrix')])])
    wf.connect([(apply_xfm, outputnode, [('out_file','T1_coregistered_2_B0')])])
    wf.connect([(antsRegistrationSyNQuick, outputnode, [('warp','B0_2_T1_SyN_defomation_field'),
                                                        ('affine_matrix','B0_2_T1_affine_matrix')])])
    wf.connect([(fsl_transf, outputnode, [('out_file','out_warp')])])
    wf.connect([(merge, datasink, [('merged_file','out_dwi')])])
    wf.connect([(flirt_b0_2_T1, datasink, [('out_matrix_file','B0_2_T1_rigid_body_matrix')])])
    wf.connect([(invertxfm, datasink, [('out_file', 'T1_2_B0_rigid_body_matrix')])])
    wf.connect([(apply_xfm, datasink, [('out_file','T1_coregistered_2_B0')])])
    wf.connect([(antsRegistrationSyNQuick, datasink, [('warp','B0_2_T1_SyN_defomation_field'),
                                                      ('affine_matrix','B0_2_T1_affine_matrix'),
                                                      ('image_warped','b0_warped_image')])])
    wf.connect([(fsl_transf, datasink, [('out_file','out_warp')])])

    return wf



def apply_all_corrections_syb(datasink_directory, name='UnwarpArtifacts'):
    """
    Combines two lists of linear transforms with the deformation field
    map obtained epi_correction by Ants.
    Additionally, computes the corresponding bspline coefficients and
    the map of determinants of the jacobian.
    """

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_sdc_syb', 'in_hmc','in_ecc', 'in_dwi', 'T1']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_warp', 'out_coeff', 'out_jacobian']),
        name='outputnode')

    split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')

    pick_ref = pe.Node(niu.Select(), name='Pick_b0')
    pick_ref.inputs.index = [0]

    flirt_b0_2_T1 = pe.Node(interface=fsl.FLIRT(dof=6), name = 'flirt_B0_2_T1')
    flirt_b0_2_T1.inputs.interp = "spline"
    flirt_b0_2_T1.inputs.cost = 'normmi'
    flirt_b0_2_T1.inputs.cost_func = 'normmi'

    invertxfm = pe.Node(interface=fsl.ConvertXFM(), name='invert_xfm')
    invertxfm.inputs.invert_xfm = True

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

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'apply_all_correction/')

    wf = pe.Workflow(name=name)
    wf.connect([(inputnode, concat_hmc_ecc, [('in_ecc','in_file2')]),
                (inputnode, concat_hmc_ecc, [('in_hmc','in_file')]),
                (concat_hmc_ecc, warps, [('out_file','premat')]),
                (inputnode, warps, [('in_sdc_syb','warp1')]),
                (inputnode, split, [('in_dwi', 'in_file')]),
                (split, pick_ref, [('out_files','inlist')]),
                (pick_ref, flirt_b0_2_T1, [('out','in_file')]),
                (inputnode, flirt_b0_2_T1, [('T1','reference')]),
                (flirt_b0_2_T1, invertxfm, [('out_matrix_file','in_file')]),
                (invertxfm, apply_xfm, [('out_file','in_matrix_file')]),
                (inputnode, apply_xfm, [('T1','in_file')]),
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

    wf.connect([(warps, datasink, [('out_file', 'out_warp')]),
                (coeffs, datasink, [('out_file', 'out_coeff')]),
                (jacobian, datasink, [('out_jacobian', 'out_jacobian')]),
                (merge, datasink, [('merged_file', 'out_file')])])

    return wf


# 
# def dwi_flirt(name='DWICoregistration', excl_nodiff=False, flirt_param={}):
#     """
#     Generates a workflow for linear registration of dwi volumes using flirt.
#
#     Inputnode
#     ---------
#     reference : FILE
#       Mandatory input. Reference data set.
#     in_file : FILE
#       Mandatory input. Moving data set.
#     ref_mask : FILE
#       Mandatory input. Binary mask of the reference volume.
#     in_xfms : FILE
#       Mandatory input. Intialisation matrices for flirt.
#     in_bval : FILE
#       Mandatory input. B values file.
#
#     """
#     from nipype.workflows.dmri.fsl.utils import _checkinitxfm
#     from nipype.workflows.dmri.fsl.utils import enhance
#
#     inputnode = pe.Node(niu.IdentityInterface(fields=['reference',
#                         'in_file', 'ref_mask', 'in_xfms', 'in_bval']),
#                         name='inputnode')
#
#     initmat = pe.Node(niu.Function(input_names=['in_bval', 'in_xfms',
#                       'excl_nodiff'], output_names=['init_xfms'],
#                                    function=_checkinitxfm), name='InitXforms')
#     initmat.inputs.excl_nodiff = excl_nodiff
#     dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
#                      args='-kernel sphere 5 -dilM'), name='MskDilate')
#     split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
#     pick_ref = pe.Node(niu.Select(), name='Pick_b0')
#     n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3), name='Bias')
#     enhb0 = pe.Node(niu.Function(input_names=['in_file', 'in_mask',
#                     'clip_limit'], output_names=['out_file'],
#                                  function=enhance), name='B0Equalize')
#     enhb0.inputs.clip_limit = 0.015
#     enhdw = pe.MapNode(niu.Function(input_names=['in_file', 'in_mask'],
#                        output_names=['out_file'], function=enhance),
#                        name='DWEqualize', iterfield=['in_file'])
#     flirt = pe.MapNode(fsl.FLIRT(**flirt_param), name='CoRegistration',
#                        iterfield=['in_file', 'in_matrix_file'])
#     thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
#                        name='RemoveNegative')
#     merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
#     outputnode = pe.Node(niu.IdentityInterface(fields=['out_file',
#                          'out_xfms', 'out_ref']), name='outputnode')
#
#     wf = pe.Workflow(name=name)
#     wf.connect([
#         (inputnode,   split,        [('in_file', 'in_file')]),
#         (inputnode,   dilate,       [('ref_mask', 'in_file')]),
#         (inputnode,   enhb0,        [('ref_mask', 'in_mask')]),
#         (inputnode,   initmat,      [('in_xfms', 'in_xfms'),
#                                      ('in_bval', 'in_bval')]),
#         (inputnode,   n4,           [('reference', 'input_image'),
#                                      ('ref_mask', 'mask_image')]),
#         (dilate,      flirt,        [('out_file', 'ref_weight'),
#                                      ('out_file', 'in_weight')]),
#         (n4,          enhb0,        [('output_image', 'in_file')]),
#         (split,       enhdw,        [('out_files', 'in_file')]),
#         (dilate,      enhdw,        [('out_file', 'in_mask')]),
#         (enhb0,       flirt,        [('out_file', 'reference')]),
#         (enhdw,       flirt,        [('out_file', 'in_file')]),
#         (initmat,     flirt,        [('init_xfms', 'in_matrix_file')]),
#         (flirt,       thres,        [('out_file', 'in_file')]),
#         (thres,       merge,        [('out_file', 'in_files')]),
#         (merge,       outputnode,   [('merged_file', 'out_file')]),
#         (enhb0,       outputnode,   [('out_file', 'out_ref')]),
#         (flirt,       outputnode,   [('out_matrix_file', 'out_xfms')])
#     ])
#     return wf


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
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,  split,      [('in_file', 'in_file')]),
        (inputnode,  dilate,     [('ref_mask', 'in_file')]),
        (inputnode,   n4,        [('reference', 'input_image'),
                                  ('ref_mask', 'mask_image')]),
#        (inputnode,  flirt,      [('ref_mask', 'reference')]),
        (n4,         flirt,      [('output_image', 'reference')]),
        (inputnode,  initmat,    [('in_xfms', 'in_xfms'),
                                  ('in_bval', 'in_bval')]),
        (dilate,     flirt,      [('out_file', 'ref_weight'),
                                  ('out_file', 'in_weight')]),
        (split,      flirt,      [('out_files', 'in_file')]),
        (initmat,    flirt,      [('init_xfms', 'in_matrix_file')]),
        (flirt,      thres,      [('out_file', 'in_file')]),
        (thres,      merge,      [('out_file', 'in_files')]),
        (merge,      outputnode, [('merged_file', 'out_file')]),
        (inputnode,      outputnode, [('reference', 'out_ref')]),
        (flirt,      outputnode, [('out_matrix_file', 'out_xfms')])
    ])
    return wf


def remove_bias(datasink_directory, name='bias_correct'):
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

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'bias_correction/')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode,    getb0,          [('in_file', 'in_file')]),
        (getb0,       n4,              [('roi_file', 'input_image')]),
        (getb0,       mask_b0,         [('roi_file', 'in_file')]),
        (mask_b0,    n4,               [('mask_file', 'mask_image')]),
        (inputnode,    split,          [('in_file', 'in_file')]),
        (n4,           mult,           [('bias_image', 'operand_files')]),
        (split,        mult,           [('out_files', 'in_file')]),
        (mult,         thres,          [('out_file', 'in_file')]),
        (thres,        merge,          [('out_file', 'in_files')]),
        (merge,        outputnode,     [('merged_file', 'out_file')]),
        (mask_b0,    outputnode,       [('mask_file', 'b0_mask')]),
        (merge,        datasink,       [('merged_file', 'out_file')]),
        (mask_b0,    datasink,         [('mask_file', 'b0_mask')]),
    ])
    return wf


# def remove_bias(datasink_directory, name='bias_correct'):
#     """
#     This workflow estimates a single multiplicative bias field from the
#     averaged *b0* image, as suggested in [Jeurissen2014]_.

#     .. warning :: Suppose b0 is the first volume.

#     .. admonition:: References

#       .. [Jeurissen2014] Jeurissen B. et al., `Multi-tissue constrained
#         spherical deconvolution for improved analysis of multi-shell diffusion
#         MRI data <http://dx.doi.org/10.1016/j.neuroimage.2014.07.061>`_.
#         NeuroImage (2014). doi: 10.1016/j.neuroimage.2014.07.061


#     Example
#     -------

#     >>> from nipype.workflows.dmri.fsl.artifacts import remove_bias
#     >>> bias = remove_bias()
#     >>> bias.inputs.inputnode.in_file = 'epi.nii'
#     >>> bias.inputs.inputnode.in_mask = 'mask.nii'
#     >>> bias.run() # doctest: +SKIP

#     """
#     import nipype.interfaces.io as nio

#     inputnode = pe.Node(niu.IdentityInterface(fields=['in_file',
#                         'in_mask']), name='inputnode')

#     outputnode = pe.Node(niu.IdentityInterface(fields=['out_file']),
#                          name='outputnode')

#    datasink = pe.Node(nio.DataSink(), name='datasink')
#    datasink.inputs.base_directory = op.join(datasink_directory, 'bias_correction/')

#     get_b0 = pe.Node(fsl.ExtractROI(args='0 1'), name='get_b0')

#     n4 = pe.Node(ants.N4BiasFieldCorrection(dimension=3,
#                  save_bias=True, bspline_fitting_distance=600), name='Bias_b0')
#     split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
#     mult = pe.MapNode(fsl.MultiImageMaths(op_string='-div %s'),
#                       iterfield=['in_file'], name='RemoveBiasOfDWIs')
#     thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
#                        name='RemoveNegative')
#     merge = pe.Node(fsl.utils.Merge(dimension='t'), name='MergeDWIs')

#     wf = pe.Workflow(name=name)
#     wf.connect([
#         (inputnode,   get_b0,       [('in_file', 'in_file')]),
#         (get_b0,      n4,           [('roi_file', 'input_image')]),
#         (inputnode,   n4,           [('in_mask', 'mask_image')]),
#         (inputnode,   split,        [('in_file', 'in_file')]),
#         (n4,          mult,         [('bias_image', 'operand_files')]),
#         (split,       mult,         [('out_files', 'in_file')]),
#         (mult,        thres,        [('out_file', 'in_file')]),
#         (thres,       merge,        [('out_file', 'in_files')]),
#         (merge,       outputnode,   [('merged_file', 'out_file')]),
#         (merge,       datasink,     [('merged_file', 'out_file')])
#     ])
#     return wf
