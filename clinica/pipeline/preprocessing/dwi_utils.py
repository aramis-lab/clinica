#!/usr/bin/python

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl

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

    assert(op.isfile(in_file))
    assert(op.isfile(in_bvals))
    assert(op.isfile(in_bvecs))
    assert(lowbval >= 0)
    
    im = nib.load(in_file)
    data = im.get_data()
    hdr = im.get_header().copy()
    bvals = np.loadtxt(in_bvals)
    bvecs = np.loadtxt(in_bvecs)

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



def b0_flirt_pipeline(name='b0coregistration', nb_b0= 5, excl_nodiff=False):
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
    from dwi_utils import merge_volumes_tdim
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
        (inputnode,  flirt,      [('ref_mask', 'reference')]),
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
