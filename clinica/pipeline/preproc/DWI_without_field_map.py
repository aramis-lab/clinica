# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:17:41 2016

@author: jacquemont
"""

def preprocess_pipeline(subject, files_directory, working_directory, datasink_directory):
    
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants
    from nipype.workflows.dmri.fsl.utils import insert_mat
    from nipype.workflows.dmri.fsl.utils import recompose_xfm
    from nipype.workflows.dmri.fsl.utils import recompose_dwi
    from nipype.workflows.dmri.fsl.artifacts import _xfm_jacobian
    from nipype.workflows.dmri.fsl.utils import extract_bval
    import os
    import os.path as op

# Inputs existence checking

    inputs=[files_directory, working_directory, datasink_directory]     
        
    for input_file in inputs:
        if not op.exists(input_file):
            raise IOError('file {} does not exist'.format(input_file))

# Utilities
    
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
    
    
    def b0_dwi_split(in_file, in_bvals, in_bvecs, lowbval=0.0):
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
        bvecs = np.loadtxt(in_bvecs).T
    
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
    
        import os.path as op
        import numpy as np 
    
        assert(op.isfile(in_b0))
        assert(op.isfile(in_dwi))
        assert(op.isfile(in_bvals))
        assert(op.isfile(in_bvecs))
        
        def merge_volumes_tdim(in_file1, in_file2):
       
            import os.path as op
            import os
    
            out_file = op.abspath('merged_files.nii.gz')
            cmd = 'fslmerge -t '+ out_file + ' ' + in_file1 + ' ' + in_file2
            os.system(cmd)
            return out_file
    
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
    
        See
        ---
        :func:`epynet.preprocessing.pipeline.hmc_pipeline`
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

# Correction pipelines definition

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
            (inputnode,     split,   [('in_file', 'in_file'),
                                      ('in_bval', 'in_bval'),
                                      ('ref_num', 'ref_num')]),
            (inputnode,  flirt,      [('in_mask', 'inputnode.ref_mask')]),
            (split,      flirt,      [('out_ref', 'inputnode.reference'),
                                      ('out_mov', 'inputnode.in_file'),
                                      ('out_bval', 'inputnode.in_bval')]),
            (flirt,      insmat,     [('outputnode.out_xfms', 'inlist')]),
            (split,      insmat,     [('volid', 'volid')]),
            (inputnode,  rot_bvec,   [('in_bvec', 'in_bvec')]),
            (insmat,     rot_bvec,   [('out', 'in_matrix')]),
            (rot_bvec,   outputnode, [('out_file', 'out_bvec')]),
            (flirt,      merged_volumes, [('outputnode.out_ref', 'in_file1')]),
            (flirt,      merged_volumes, [('outputnode.out_file', 'in_file2')]),
            (merged_volumes,      outputnode, [('out_file', 'out_file')]),
            (insmat,     outputnode, [('out', 'out_xfms')]),
            (merged_volumes,      datasink, [('out_file', 'out_file')]),
            (insmat,     datasink, [('out', 'out_xfms')]),
            (rot_bvec,   datasink, [('out_file', 'out_bvec')])
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
    
        from nipype.workflows.data import get_flirt_schedule
        params = dict(dof=12, no_search=True, interp='spline', bgvalue=0,
                      schedule=get_flirt_schedule('ecc'))
    
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['in_file', 'in_bval', 'in_mask', 'in_xfms']), name='inputnode')
            
        getb0 = pe.Node(fsl.ExtractROI(t_min=0, t_size=1), name='get_b0')
        
        pick_dws = pe.Node(niu.Function(
            input_names=['in_dwi', 'in_bval', 'b'], output_names=['out_file'],
            function=extract_bval), name='ExtractDWI')
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
            
        datasink = pe.Node(nio.DataSink(), name='datasink')
        datasink.inputs.base_directory = op.join(datasink_directory, 'ecc_correction/')
    
        outputnode = pe.Node(niu.IdentityInterface(
            fields=['out_file', 'out_xfms']), name='outputnode')
    
        wf = pe.Workflow(name=name)
        wf.connect([
            (inputnode,  getb0,     [('in_file', 'in_file')]),
            (inputnode,  pick_dws,   [('in_file', 'in_dwi'),
                                      ('in_bval', 'in_bval')]),
            (inputnode,  merge,      [('in_file', 'in_dwi'),
                                      ('in_bval', 'in_bval')]),
            (inputnode,  flirt,      [('in_mask', 'inputnode.ref_mask'),
                                      ('in_xfms', 'inputnode.in_xfms'),
                                      ('in_bval', 'inputnode.in_bval')]),
            (inputnode,  get_mat,    [('in_bval', 'in_bval')]),
            (getb0,     flirt,      [('roi_file', 'inputnode.reference')]),
            (pick_dws,   flirt,      [('out_file', 'inputnode.in_file')]),
            (flirt,      get_mat,    [('outputnode.out_xfms', 'in_xfms')]),
            (flirt,      mult,       [(('outputnode.out_xfms', _xfm_jacobian),
                                       'operand_value')]),
            (flirt,      split,      [('outputnode.out_file', 'in_file')]),
            (split,      mult,       [('out_files', 'in_file')]),
            (mult,       thres,      [('out_file', 'in_file')]),
            (thres,      merge,      [('out_file', 'in_corrected')]),
            (get_mat,    outputnode, [('out_files', 'out_xfms')]),
            (merge,      outputnode, [('out_file', 'out_file')]),
            (get_mat,    datasink, [('out_files', 'out_xfms')]),
            (merge,      datasink, [('out_file', 'out_file')])
        ])
        return wf
    
    def epi_pipeline():
        
        """
        This workflow allows to correct for echo-planare induced susceptibility artifacts without fieldmap 
        (e.g. ADNI Database) by elastically register DWIs to their respective baseline T1-weighted 
        structural scans using an inverse consistent registration algorithm with a mutual information cost 
        function (SyN algorithm). This workflow allows also a coregistration of DWIs with their respective 
        baseline T1-weighted structural scans in order to latter combine tracks and cortex parcelation.
        ..  warning:: This workflow rotates the `b`-vectors' 
        .. References
          .. Nir et al. (Neurobiology of Aging 2015)- Connectivity network measures predict volumetric atrophy in mild cognitive impairment
            
            Leow et al. (IEEE Trans Med Imaging 2007)- Statistical Properties of Jacobian Maps and the Realization of Unbiased Large Deformation Nonlinear Image Registration

        Inputnode
        ---------
        DWI : FILE
          Mandatory input. Input dwi file.
        bvec : FILE
          Mandatory input. Vector file of the diffusion directions of the dwi dataset.
        T1 : FILE
          Mandatory input. Input T1 file.

        Outputnode
        ----------
    
        outputnode.out_dwi - corrected dwi file
        outputnode.out_bvec - rotated gradient vectors table
        outputnode.B0_2_T1_rigid_body_matrix - B0 to T1 image FLIRT rigid body fsl coregistration matrix
        outputnode.B0_2_T1_affine_matrix - B0 to T1 image ANTs affine itk coregistration matrix
        outputnode.B0_2_T1_SyN_defomation_field - B0 to T1 image ANTs SyN itk warp
        outputnode.out_warp - Out warp allowing DWI to T1 registration and susceptibilty induced artifacts correction

        Example
        -------
        >>> epi = epi_pipeline()
        >>> epi.inputs.inputnode.DWI = 'DWI.nii'
        >>> epi.inputs.inputnode.bvec = 'bvec.txt'
        >>> epi.inputs.inputnode.T1 = 'T1.nii'
        >>> epi.run() # doctest: +SKIP
        """
        
        import nipype.interfaces.io as nio
        import nipype.interfaces.ants as ants
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as niu
        import nipype.interfaces.fsl as fsl
        
        def expend_matrix_list(in_matrix, in_bvec):
            
            import numpy as np
            
            bvecs = np.loadtxt(in_bvec).T
            out_matrix_list = [in_matrix]
            
            out_matrix_list = out_matrix_list * len(bvecs)
        
            return out_matrix_list
        
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
    
        inputnode = pe.Node(niu.IdentityInterface(fields=['T1', 'DWI', 'bvec']), name='inputnode')
        
        split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
        pick_ref = pe.Node(niu.Select(), name='Pick_b0') 
        pick_ref.inputs.index = [0]
        
        flirt_b0_2_T1 = pe.Node(interface=fsl.FLIRT(dof=6), name = 'flirt_B0_2_T1')
        flirt_b0_2_T1.inputs.interp = "spline"
        flirt_b0_2_T1.inputs.cost = 'normmi'
        flirt_b0_2_T1.inputs.cost_func = 'normmi'
    
        apply_xfm = pe.Node(interface=fsl.ApplyXfm(), name='apply_xfm')
        apply_xfm.inputs.apply_xfm = True
        
        expend_matrix = pe.Node(interface=niu.Function(input_names=['in_matrix', 'in_bvec'], output_names=['out_matrix_list'], function=expend_matrix_list), name='expend_matrix')    
        
        rot_bvec = pe.Node(niu.Function(input_names=['in_matrix','in_bvec'],
                           output_names=['out_file'], function=rotate_bvecs),
                           name='Rotate_Bvec')
        
        antsRegistrationSyNQuick = pe.Node(interface=niu.Function(input_names=['fixe_image', 'moving_image'], output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'], 
                                                                  function=antsRegistrationSyNQuick), name='antsRegistrationSyNQuick')
        
        merge_transform = pe.Node(niu.Merge(2), name='MergeTransforms')
        
        combin_warp = pe.Node(interface=niu.Function(input_names=['in_file', 'transforms_list', 'reference'], output_names=['out_warp'], 
                                                        function=antscombintransform), name='combin_warp')
        
        coeffs = pe.Node(fsl.WarpUtils(out_format='spline'), name='CoeffComp')
        
        fsl_transf = pe.Node(fsl.WarpUtils(out_format='field'), name='fsl_transf')
        
        warp_epi = pe.Node(fsl.ConvertWarp(), name='warp_epi')
        
        apply_warp = pe.MapNode(interface=fsl.ApplyWarp(), iterfield=['in_file'],name='apply_warp')
        apply_warp.inputs.interp = 'spline'
        
        thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                           name='RemoveNegative')
        
        merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
        
        outputnode = pe.Node(niu.IdentityInterface(fields=['B0_2_T1_rigid_body_matrix', 'out_bvec',
                                                           'B0_2_T1_SyN_defomation_field', 'B0_2_T1_affine_matrix',
                                                           'out_dwi', 'out_warp']), name='outputnode')
        
        datasink = pe.Node(nio.DataSink(), name='datasink')
        datasink.inputs.base_directory = op.join(datasink_directory,'epi_correction/')
        
        wf = pe.Workflow(name='epi_pipeline')
        
        wf.connect([(inputnode, split,[('DWI','in_file')])])
        wf.connect([(split, pick_ref, [('out_files','inlist')])])
        wf.connect([(pick_ref, flirt_b0_2_T1, [('out','in_file')])])
        wf.connect([(inputnode, flirt_b0_2_T1, [('T1','reference')])])
        wf.connect([(inputnode, rot_bvec, [('bvec', 'in_bvec')])])
        wf.connect([(flirt_b0_2_T1, expend_matrix, [('out_matrix_file','in_matrix')])])
        wf.connect([(inputnode, expend_matrix, [('bvec','in_bvec')])])
        wf.connect([(expend_matrix, rot_bvec, [('out_matrix_list','in_matrix')])])
        wf.connect([(inputnode, antsRegistrationSyNQuick, [('T1','fixe_image')])])
        wf.connect([(flirt_b0_2_T1, antsRegistrationSyNQuick,[('out_file','moving_image')])])
        wf.connect([(antsRegistrationSyNQuick, merge_transform, [('affine_matrix','in2'), ('warp','in1')])])
        wf.connect([(flirt_b0_2_T1, combin_warp, [('out_file','in_file')])])
        wf.connect([(merge_transform, combin_warp, [('out','transforms_list')])])
        wf.connect([(inputnode, combin_warp, [('T1','reference')])])   
        wf.connect([(inputnode, coeffs, [('T1', 'reference')])])
        wf.connect([(combin_warp, coeffs, [('out_warp', 'in_file')])])
        wf.connect([(coeffs, fsl_transf, [('out_file', 'in_file')])])
        wf.connect([(inputnode, fsl_transf, [('T1', 'reference')])])
        wf.connect([(inputnode, warp_epi, [('T1','reference')])])
        wf.connect([(flirt_b0_2_T1, warp_epi, [('out_matrix_file','premat')])])
        wf.connect([(fsl_transf, warp_epi, [('out_file','warp1')])])
        wf.connect([(warp_epi, apply_warp, [('out_file','field_file')])])
        wf.connect([(split, apply_warp, [('out_files','in_file')])])
        wf.connect([(inputnode, apply_warp, [('T1','ref_file')])])
        wf.connect([(apply_warp, thres, [('out_file','in_file')])])
        wf.connect([(thres, merge, [('out_file','in_files')])])
        wf.connect([(merge, outputnode, [('merged_file','out_dwi')])])
        wf.connect([(flirt_b0_2_T1, outputnode, [('out_matrix_file','B0_2_T1_rigid_body_matrix')])])
        wf.connect([(antsRegistrationSyNQuick, outputnode, [('warp','B0_2_T1_SyN_defomation_field'),
                                                            ('affine_matrix','B0_2_T1_affine_matrix')])])
        wf.connect([(warp_epi, outputnode, [('out_file','out_warp')])])
        wf.connect([(rot_bvec, outputnode, [('out_file','out_bvec')])])
        wf.connect([(merge, datasink, [('merged_file','out_dwi')])])
        wf.connect([(flirt_b0_2_T1, datasink, [('out_matrix_file','B0_2_T1_rigid_body_matrix'), ('out_file','flirt_out_file')])])
        wf.connect([(antsRegistrationSyNQuick, datasink, [('warp','B0_2_T1_SyN_defomation_field'),
                                                          ('affine_matrix','B0_2_T1_affine_matrix'),
                                                          ('image_warped','b0_warped_image')])])
        wf.connect([(warp_epi, datasink, [('out_file','out_warp')])])
        wf.connect([(rot_bvec, datasink, [('out_file','out_bvec')])])
    
        return wf
    
    def apply_all_corrections(name='UnwarpArtifacts'):
        """
        Combines two lists of linear transforms with the deformation field
        map obtained epi_correction by Ants.
        Additionally, computes the corresponding bspline coefficients and
        the map of determinants of the jacobian.
        """
            
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['in_epi', 'in_hmc','in_ecc', 'in_dwi', 'T1']), name='inputnode')
        outputnode = pe.Node(niu.IdentityInterface(
            fields=['out_file', 'out_warp', 'out_coeff', 'out_jacobian']),
            name='outputnode')
            
        split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
                    
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
                    (inputnode, warps, [('in_epi','warp1')]),
                    (inputnode, warps, [('T1','reference')]),
                    (inputnode, split, [('in_dwi', 'in_file')]),
                    (warps, unwarp, [('out_file', 'field_file')]),
                    (split, unwarp, [('out_files', 'in_file')]),
                    (inputnode, unwarp, [('T1', 'ref_file')]),
                    (inputnode, coeffs, [('T1', 'reference')]),
                    (warps, coeffs, [('out_file', 'in_file')]),
                    (inputnode, jacobian, [('T1', 'reference')]),
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
    
    def orig_flirt(name='DWICoregistration', excl_nodiff=False,
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
        from nipype.workflows.data import get_flirt_schedule
    
        inputnode = pe.Node(niu.IdentityInterface(fields=['b0_avg_reference',
                            'in_file', 'ref_mask']),
                            name='inputnode')
        
        flirt_param = dict(dof=6, interp='spline', cost='normmi', cost_func='normmi', bins=50, save_log=True, padding_size=10,
                      schedule=get_flirt_schedule('hmc'),
                      searchr_x=[-4, 4], searchr_y=[-4, 4], searchr_z=[-4, 4], fine_search=1, coarse_search=10 )
        
        dilate = pe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                         args='-kernel sphere 5 -dilM'), name='MskDilate')
        split = pe.Node(fsl.Split(dimension='t'), name='SplitDWIs')
    
        flirt = pe.MapNode(fsl.FLIRT(**flirt_param), name='CoRegistration',
                           iterfield=['in_file'])
        thres = pe.MapNode(fsl.Threshold(thresh=0.0), iterfield=['in_file'],
                           name='RemoveNegative')
        merge = pe.Node(fsl.Merge(dimension='t'), name='MergeDWIs')
        outputnode = pe.Node(niu.IdentityInterface(fields=['out_file',
                             'out_xfms', 'out_ref']), name='outputnode')
        wf = pe.Workflow(name=name)
        wf.connect([
            (inputnode,  split,      [('in_file', 'in_file')]),
            (inputnode,  dilate,     [('ref_mask', 'in_file')]),
            (inputnode,  flirt,      [('b0_avg_reference', 'reference')]),
            (dilate,     flirt,      [('out_file', 'ref_weight'),
                                      ('out_file', 'in_weight')]),
            (split,      flirt,      [('out_files', 'in_file')]),                         
            (flirt,      thres,      [('out_file', 'in_file')]),
            (thres,      merge,      [('out_file', 'in_files')]),
            (merge,      outputnode, [('merged_file', 'out_file')]),
            (inputnode,      outputnode, [('b0_avg_reference', 'out_ref')]),
            (flirt,      outputnode, [('out_matrix_file', 'out_xfms')])
        ])
        return wf

# Higher level pipeline definition

    # Datagrabbing and input node    

    subject_list = [subject]
    
    infosource = pe.Node(interface=niu.IdentityInterface(fields=['subject_id']), name="infosource")
    infosource.iterables = ('subject_id', subject_list)
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'], outfields=['dwi_image','bvectors_directions','bvalues','T1_image']), name='datasource')
    datasource.inputs.base_directory = files_directory
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(dwi_image='%s/DWI.nii',
                                            bvalues='%s/b_values.txt',
                                            bvectors_directions='%s/b_vectors.txt',
                                            T1_image='%s/T1.nii')
    datasource.inputs.template_args = dict(dwi_image=[['subject_id']],
                                           bvalues=[['subject_id']],
                                           bvectors_directions=[['subject_id']],
                                           T1_image=[['subject_id']])
    datasource.inputs.sort_filelist = True
    
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=["dwi_image", "bvectors_directions", "bvalues", 'T1_image']), name="inputnode")
    
    datasink_prep = pe.Node(nio.DataSink(), name='datasink_prep')
    datasink_prep.inputs.base_directory = op.join(datasink_directory, 'pre_preprocess/')

    datasink_tracto = pe.Node(nio.DataSink(), name='datasink_tracto')
    datasink_tracto.inputs.base_directory = op.join(datasink_directory, 'Outputs_for_Tractography/')

    # Pre-preprocess' nodes definition
    
    b0_dwi_split = pe.Node(niu.Function(input_names=['in_file', 'in_bvals', 'in_bvecs'], output_names=['out_b0', 'out_dwi', 'out_bvals', 'out_bvecs'], function=b0_dwi_split), name='b0_dwi_split')
    b0_flirt = b0_flirt_pipeline(name='b0_co_registration')            
    b0_avg = pe.Node(niu.Function(input_names=['in_file'], output_names=['out_file'], function=b0_average), name='b0_average')
    mask_b0 = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True), name='mask_b0')
    insert_b0_into_dwi = pe.Node(niu.Function(input_names=['in_b0', 'in_dwi', 'in_bvals', 'in_bvecs'], output_names=['out_dwi', 'out_bvals', 'out_bvecs'], function=insert_b0_into_dwi), name='insert_b0avg_into_dwi')

    # Higher level pipeline construction
    
    hmc = hmc_pipeline()
    ecc = ecc_pipeline()
    epi = epi_pipeline()
    bias = remove_bias()
    
    aac = apply_all_corrections()
    
    preprocess = pe.Workflow(name='preprocess')
    preprocess.base_dir = working_directory
    
    preprocess.connect([(infosource, datasource, [('subject_id','subject_id')])])
    preprocess.connect([(datasource, inputnode, [('dwi_image','dwi_image'), ('bvalues','bvalues'), ('bvectors_directions','bvectors_directions'), ('T1_image','T1_image')])])
    preprocess.connect([(inputnode, b0_dwi_split, [('bvalues', 'in_bvals'),('bvectors_directions', 'in_bvecs'), ('dwi_image', 'in_file')])])
    preprocess.connect([(b0_dwi_split, b0_flirt, [('out_b0', 'inputnode.in_file')])])
    preprocess.connect([(b0_flirt, b0_avg, [('outputnode.out_file', 'in_file')])])
    preprocess.connect([(b0_avg, insert_b0_into_dwi, [('out_file', 'in_b0')])])
    preprocess.connect([(b0_avg, mask_b0,[('out_file', 'in_file')])])
    preprocess.connect([(b0_dwi_split, insert_b0_into_dwi, [('out_dwi','in_dwi'), ('out_bvals','in_bvals'), ('out_bvecs','in_bvecs')])])
    preprocess.connect([(insert_b0_into_dwi, datasink_prep, [('out_dwi','dwi_b0_merge'), ('out_bvals','out_bvals_dwi_b0_merge'),('out_bvecs','out_bvecs_dwi_b0_merge')])])
    preprocess.connect([(mask_b0, datasink_prep, [('mask_file','mask_b0')])])
    preprocess.connect([(b0_avg, datasink_prep, [('out_file','b0_average')])])
    
    preprocess.connect([(insert_b0_into_dwi, hmc,[('out_dwi','inputnode.in_file'), ('out_bvals','inputnode.in_bval'), ('out_bvecs','inputnode.in_bvec')])])
    preprocess.connect([(mask_b0, hmc, [('mask_file','inputnode.in_mask')])])
    preprocess.connect([(hmc, ecc, [('outputnode.out_xfms','inputnode.in_xfms'),('outputnode.out_file','inputnode.in_file')])])
    preprocess.connect([(insert_b0_into_dwi, ecc, [('out_bvals','inputnode.in_bval')])])
    preprocess.connect([(mask_b0, ecc, [('mask_file','inputnode.in_mask')])])
    preprocess.connect([(ecc, epi, [('outputnode.out_file','inputnode.DWI')])])
    preprocess.connect([(inputnode, epi, [('T1_image','inputnode.T1')])])
    preprocess.connect([(hmc, epi, [('outputnode.out_bvec','inputnode.bvec')])])
    
    preprocess.connect([(insert_b0_into_dwi, aac, [('out_dwi', 'inputnode.in_dwi')])])
    preprocess.connect([(hmc, aac, [('outputnode.out_xfms', 'inputnode.in_hmc')])])
    preprocess.connect([(ecc, aac, [('outputnode.out_xfms', 'inputnode.in_ecc')])])
    preprocess.connect([(epi, aac, [('outputnode.out_warp', 'inputnode.in_epi')])])
    preprocess.connect([(inputnode, aac, [('T1_image','inputnode.T1')])])
    
    
    preprocess.connect([(aac, bias, [('outputnode.out_file','inputnode.in_file')])])
    
    preprocess.connect([(aac, datasink_prep, [('outputnode.out_file','DWI_hmc_ecc_epi_corrected'),
                                                                ('outputnode.out_jacobian','out_jacobian'),
                                                                ('outputnode.out_coeff','out_coeff'),
                                                                ('outputnode.out_warp','out_warp')])])
    
    preprocess.connect([(bias, datasink_prep, [('outputnode.out_file','DWI_hmc_ecc_epi_bias_corrected')])])
    
    preprocess.connect([(bias, datasink_tracto, [('outputnode.out_file','DWI_hmc_ecc_epi_bias_corrected')])])
    preprocess.connect([(epi, datasink_tracto, [('outputnode.out_bvec', 'out_bvecs')])])
    preprocess.connect([(insert_b0_into_dwi, datasink_tracto, [('out_bvals','out_bvals')])])
    preprocess.connect([(bias, datasink_tracto, [('outputnode.b0_mask','b0_mask')])])
