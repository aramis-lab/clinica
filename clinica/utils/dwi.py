# coding: utf8


"""This module contains utilities for DWI handling."""


def merge_volumes_tdim(in_file1, in_file2):
    """
    Merge 'in_file1' and 'in_file2' in the t dimension.

    Args:
        in_file1 (str): First set of volumes.
        in_file2 (str): Second set of volumes.

    Returns:
        out_file (str): The two sets of volumes merged.
    """
    import os.path as op
    import os

    out_file = op.abspath('merged_files.nii.gz')
    cmd = 'fslmerge -t %s %s %s ' % (out_file, in_file1, in_file2)
    os.system(cmd)
    return out_file


def count_b0s(in_bval, low_bval=5.0):
    """
    Count the number of volumes where b<=low_bval.

    Args:
        in_bval (str): bval file.
        low_bval (Optional[int]): Define the b0 volumes as all volume
            bval <= lowbval. (Default=5.0)

    Returns:
        num_b0s: Number of b0s.
    """
    import numpy as np

    bvals = np.loadtxt(in_bval)
    num_b0s = len(np.where(bvals <= low_bval)[0])

    return num_b0s


def b0_average(in_file, out_file=None):
    """
    Average the b0 volumes.

    Args:
        in_file (str): The b0 volumes already registered.
        out_file (optional[str]): Name of the file.

    Returns:
        The mean of the b0 volumes.

    Warnings:
        The b0 volumes must be registered.
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


def b0_dwi_split(in_dwi, in_bval, in_bvec, low_bval=5.0):
    """
    Split the DWI volumes into two datasets :
     - the first dataset contains the set of b<=low_bval volumes.
     - the second dataset contains the set of DWI volumes.

    Args:
        in_dwi (str): DWI dataset.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.
        low_bval (Optional[int]): Define the b0 volumes as all volume
            bval <= lowbval. (Default=5.0)

    Returns:
        out_b0 (str): The set of b<=low_bval volumes.
        out_dwi (str): Output. The set of b>low_bval volumes.
        out_bvals (str): The b-values corresponding to the out_dwi.
        out_bvecs (str): The b-vecs corresponding to the out_dwi.
    """
    import numpy as np
    import nibabel as nib
    import os.path as op
    import warnings

    assert(op.isfile(in_dwi))
    assert(op.isfile(in_bval))
    assert(op.isfile(in_bvec))
    assert(low_bval >= 0)

    im = nib.load(in_dwi)
    data = im.get_data()
    hdr = im.get_header().copy()
    bvals = np.loadtxt(in_bval)
    bvecs = np.loadtxt(in_bvec)

    if bvals.shape[0] == bvecs.shape[0]:
        warnings.warn('Warning: The b-vectors file should be column-wise. '
                      + 'The b-vectors will be transposed',
                      UserWarning)
        bvecs = bvecs.T

    lowbs = np.where(bvals <= low_bval)[0]

    fname_b0, ext_b0 = op.splitext(op.basename(in_dwi))
    if ext_b0 == ".gz":
        fname_b0, ext2 = op.splitext(fname_b0)
        ext_b0 = ext2 + ext_b0
    out_b0 = op.abspath("%s_b0%s" % (fname_b0, ext_b0))
    # out_b0 = op.abspath('b0.nii.gz')
    b0data = data[..., lowbs]
    hdr.set_data_shape(b0data.shape)
    nib.Nifti1Image(b0data, im.get_affine(), hdr).to_filename(out_b0)

    dwi_bvals = np.where(bvals > low_bval)[0]
    out_dwi = op.abspath('dwi.nii.gz')
    dwi_data = data[..., dwi_bvals]
    hdr.set_data_shape(dwi_data.shape)
    nib.Nifti1Image(dwi_data, im.get_affine(), hdr).to_filename(out_dwi)

    bvals_dwi = bvals[dwi_bvals]
    out_bvals = op.abspath('bvals')
    np.savetxt(out_bvals, bvals_dwi, fmt='%d', delimiter=' ')

    bvecs_dwi = np.array([bvecs[0][dwi_bvals].tolist(),
                          bvecs[1][dwi_bvals].tolist(),
                          bvecs[2][dwi_bvals].tolist()])
    out_bvecs = op.abspath('bvecs')
    np.savetxt(out_bvecs, bvecs_dwi, fmt='%10.5f', delimiter=' ')

    return out_b0, out_dwi, out_bvals, out_bvecs


def insert_b0_into_dwi(in_b0, in_dwi, in_bval, in_bvec):
    """
    This function inserts a b0 volume into the dwi dataset as the first volume
    and updates the bvals and bvecs files.

    Args:
        in_b0 (str): One b=0 volume (could be the average of a b0 dataset).
        in_dwi (str): The set of DWI volumes.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.

    Returns:
        out_dwi (str): Diffusion dataset : b0 volume + dwi volumes.
        out_bval (str): B-values update.
        out_bvec (str): Directions of diffusion update.
    """
    from clinica.utils.dwi import merge_volumes_tdim
    import os.path as op
    import numpy as np

    assert(op.isfile(in_b0))
    assert(op.isfile(in_dwi))
    assert(op.isfile(in_bval))
    assert(op.isfile(in_bvec))

    out_dwi = merge_volumes_tdim(in_b0, in_dwi)

    lst = np.loadtxt(in_bval).tolist()
    lst.insert(0, 0)
    out_bvals = op.abspath('bvals')
    np.savetxt(out_bvals, np.matrix(lst), fmt='%d', delimiter=' ')

    bvecs = np.loadtxt(in_bvec)
    bvecs_0 = bvecs[0].tolist()
    bvecs_0.insert(0, 0.0)
    bvecs_1 = bvecs[1].tolist()
    bvecs_1.insert(0, 0.0)
    bvecs_2 = bvecs[2].tolist()
    bvecs_2.insert(0, 0.0)
    bvecs_dwi = np.array([bvecs_0, bvecs_1, bvecs_2])
    out_bvecs = op.abspath('bvecs')
    np.savetxt(out_bvecs, bvecs_dwi, fmt='%10.5f', delimiter=' ')

    return out_dwi, out_bvals, out_bvecs


def prepare_reference_b0(in_dwi, in_bval, in_bvec, low_bval=5):
    """
    This function prepare the data for further corrections. It coregisters
    the B0 images and then average it in order to obtain only
    one average B0 images.

    Args:
        in_dwi (str): Input DWI file.
        in_bvec (str): Vector file of the diffusion directions of
            the dwi dataset.
        in_bval (str): B-values file.
        low_bval (optional[int]):

    Returns:
        out_reference_b0 (str): Average of the B0 images or the only B0 image.
        out_b0_mask (str): Binary mask obtained from the average of
            the B0 images.
        out_b0_dwi_merge (str): Average of B0 images merged to the DWIs.
        out_updated_bval (str): Updated gradient values table.
        out_updated_bvec (str): Updated gradient vectors table.
    """
    from clinica.utils.dwi import (insert_b0_into_dwi, b0_dwi_split,
                                   count_b0s, b0_average)
    from clinica.workflows.dwi_preprocessing import b0_flirt_pipeline
    from clinica.utils.stream import cprint

    import os.path as op

    import tempfile

    # Count the number of b0s
    nb_b0s = count_b0s(in_bval=in_bval, low_bval=low_bval)
    # cprint('Found %s b0 for %s' % (nb_b0s, in_dwi))

    # Split dataset into two datasets: the b0 and the b>low_bval datasets
    [extracted_b0, out_split_dwi, out_split_bval, out_split_bvec] = \
        b0_dwi_split(
            in_dwi=in_dwi, in_bval=in_bval, in_bvec=in_bvec, low_bval=low_bval)

    if nb_b0s == 1:
        # The reference b0 is the extracted b0
        # cprint('Only one b0 for %s' % in_dwi)
        out_reference_b0 = extracted_b0
    elif nb_b0s > 1:
        # Register the b0 onto the first b0
        b0_flirt = b0_flirt_pipeline(num_b0s=nb_b0s)
        b0_flirt.inputs.inputnode.in_file = extracted_b0
        # BUG: Nipype does allow to extract the output after running the
        # workflow: we need to 'guess' where the output will be generated
        tmp_dir = tempfile.mkdtemp()
        b0_flirt.base_dir = tmp_dir
        b0_flirt.run()
        # out_node = b0_flirt.get_node('outputnode')
        registered_b0s = op.abspath(op.join(
            tmp_dir, 'b0_coregistration', 'concat_ref_moving',
            'merged_files.nii.gz'))
        # cprint('B0 s will be averaged (file = ' + registered_b0s + ')')
        # Average the b0s to obtain the reference b0
        out_reference_b0 = b0_average(in_file=registered_b0s)
    else:
        raise ValueError(
            'The number of b0s should be strictly positive (b-val file =%s).'
            % in_bval)

    # Merge datasets such that bval(DWI) = (0 b1 ... bn)
    [out_b0_dwi_merge, out_updated_bval, out_updated_bvec] = \
        insert_b0_into_dwi(in_b0=out_reference_b0, in_dwi=out_split_dwi,
                           in_bval=out_split_bval, in_bvec=out_split_bvec)

    return out_reference_b0, out_b0_dwi_merge, out_updated_bval, out_updated_bvec


def hmc_split(in_file, in_bval, ref_num=0, lowbval=5.0):
    """
    Selects the reference ('out_ref') and moving ('out_mov') volumes
    from a dwi dataset for the purpose of head motion correction (HMC).

    Args:
        in_file (str): DWI dataset.
        in_bval (str): File describing the b-values of the DWI dataset.
        ref_num (Optional[str]): The reference volume in the dwi dataset.
            Default ref_num= 0.
        lowbval (Optional[int]): Define the b0 volumes as all volume
            bval <= lowbval. (Default=5.0)
    Returns:
        out_ref (str): The reference volume.
        out_mov (str): The moving volume to align to the reference volume.
        out_bval (str): The b-values corresponding to the moving volume.
        volid (int): Index of the reference volume.
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


def check_dwi_volume(in_dwi, in_bvec, in_bval):
    """
    Check that # DWI = # B-val = # B-vec.

    Raises
        IOError
    """
    import nibabel as nib
    import numpy as np

    bvals = np.loadtxt(in_bval)
    num_b_vals = len(bvals)

    bvecs = np.loadtxt(in_bvec)
    _, num_b_vecs = bvecs.shape

    img = nib.load(in_dwi)
    _, _, _, num_dwis = img.shape

    if not (num_b_vals == num_b_vecs == num_dwis):
        raise IOError('Number of DWIs, b-vals and b-vecs mismatch '
                      '(# DWI = %s, # B-vec = %s, #B-val = %s) ' %
                      (num_dwis, num_b_vecs, num_b_vals))
