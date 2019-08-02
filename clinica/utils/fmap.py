# coding: utf8


def rads2hz(in_file, delta_te, out_file=None):
    """
    Converts input phase difference map to Hz
    """
    import numpy as np
    import nibabel as nb
    import os.path as op
    import math
    from nipype.utils import NUMPY_MMAP

    if out_file is None:
        fname, fext = op.splitext(op.basename(in_file))
        if fext == '.gz':
            fname, _ = op.splitext(fname)
        out_file = op.abspath('./%s_radsec.nii.gz' % fname)

    im = nb.load(in_file, mmap=NUMPY_MMAP)
    data = im.get_data().astype(np.float32) * (1.0 / (delta_te * 2 * math.pi))
    nb.Nifti1Image(data, im.affine, im.header).to_filename(out_file)
    return out_file


def resample_fmap_to_b0(in_fmap, in_b0, out_file=None):
    """
    Resample fieldmap onto the b0 image.

    Warnings:
        The fieldmap should already be aligned on the b0.

    Args:
        in_fmap(str): Fieldmap image.
        in_b0(str): B0 image.
        out_file(optional[str]): Filename (default: <in_fmap>_space-b0.nii.gz

    Returns:
        Fieldmap image resampled on b0.
    """
    import os.path as op
    import nibabel
    from nilearn.image import resample_to_img

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_fmap))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_resampled_fmap = op.abspath("%s_space-b0%s" % (fname, ext))
    else:
        out_resampled_fmap = out_file

    resampled_fmap = resample_to_img(
        source_img=in_fmap, target_img=in_b0,
        interpolation='continuous')

    nibabel.nifti1.save(resampled_fmap, out_resampled_fmap)

    return out_resampled_fmap


def remove_filename_extension(in_file):
    """
    Remove extension from filename
    """
    from nipype.utils.filemanip import split_filename
    import os

    path, source_file_dwi, _ = split_filename(in_file)
    file_without_extension = os.path.join(path, source_file_dwi)

    return file_without_extension
