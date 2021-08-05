# coding: utf8
from typing import Optional


def convert_phase_in_radians(in_file: str, out_file: Optional[str] = None):
    """Convert phase image in radians."""
    import math
    import os
    import os.path as op

    import nibabel as nb
    import numpy as np

    assert op.isfile(in_file)

    img = nb.load(in_file)
    imin = np.amin(img.get_fdata(dtype="float32"))
    imax = np.amax(img.get_fdata(dtype="float32"))

    out_file = out_file or op.abspath("phase_in_rad.nii.gz")

    cmd = "fslmaths %s -mul %s -div %s %s -odt float" % (
        in_file,
        2.0 * math.pi,
        imax,
        out_file,
    )
    os.system(cmd)

    data = (
        img2.get_fdata(dtype="float32").astype(np.float32)
        - img1.get_fdata(dtype="float32").astype(np.float32)
    ) * (1.0 / delta_te)
    nb.Nifti1Image(data, img.get_affine(), img.get_header()).to_filename(out_file)

    return out_file


def create_phase_in_radsec(
    in_phase1, in_phase2, delta_te, out_file: Optional[str] = None
):
    """Converts input (unwarpped) phase1 and phase2 map to into a fieldmap inrads.

    Warning:
        delta_te should be in seconds.
    """
    import os.path as op

    import nibabel as nb
    import numpy as np

    out_file = out_file or op.abspath("fmap_radsec.nii.gz")

    img1 = nb.load(in_phase1)
    img2 = nb.load(in_phase2)
    data = (
        img2.get_fdata(dtype="float32").astype(np.float32)
        - img1.get_fdata(dtype="float32").astype(np.float32)
    ) * (1.0 / delta_te)
    nb.Nifti1Image(data, img1.get_affine(), img1.get_header()).to_filename(out_file)
    return out_file


def resample_fmap_to_b0(in_fmap, in_b0, out_file: Optional[str] = None):
    """Resample fieldmap onto the b0 image.

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

    if not out_file:
        fname, ext = op.splitext(op.basename(in_fmap))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_resampled_fmap = op.abspath(f"%s_space-b0%s" % (fname, ext))
    else:
        out_resampled_fmap = out_file

    resampled_fmap = resample_to_img(
        source_img=in_fmap, target_img=in_b0, interpolation="continuous"
    )

    nibabel.nifti1.save(resampled_fmap, out_resampled_fmap)

    return out_resampled_fmap
