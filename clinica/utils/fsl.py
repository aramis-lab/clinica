# coding: utf8


"""This module contains FSL utilities."""


def standard_space_roi(in_t1, in_mask=None):
    """
    Pre-mask the structural image to standard space.

    The purpose of the FSL command standard_space_roi is to do an initial
    masking of the image before performing a brain extraction with FSL BET.
    This uses flirt to register the input image to a standard space whole-head
    image; the resulting transform is then inverted and a standard-space brain
    mask is transformed into the space of the input image, and then applied to
    this to create the output.

    Args:
        in_t1 (str): File containing the T1-weighted image.
        in_mask (Optional[str]): Mask output using transformed standard space
            mask (default=None, i.e. it will use the 2mm dilated MNI mask
            located in ${FSLDIR}/data/standard).

    Returns:
        out_pre_mask (str): Pre-mask the structural image.
    """
    import os.path as op
    import os

    assert(op.isfile(in_t1))

    out_pre_mask = op.abspath('T1_pre_bet.nii.gz')

    cmd = 'standard_space_roi %s %s -roiNONE' % (in_t1, out_pre_mask)
    if in_mask is None:
        cmd = cmd + ' -b'
    else:
        cmd = cmd + ' -maskMASK ' + in_mask
    os.system(cmd)
    return out_pre_mask
