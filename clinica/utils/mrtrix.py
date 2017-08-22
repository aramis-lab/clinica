"""This module contains MRtrix utilities."""


def dilate_mask(in_mask, npass=4, nthreads=2):
    """
    Dilate the mask.

    This function performs the dilation of a given binary (brain) mask.

    Args:
        in_mask (str): Binary mask.
        npass (Optional[int]): Number of times to repeatedly apply the dilation
            (default=4).
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 0 disables multi-threading).

    Returns:
        out_dilated_mask : FILE
          Output. Dilated mask.
    """
    import os.path as op
    import os

    assert(op.isfile(in_mask))

    out_dilated_mask = op.abspath('dilated_mask.nii')
    cmd = 'maskfilter -npass ' + str(npass) + ' -nthreads ' + str(nthreads) + ' ' + in_mask + ' dilate ' + out_dilated_mask
    os.system(cmd)
    return out_dilated_mask
