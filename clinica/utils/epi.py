# coding: utf-8


def bids_dir_to_fsl_dir(bids_dir):
    """
    Converts BIDS PhaseEncodingDirection parameters (i,j,k,i-,j-,k-) to
    FSL direction (x,y,z,x-,y-,z-).
    """
    fsl_dir = bids_dir.lower()
    if fsl_dir == 'i-':
        return 'x-'
    if fsl_dir == 'i':
        return 'x'
    if fsl_dir == 'j-':
        return 'y-'
    if fsl_dir == 'j':
        return 'y'
    if fsl_dir == 'k-':
        return 'z-'
    if fsl_dir == 'k':
        return 'z'

    raise RuntimeError("PhaseEncodingDirection " + fsl_dir + " is unknown, it should be a value in (x,y,z,x-,y-,z-)")

    return fsl_dir
