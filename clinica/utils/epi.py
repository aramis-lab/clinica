# coding: utf-8


def bids_dir_to_fsl_dir(bids_dir):
    """
    Converts BIDS PhaseEncodingDirection parameters (i,j,k,i-,j-,k-) to
    FSL direction (x,y,z,x-,y-,z-).
    """
    from clinica.utils.stream import cprint

    fsl_dir = bids_dir.lower()
    if fsl_dir == 'i-':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes x- for FSL Fugue")
        return 'x-'
    if fsl_dir == 'i':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes x for FSL Fugue")
        return 'x'
    if fsl_dir == 'j-':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes y- for FSL Fugue")
        return 'y-'
    if fsl_dir == 'j':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes y for FSL Fugue")
        return 'y'
    if fsl_dir == 'k-':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes z- for FSL Fugue")
        return 'z-'
    if fsl_dir == 'k':
        cprint("PhaseEncodingDirection " + fsl_dir + " becomes z for FSL Fugue")
        return

    raise RuntimeError("PhaseEncodingDirection " + fsl_dir + " is unknown, it should be a value in (x,y,z,x-,y-,z-)")

    return fsl_dir
