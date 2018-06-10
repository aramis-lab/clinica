import sys
import numpy as np
import os.path

def roundbval(bvalsFilename, newbvalsFilename, bStep):
    """round the bvals in the .bval file.

    If required, b-values can be rounded up to a specific threshold (bStep parameter).

    Parameters
    ----------
    :param str bvalsFilename: The path to bval file.
    :param str newbvalsFilename: The path to output scheme file (optional).
    :param float or list or np.bStep: If bStep is a scalar, round b-values to nearest integer multiple of bStep. If bStep is a list, it is treated as an array of shells in increasing order. B-values will be forced to the nearest shell value.
    """

    if not os.path.exists(bvalsFilename):
        raise RuntimeError( 'bvals file not exist:' + bvalsFilename )

    #if newbvalsFilename is None:
        #newbvalsFilename = os.path.splitext(bvalsFilename)[0]+".scheme"

    # load files and check size
    bvals = np.loadtxt(bvalsFilename)
    
    # convert bStep(str) to tuple to the np.array
    bStep = eval(bStep)

    # if requested, round the b-values
    bStep = np.array(bStep, dtype = np.float)
    print bStep 
    print bStep.size
    if bStep.size == 1 and bStep > 1.0:
        print "-> Rounding b-values to nearest multiple of %s" % np.array_str(bStep)
        bvals = np.round(bvals/bStep) * bStep
    elif bStep.size > 1:
        print "-> Setting b-values to the closest shell in %s" % np.array_str(bStep)
        for i in range(0, bvals.size):
            diff = min(abs(bvals[i] - bStep))
            ind = np.argmin(abs(bvals[i] - bStep))

            # warn if b > 99 is set to 0, possible error
            if (bStep[ind] == 0.0 and diff > 100) or (bStep[ind] > 0.0 and diff > bStep[ind] / 20.0):
                # For non-zero shells, warn if actual b-value is off by more than 5%. For zero shells, warn above 50. Assuming s / mm^2
                print "   Warning: measurement %d has b-value %d, being forced to %d\n'" % i, bvals[i], bStep[ind]

            bvals[i] = bStep[ind]

    np.savetxt( newbvalsFilename, bvals.T, fmt="%.06f", newline=' ' )
    print "-> Writing new rounded bval file to [ %s ]" % newbvalsFilename
    return newbvalsFilename


if __name__ == '__main__':
    bvalsFilename = sys.argv[1]
    newbvalsFilename = sys.argv[2]
    bStep = sys.argv[3] 
    roundbval( bvalsFilename, newbvalsFilename, bStep)