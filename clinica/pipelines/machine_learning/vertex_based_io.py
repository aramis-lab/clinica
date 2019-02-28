# coding: utf8

import numpy as np
import pandas as pd
import nibabel as nib
import os

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@inria.fr"
__status__ = "Development"


def load_data(mgh_list):
    """

    Args: mgh_list : list of mgh files. Each element contains as many paths as
    needed (each element must be associated to a single subject). Surfaces must
    have the same number of vertices accross subjects.

    Returns: data : matrix of raw data

    """
    import nibabel as nib
    import numpy as np

    # Construct 0-matrix with the good size, based on the size of the surfaces
    # provided by the first subject
    N_vertex = []  # array containing the surface size of the different surfaces of a subject
    sample = mgh_list[0]
    for i in range(len(sample)):
        N_vertex.append(np.max(nib.load(sample[i]).header.get_data_shape()))
    data = np.zeros((len(mgh_list), np.sum(N_vertex)))

    # Fill data matrix
    N_cumul = np.concatenate(([0], np.cumsum(N_vertex)))
    for s in range(len(mgh_list)):
        for h in range(len(mgh_list[s])):
            current_data = np.squeeze(nib.load(mgh_list[s][h]).get_data())
            data[s, N_cumul[h]: N_cumul[h+1]] = current_data
    return data
