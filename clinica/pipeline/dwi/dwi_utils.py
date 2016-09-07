#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the post-processing pipeline."""

def dti_atlas_scalar_analysis(input_image, atlas_labels_image, name_output_file=None):
    """
    Compute statistics.

    Given a set a labeled tracts, this function compute statistics of a scalar image from DTI (e.g. FA, MD, etc.) in each tract.

    Args:
        input_image (str): File containing a scalar image (e.g. FA, MD, etc.).
        atlas_labels_image (str): File containing labels. These labels are used to compute statistics
        name_output_file (Optional[str]): Name of the output statistics file (default=scalar_stats.csv).

    Returns:
        outfile (str): CSV files containing the statistics (content of the columns: #Label, Mean, Standard deviation, #Voxels).
    """

    import nibabel as nib
    import numpy as np
    import pandas as pd
    import os.path as op

    if name_output_file is None:
        outfile = op.abspath('scalar_stats.csv')
    else:
        outfile = op.abspath(name_output_file);

    dti_atlas = nib.load(atlas_labels_image)
    atlas_image_data = dti_atlas.get_data()

    labels = list(set(atlas_image_data.ravel()))
    stats_scalar = np.zeros((len(labels),4))

    in_image = nib.load(input_image)
    scalar_image_data = in_image.get_data()

    for index, label in enumerate(labels):
        stats_scalar[index, 0] = label
        atlas_label_index = np.array(np.where(atlas_image_data==label))
        nb_voxel = atlas_label_index.shape[1]
        stats_scalar[index, 3] = nb_voxel
        labeled_voxel = labeled_voxel = scalar_image_data[atlas_label_index[0,:], atlas_label_index[1,:], atlas_label_index[2,:]]
        mean_scalar = labeled_voxel.mean()
        stats_scalar[index, 1] = mean_scalar
        std_scalar = labeled_voxel.std()
        stats_scalar[index, 2] = std_scalar

    columns = np.array(['Label', 'Mean scalar', 'Std scalar', 'Nb of voxel'])
    data = pd.DataFrame(stats_scalar, columns=columns)

    data.to_csv(outfile, sep=',', index=False)

    return outfile
