# coding: utf8

"""Statistics_Volume_Correction - Clinica Utilities.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""

def peak_correction(t_map, t_threshold, output_name=None):
    import nibabel as nib
    from os.path import join, basename
    import numpy as np

    original_nifti = nib.load(t_map)
    data = original_nifti.get_data()
    data[data < t_threshold] = 0
    new_data = nib.Nifti1Image(data, affine=original_nifti.affine, header=original_nifti.header)
    if output_name:
        filename = output_name
    else:
        filename = join('./peak_corrected_' + str(t_threshold) + basename(t_map)
    nib.save(new_data, filename)
    return filename

def cluster_correction(t_map, t_thresh, c_thresh, output_name=None):
    import nibabel as nib
    from os.path import join, basename
    import numpy as np
    from scipy.ndimage.measurements import label

    original_nifti = nib.load(t_map)
    data = original_nifti.get_data()
    data[data < t_thresh] = 0
    labeled_mask, num_features = label(data)
    for i in range(1, num_features + 1):
        if np.sum(labeled_mask == i) < c_thresh:
            print('Label number ' + str(i) + ' cluster size is: ' + str(np.sum(labeled_mask == i)) + ' so it is removed')
            data[labeled_mask == i] = 0
    new_data = nib.Nifti1Image(data, affine=original_nifti.affine, header=original_nifti.header)
    if output_name:
        filename = output_name
    else:
        filename = join('./cluster_corrected_t-' + str(t_thresh) + '_c-' + str(c_thresh) + basename(t_map))
    nib.save(new_data, filename)
    return filename