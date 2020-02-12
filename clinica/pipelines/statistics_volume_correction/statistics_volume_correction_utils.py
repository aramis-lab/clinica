# coding: utf8

"""
Statistics_Volume_Correction - Clinica Utilities.
"""

def peak_correction(t_map, t_threshold, output_name=None):
    import nibabel as nib
    from os.path import join, basename, abspath
    import numpy as np

    original_nifti = nib.load(t_map)
    data = original_nifti.get_data()
    data[data < t_threshold] = 0
    new_data = nib.Nifti1Image(data, affine=original_nifti.affine, header=original_nifti.header)
    if output_name:
        filename = output_name
    else:
        filename = join('./peak_corrected_' + str(t_threshold) + basename(t_map))
    nib.save(new_data, filename)
    return abspath(filename)


def cluster_correction(t_map, t_thresh, c_thresh, output_name=None):
    import nibabel as nib
    from os.path import join, basename, abspath
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
    return abspath(filename)


def produce_figures(nii_file, template, type_of_correction, t_thresh, c_thresh, n_cuts):
    from nilearn import plotting
    import numpy as np
    from os.path import abspath

    assert type_of_correction in ['FWE', 'FDR'], 'Type of correction must be FWE or FDR'
    if not np.isnan(c_thresh):
        correction = 'Cluster'
    else:
        correction = 'Peak'

    my_title = correction + ' correction ' + type_of_correction + ' Threshold = ' + str(t_thresh)
    if not np.isnan(c_thresh):
        my_title = my_title + ' - min cluster size = ' + str(c_thresh),

    plotting.plot_glass_brain(nii_file,
                              #title=my_title,
                              output_file='./glass_brain.png')

    plotting.plot_stat_map(nii_file,
                           display_mode='x',
                           cut_coords=np.linspace(-70, 67, n_cuts),
                           bg_img=template,
                           colorbar=False,
                           draw_cross=True,
                           output_file='./statmap_x.png')

    plotting.plot_stat_map(nii_file,
                           display_mode='y',
                           cut_coords=np.linspace(-104, 69, n_cuts),
                           bg_img=template,
                           colorbar=False,
                           draw_cross=True,
                           output_file='./statmap_y.png')

    plotting.plot_stat_map(nii_file,
                           display_mode='z',
                           cut_coords=np.linspace(-45, 78, n_cuts),
                           bg_img=template,
                           colorbar=False,
                           draw_cross=True,
                           output_file='./statmap_z.png')

    return [abspath('./glass_brain.png'), abspath('./statmap_x.png'), abspath('./statmap_y.png'), abspath('./statmap_z.png')]


def generate_output(t_map, figs, name):
    from os import makedirs
    from os.path import join, dirname, basename, splitext
    from shutil import copyfile

    outfolder = join(dirname(t_map), 'corrected_results_' + splitext(basename(t_map))[0], name)
    makedirs(outfolder)
    copyfile(figs[0], join(outfolder, 'glass_brain.png'))
    copyfile(figs[1], join(outfolder, 't_statistics_thresholded_x.png'))
    copyfile(figs[2], join(outfolder, 't_statistics_thresholded_y.png'))
    copyfile(figs[3], join(outfolder, 't_statistics_thresholded_z.png'))

