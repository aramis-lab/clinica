def peak_correction(t_map, t_threshold, output_name=None):
    """
    Threshold the t_map with t_threshold. Pixel intensities that are less than t_threshold are set to 0, other values
    are left unchanged.

    Args:
        t_map: (str) path to t-statistics nifti map
        t_threshold: (float) threshold on t value
        output_name: (str) optional output name

    Returns:
        path to the generated file.
    """
    from os.path import abspath, basename, join

    import nibabel as nib

    original_nifti = nib.load(t_map)
    data = original_nifti.get_fdata(dtype="float32")
    data[data < t_threshold] = 0
    new_data = nib.Nifti1Image(
        data, affine=original_nifti.affine, header=original_nifti.header
    )
    if output_name:
        filename = output_name
    else:
        filename = join("./peak_corrected_" + str(t_threshold) + basename(t_map))
    nib.save(new_data, filename)
    return abspath(filename)


def cluster_correction(t_map, t_thresh, c_thresh, output_name=None):
    """
    Performs cluster correction. First t_map is thresholded with t_thresh (like in peak_correction()). Then, clusters
    that have a size less than c_thresh are removed
    Args:
        t_map: (str) path to t-statistics nifti map
        t_thresh: (float) threshold on t value
        c_thresh: (int) minimal size of clusters after thresholding
        output_name: (str) optional output name

    Returns:
        path to the generated file.
    """
    from os.path import abspath, basename, join

    import nibabel as nib
    import numpy as np
    from scipy.ndimage.measurements import label

    original_nifti = nib.load(t_map)
    data = original_nifti.get_fdata(dtype="float32")
    data[data < t_thresh] = 0
    labeled_mask, num_features = label(data)
    for i in range(1, num_features + 1):
        if np.sum(labeled_mask == i) < c_thresh:
            print(
                "Label number "
                + str(i)
                + " cluster size is: "
                + str(np.sum(labeled_mask == i))
                + " so it is removed"
            )
            data[labeled_mask == i] = 0
    new_data = nib.Nifti1Image(
        data, affine=original_nifti.affine, header=original_nifti.header
    )
    if output_name:
        filename = output_name
    else:
        filename = join(
            "./cluster_corrected_t-"
            + str(t_thresh)
            + "_c-"
            + str(c_thresh)
            + basename(t_map)
        )
    nib.save(new_data, filename)
    return abspath(filename)


def produce_figures(nii_file, template, type_of_correction, t_thresh, c_thresh, n_cuts):
    """
    Produce the output figures

    Args:
        nii_file: (str) path to the nifti file (generated at previous steps)
        template: (str) path to template used for the stat map plot
        type_of_correction: (str) Can be either FWE or FDR (used only in potential figure titles)
        t_thresh: (str) t value threshold used (used only in potential figure titles)
        c_thresh: (int) cluster minimal size used (used only in potential figure titles)
        n_cuts: (int) number of cuts in fig

    Returns:
        List of path to image files: glass brain, statmap along x, statmap along y, statmap along z
    """
    from os.path import abspath

    import numpy as np
    from nilearn import plotting

    assert type_of_correction in ["FWE", "FDR"], "Type of correction must be FWE or FDR"
    if not np.isnan(c_thresh):
        correction = "Cluster"
    else:
        correction = "Peak"

    my_title = (
        correction
        + " correction "
        + type_of_correction
        + " Threshold = "
        + str(t_thresh)
    )
    if not np.isnan(c_thresh):
        my_title = (my_title + " - min cluster size = " + str(c_thresh),)

    plotting.plot_glass_brain(nii_file, output_file="./glass_brain.png")

    plotting.plot_stat_map(
        nii_file,
        display_mode="x",
        cut_coords=np.linspace(-70, 67, n_cuts),
        bg_img=template,
        colorbar=False,
        draw_cross=True,
        output_file="./statmap_x.png",
    )

    plotting.plot_stat_map(
        nii_file,
        display_mode="y",
        cut_coords=np.linspace(-104, 69, n_cuts),
        bg_img=template,
        colorbar=False,
        draw_cross=True,
        output_file="./statmap_y.png",
    )

    plotting.plot_stat_map(
        nii_file,
        display_mode="z",
        cut_coords=np.linspace(-45, 78, n_cuts),
        bg_img=template,
        colorbar=False,
        draw_cross=True,
        output_file="./statmap_z.png",
    )

    return [
        abspath("./glass_brain.png"),
        abspath("./statmap_x.png"),
        abspath("./statmap_y.png"),
        abspath("./statmap_z.png"),
    ]


def generate_output(t_map, figs, name):
    """
        Produce output
    Args:
        t_map: (str) path to t-map on which whole pipeline was based
        figs: (list of str) paths to figs to save
        name: (str) name of the correction (ex: cluster_correction_FWE)

    Returns:
        Nothing
    """
    from os import makedirs
    from os.path import basename, dirname, join, splitext
    from shutil import copyfile

    # Will extract group-GroupTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics from TStatistics file
    t_map_basename = splitext(basename(t_map))[0]

    out_folder = join(dirname(t_map), t_map_basename.replace("TStatistics", name))
    makedirs(out_folder)
    copyfile(
        figs[0],
        join(
            out_folder,
            t_map_basename.replace("TStatistics", "desc-" + name + "_GlassBrain.png"),
        ),
    )
    copyfile(
        figs[1],
        join(
            out_folder,
            t_map_basename.replace(
                "TStatistics", "desc-" + name + "_axis-x_TStatistics.png"
            ),
        ),
    )
    copyfile(
        figs[2],
        join(
            out_folder,
            t_map_basename.replace(
                "TStatistics", "desc-" + name + "_axis-y_TStatistics.png"
            ),
        ),
    )
    copyfile(
        figs[3],
        join(
            out_folder,
            t_map_basename.replace(
                "TStatistics", "desc-" + name + "_axis-z_TStatistics.png"
            ),
        ),
    )
