from typing import List


def peak_correction(t_map: str, t_threshold: float, output_name: str = None) -> str:
    """Threshold the t_map with t_threshold.

    Pixel intensities that are less than t_threshold are set to 0, other values
    are left unchanged.

    Parameters
    ----------
    t_map: str
        Path to t-statistics nifti map
    t_threshold: float
        Threshold on t value
    output_name: str, optional
        Optional output name

    Returns
    -------
    str:
        Path to the generated file.
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


def cluster_correction(
    t_map: str, t_thresh: float, c_thresh: int, output_name: str = None
) -> str:
    """Performs cluster correction.

    First t_map is thresholded with t_thresh (like in peak_correction()). Then, clusters
    that have a size less than c_thresh are removed

    Parameters
    ----------
    t_map: str
        Path to t-statistics nifti map
    t_thresh: float
        Threshold on t value
    c_thresh: int
        Minimal size of clusters after thresholding
    output_name: str, optional
        Optional output name

    Returns
    -------
    str:
        Path to the generated file.
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


def produce_figures(
    nii_file: str,
    template: str,
    type_of_correction: str,
    t_thresh: str,
    c_thresh: int,
    n_cuts: int,
) -> list:
    """Produce the output figures.

    Parameters
    ----------
    nii_file: str
        Path to the nifti file (generated at previous steps)
    template: str
        Path to template used for the stat map plot
    type_of_correction: str
        Can be either FWE or FDR (used only in potential figure titles)
    t_thresh: str
        T value threshold used (used only in potential figure titles)
    c_thresh: int
        Cluster minimal size used (used only in potential figure titles)
    n_cuts: int
        Number of cuts in fig

    Returns
    -------
    list of string:
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


def generate_output(t_map: str, figs: list, correction_name: str) -> None:
    """Extract the output generated and copy it to the output folder

    Parameters
    ----------
    t_map: str
        Path to t-map on which whole pipeline was based
    figs: list of str
        Paths to figs to save
    correction_name: str
        Name of the correction (ex: cluster_correction_FWE)
    """
    from os import makedirs
    from os.path import basename, dirname, join, splitext
    from shutil import copyfile

    # Will extract group-GroupTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics from TStatistics file
    t_map_basename = splitext(basename(t_map))[0]

    out_folder = join(
        dirname(t_map), t_map_basename.replace("TStatistics", correction_name)
    )
    makedirs(out_folder)
    suffixes = (
        "GlassBrain.png",
        "axis-x_TStatistics.png",
        "axis-y_TStatistics.png",
        "axis-z_TStatistics.png",
    )
    for fig, suffix in zip(figs, suffixes):
        copyfile(
            fig,
            join(
                out_folder,
                t_map_basename.replace(
                    "TStatistics", f"desc-{correction_name}_{suffix}"
                ),
            ),
        )
