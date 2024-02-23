def atlas_statistics(in_image, atlas_list):
    """Generate regional measure from atlas_list in TSV files.

    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Args:
        in_image: A Nifti image
        atlas_list: List of names of atlas to be applied

    Returns:
        List of paths to TSV files
    """
    from pathlib import Path

    from nipype.utils.filemanip import split_filename

    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.statistics import statistics_on_atlas
    from clinica.utils.ux import print_end_image

    subject_id = get_subject_id(in_image)

    orig_dir, base, ext = split_filename(in_image)
    atlas_statistics_list = []
    for atlas in atlas_list:
        out_atlas_statistics = Path(
            f"./{base}_space-{atlas}_map-graymatter_statistics.tsv"
        ).resolve()
        statistics_on_atlas(in_image, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)
    print_end_image(subject_id)
    return atlas_statistics_list
