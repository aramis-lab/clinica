# coding: utf8

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


def atlas_statistics(in_image, atlas_list):
    """
    For each atlas name provided it calculates for the input image the mean for each region in the atlas and saves it to a TSV file.

    Args:
        in_image: A Nifti image
        atlas_list: List of names of atlas to be applied

    Returns:
        List of paths to TSV files
    """
    from os.path import abspath, join
    from nipype.utils.filemanip import split_filename
    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas
    from clinica.utils.stream import cprint

    orig_dir, base, ext = split_filename(in_image)
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics_list = []
    for atlas in atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(
                    join('./' + base + '_space-' + atlas + '_map-graymatter_statistics.tsv'))
                cprint(out_atlas_statistics)
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
    return atlas_statistics_list
