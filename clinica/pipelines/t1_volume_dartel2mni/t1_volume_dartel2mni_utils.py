# coding: utf8

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def prepare_flowfields(flow_fields, tissues):
    """

    :param flow_fields:
    :param tissues:
    :return:
    """
    return [[f] * len(tissues) for f in flow_fields]


def prepare_existing_dartel_flowfields(flow_fields, tissues):
    """

    :param flow_fields:
    :param tissues:
    :return:
    """
    return [f * len(tissues) for f in flow_fields]


def join_smoothed_files(smoothed_normalized_files):
    """
    Joins outputs
    :param smoothed_normalized_files:
    :return:
    """

    return [[x for smooth in subject for x in smooth] for subject in zip(*smoothed_normalized_files)]


def atlas_statistics(in_image, in_atlas_list):
    """
    For each atlas name provided it calculates for the input image the mean for each region in the atlas and saves it to a tsv file.
    :param in_image: A Nifti image
    :param in_atlas_list: List of names of atlas to be applied
    :return: List of paths to tsv files
    """
    from os import getcwd
    from os.path import abspath, join
    from nipype.utils.filemanip import split_filename
    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas

    orig_dir, base, ext = split_filename(in_image)
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics_list = []
    for atlas in in_atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(join(getcwd(), base + '_space-' + atlas
                                                    + '_map-graymatter_statistics.tsv'))
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break

    return atlas_statistics_list


def select_gm_images(in_images):
    """
    Selects only
    :param in_images:
    :return:
    """
    return [image for subject in in_images for image in subject if ('c1' in image or 'graymatter' in image)]
