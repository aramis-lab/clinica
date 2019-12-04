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
    return [[f] * len(tissues) for f in flow_fields]


def prepare_existing_dartel_flowfields(flow_fields, tissues):
    return [f * len(tissues) for f in flow_fields]


def join_smoothed_files(smoothed_normalized_files):
    """
    Joins outputs
    """
    return [[x for smooth in subject for x in smooth] for subject in zip(*smoothed_normalized_files)]


def select_gm_images(in_images):
    """
    Selects only
    """
    return [image for subject in in_images for image in subject if ('c1' in image or 'graymatter' in image)]
