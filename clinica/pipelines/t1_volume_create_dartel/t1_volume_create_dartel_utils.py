# coding: utf8

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def get_class_images(class_images, index_list):
    """
    Utility method to extract class images from a multi session <class_images> set.

    Tissue types are:
        - 1: Grey Matter
        - 2: White Matter
        - 3: CerebroSpinal Fluid
        - 4: Skull
        - 5: Out-of-brain soft tissue
        - 6: Head surrounding

    Args:
        class_images: image set from which to extract images.
        index_list: index list of the classes to extract.

    Returns:
        Extracted images in a list of lists (without empty lists).

    Example:
    >>> class_n_images = get_class_images(class_images,[1,2])
    """
    # Declare class images list
    class_n_images = {}
    for idx in index_list:
        class_n_images[idx] = []

    for session in class_images:
        for idx in index_list:
            class_n_images[idx].extend(session[idx-1])

    result = []
    for idx in index_list:
        if class_n_images[idx]:
            result.append(class_n_images[idx])

    return result
