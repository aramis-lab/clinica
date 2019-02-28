# coding: utf8

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def read_psf(json_path):
    """readpsf is able to read the json file associated with the PET volume, and extract information needed to perform
    partial volume correction : the effective resolution in plane, and the effective resolution axial.

        Args:
            (string) json_path : Path to the json file containing the information

        Returns:
            (float) The effective resolution in plane
            (float) The effective axial resolution
        """
    import json
    import os

    if not os.path.exists(json_path):
        raise Exception('File ' + json_path + ' does not exist. Point spread function information must be provided.')

    # This is a standard way of reading data in a json
    with open(json_path) as df:
        data = json.load(df)
    in_plane = data['Psf'][0]['EffectiveResolutionInPlane']
    axial = data['Psf'][0]['EffectiveResolutionAxial']
    return in_plane, axial


def create_binary_mask(tissues, threshold=0.3):
    """

    Args:
        tissues:
        threshold:

    Returns:

    """
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join, basename

    if len(tissues) == 0:
        raise RuntimeError('The length of the list of tissues must be greater than zero.')

    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_data().shape)

    data = np.zeros(shape=shape)

    for image in tissues:
        data = data + nib.load(image).get_data()

    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + '_brainmask.nii')

    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def apply_binary_mask(image, binary_mask):
    import nibabel as nib
    from os import getcwd
    from os.path import join, basename

    original_image = nib.load(image)
    mask = nib.load(binary_mask)

    data = original_image.get_data() * mask.get_data()

    masked_image_path = join(getcwd(), 'masked_' + basename(image))
    masked_image = nib.Nifti1Image(data, original_image.affine, header=original_image.header)
    nib.save(masked_image, masked_image_path)
    return masked_image_path


def create_pvc_mask(tissues):

    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join

    if len(tissues) == 0:
        raise RuntimeError('The length of the list of tissues must be greater than zero.')

    img_0 = nib.load(tissues[0])
    shape = img_0.get_data().shape
    background = np.zeros(shape=shape)

    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)

    for i in range(len(tissues)):
        image = nib.load(tissues[i])
        data[..., i] = np.array(image.get_data())
        background = background + image.get_data()

    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)

    out_mask = join(getcwd(), 'pvc_mask.nii')
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def pet_pvc_name(pet_image, pvc_method):
    from os.path import basename
    pet_pvc_path = 'pvc-' + pvc_method.lower() + '_' + basename(pet_image)
    return pet_pvc_path


def normalize_to_reference(pet_image, region_mask):
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import basename, join

    pet = nib.load(pet_image)
    ref = nib.load(region_mask)

    region = np.multiply(pet.get_data(), ref.get_data())
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    data = pet.get_data() / region_mean

    suvr_pet_path = join(getcwd(), 'suvr_' + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def atlas_statistics(in_image, in_atlas_list):
    """

    Args:
        in_image:
        in_atlas_list:

    Returns:

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
                out_atlas_statistics = abspath(join(getcwd(), base + '_space-' + atlas + '_statistics.tsv'))
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break

    return atlas_statistics_list


def pet_container_from_filename(pet_filename):
    """

    Args:
        pet_filename:

    Returns:

    """
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', pet_filename)

    if m is None:
        raise ValueError('Input filename is not in a BIDS or CAPS compliant format. It doesn\'t contain the subject' +
                         ' and session informations.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session, 'pet/preprocessing')


def expand_into_list(in_field, n_tissues):
    return [in_field] * n_tissues


def get_from_list(in_list, index):
    return in_list[index]
