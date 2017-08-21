#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the registration pipeline or other pipelines."""


def convert_flirt_transformation_to_mrtrix_transformation(
        in_source_image, in_reference_image, in_flirt_matrix, name_output_matrix=None):
    """
    Convert flirt matrix to mrtrix matrix.

    This function converts a transformation matrix produced by FSL's flirt command into a format usable by MRtrix.
    The output of this function is usually for the mrtransform command.

    Args:
        in_source_image (str): File containing the source image used in FSL flirt with the -in flag.
        in_reference_image (str): File containing the reference image used in FSL flirt with the -ref flag.
        in_flirt_matrix (str): File containing the transformation matrix obtained by FSL flirt.
        name_output_matrix (Optional[str]): Name of the output matrix (default=deformed_image.nii.gz).

    Returns:
        out_mrtrix_matrix (str): Transformation matrix in MRtrix format.
    """
    import os.path as op
    import os

    assert(op.isfile(in_source_image))
    assert(op.isfile(in_reference_image))
    assert(op.isfile(in_flirt_matrix))

    if name_output_matrix is None:
        out_mrtrix_matrix = op.abspath('mrtrix_matrix.mat')
    else:
        out_mrtrix_matrix = op.abspath(name_output_matrix)

    cmd = 'transformconvert ' + in_flirt_matrix + ' ' + in_source_image + ' ' + in_reference_image + ' flirt_import ' + out_mrtrix_matrix
    os.system(cmd)

    return out_mrtrix_matrix


def apply_mrtrix_transform_without_resampling(in_image, in_mrtrix_matrix, name_output_image=None):
    """
    Apply a transformation without resampling.

    This function applies a linear transform on the input image without
    reslicing: it only modifies the transform matrix in the image header.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_mrtrix_matrix (str): File containing the transformation matrix obtained by the MRtrix transformconvert command.
        name_output_image (Optional[str]): Name of the output image (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according
            to in_mrtrix_matrix transformation.

    Example:
        >>> from clinica.utils.mri_registration import apply_mrtrix_transform_without_resampling
        >>> apply_mrtrix_transform_without_resampling('t1_image.nii', 't1_to_diffusion.txt', 't1_in_diffusion_space.nii')
    """
    import os.path as op
    import os

    assert(op.isfile(in_image))
    assert(op.isfile(in_mrtrix_matrix))

    if name_output_image is None:
        out_deformed_image = op.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = op.abspath(name_output_image)

    cmd = 'mrtransform -linear ' + in_mrtrix_matrix + ' ' + in_image + ' ' + out_deformed_image
    os.system(cmd)

    return out_deformed_image


def apply_ants_registration_syn_quick_transformation(
        in_image, in_reference_image, in_affine_transformation, in_bspline_transformation, name_output_image=None):
    """
    Apply a transformation obtained with antsRegistrationSyNQuick.sh.

    This function applies a rigid & deformable B-Spline syn transformation which has been estimated previously with
    antsRegistrationSyNQuick script.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_reference_image (str): File defining the spacing, origin, size, and direction of the output warped image.
        in_affine_transformation (str): File containing the transformation matrix obtained by antsRegistrationSyNQuick
            (expected file: [Prefix]0GenericAffine.mat).
        in_bspline_transformation (str): File containing the transformation matrix obtained by antsRegistrationSyNQuick
            (expected file: [Prefix]1Warp.nii.gz).
        name_output_image (Optional[str]): Name of the output image (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according to in_affine_transformation and
            in_bspline_transformation transformations.

    Example:
        >>> from clinica.utils.mri_registration import apply_ants_registration_syn_quick_transformation
        >>> apply_ants_registration_syn_quick_transformation('my_image.nii.gz', 'output0GenericAffine.mat', 'output1Warp.nii.gz', 'my_deformed_image.nii.gz')
    """
    import os.path as op
    import os

    assert(op.isfile(in_image))
    assert(op.isfile(in_affine_transformation))
    assert(op.isfile(in_bspline_transformation))

    if name_output_image is None:
        out_deformed_image = op.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = op.abspath(name_output_image)

    cmd = 'antsApplyTransforms -d 3 -e 0 -i ' + in_image + ' -o ' + out_deformed_image + ' -t ' + in_bspline_transformation + ' -t ' + in_affine_transformation + ' -r ' + in_reference_image + ' --interpolation BSpline'
    os.system(cmd)

    return out_deformed_image


def ants_registration_syn_quick(fixe_image, moving_image, prefix_output=None):
    """
    TODO.

    TODO

    Args:

    Returns:

    Example:
        >>> from clinica.utils.mri_registration import apply_ants_registration_syn_quick_transformation
        >>> apply_ants_registration_syn_quick_transformation('my_image.nii.gz', 'output0GenericAffine.mat', 'output1Warp.nii.gz', 'my_deformed_image.nii.gz')
    """
    import subprocess
    import os.path as op

    image_warped = op.abspath('SyN_QuickWarped.nii.gz')
    affine_matrix = op.abspath('SyN_Quick0GenericAffine.mat')
    warp = op.abspath('SyN_Quick1Warp.nii.gz')
    inverse_warped = op.abspath('SyN_QuickInverseWarped.nii.gz')
    inverse_warp = op.abspath('SyN_Quick1InverseWarp.nii.gz')

    cmd = 'antsRegistrationSyNQuick.sh -t br -d 3 -f ' + fixe_image + ' -m ' + moving_image + ' -o SyN_Quick'
    subprocess.call([cmd], shell=True)

    return image_warped, affine_matrix, warp, inverse_warped, inverse_warp


def ants_combine_transform(in_file, transforms_list, reference):
    """
    Apply a transformation obtained with antsRegistrationSyNQuick.sh.

    This function applies a rigid & deformable B-Spline syn transformation which has been estimated previously with
    antsRegistrationSyNQuick script.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_reference_image (str): File defining the spacing, origin, size, and direction of the output warped image.
        in_affine_transformation (str): File containing the affine transformation matrix obtained by antsRegistrationSyNQuick
            (expected file: [Prefix]0GenericAffine.mat).
        in_bspline_transformation (str): File containing the BSpline transformation obtained by antsRegistrationSyNQuick
            (expected file: [Prefix]1Warp.nii.gz).
        name_output_image (Optional[str]): Name of the output image (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according to in_affine_transformation and
            in_bspline_transformation transformations.

    Example:
        >>> from clinica.utils.mri_registration import apply_ants_registration_syn_quick_transformation
        >>> apply_ants_registration_syn_quick_transformation('my_image.nii.gz', 'output0GenericAffine.mat', 'output1Warp.nii.gz', 'my_deformed_image.nii.gz')
    """
    import os
    import os.path as op

    out_warp = op.abspath('out_warp.nii.gz')

    transforms = ""
    for trans in transforms_list:
        transforms += " " + trans
    cmd = 'antsApplyTransforms -o [out_warp.nii.gz,1] -i ' + in_file + ' -r ' + reference + ' -t' + transforms
    os.system(cmd)

    return out_warp
