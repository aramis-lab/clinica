# coding: utf8


"""This module contains functions used for the registration aspects."""


def convert_flirt_transformation_to_mrtrix_transformation(
        in_source_image,
        in_reference_image,
        in_flirt_matrix,
        name_output_matrix=None):
    """
    Convert flirt matrix to mrtrix matrix.

    This function converts a transformation matrix produced by FSL's flirt
    command into a format usable by MRtrix. The output of this function
    is usually for the mrtransform command.

    Args:
        in_source_image (str): File containing the source image used in
            FSL flirt with the -in flag.
        in_reference_image (str): File containing the reference image used in
            FSL flirt with the -ref flag.
        in_flirt_matrix (str): File containing the transformation matrix
            obtained by FSL flirt.
        name_output_matrix (Optional[str]): Name of the output matrix
            (default=deformed_image.nii.gz).

    Returns:
        out_mrtrix_matrix (str): Transformation matrix in MRtrix format.
    """
    import os

    assert(os.path.isfile(in_source_image))
    assert(os.path.isfile(in_reference_image))
    assert(os.path.isfile(in_flirt_matrix))

    if name_output_matrix is None:
        out_mrtrix_matrix = os.path.abspath('mrtrix_matrix.mat')
    else:
        out_mrtrix_matrix = os.path.abspath(name_output_matrix)

    cmd = 'transformconvert %s %s %s flirt_import %s' \
          % (in_flirt_matrix, in_source_image, in_reference_image,
             out_mrtrix_matrix)
    os.system(cmd)

    return out_mrtrix_matrix


def apply_mrtrix_transform_without_resampling(
        in_image, in_mrtrix_matrix, name_output_image=None):
    """
    Apply a transformation without resampling.

    This function applies a linear transform on the input image without
    reslicing: it only modifies the transform matrix in the image header.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_mrtrix_matrix (str): File containing the transformation matrix
            obtained by the MRtrix transformconvert command.
        name_output_image (Optional[str]): Name of the output image
            (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according
            to in_mrtrix_matrix transformation.
    """
    import os

    assert(os.path.isfile(in_image))
    assert(os.path.isfile(in_mrtrix_matrix))

    if name_output_image is None:
        out_deformed_image = os.path.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = os.path.abspath(name_output_image)

    cmd = 'mrtransform -linear %s %s %s' \
          % (in_mrtrix_matrix, in_image, out_deformed_image)
    os.system(cmd)

    return out_deformed_image


def apply_ants_registration_syn_quick_transformation(
        in_image,
        in_reference_image,
        in_affine_transformation,
        in_bspline_transformation,
        name_output_image=None):
    """
    Apply a transformation obtained with antsRegistrationSyNQuick.sh.

    This function applies a rigid & deformable B-Spline syn transformation
    which has been estimated previously with antsRegistrationSyNQuick script.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_reference_image (str): File defining the spacing, origin, size,
            and direction of the output warped image.
        in_affine_transformation (str): File containing the transformation
            matrix obtained by antsRegistrationSyNQuick (expected file:
            [Prefix]0GenericAffine.mat).
        in_bspline_transformation (str): File containing the transformation
            matrix obtained by antsRegistrationSyNQuick (expected file:
            [Prefix]1Warp.nii.gz).
        name_output_image (Optional[str]): Name of the output image
            (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according
            to in_affine_transformation and in_bspline_transformation
            transformations.
    """
    import os

    assert(os.path.isfile(in_image))
    assert(os.path.isfile(in_affine_transformation))
    assert(os.path.isfile(in_bspline_transformation))

    if name_output_image is None:
        out_deformed_image = os.path.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = os.path.abspath(name_output_image)

    cmd = 'antsApplyTransforms -d 3 -e 0 -i %s -o %s -t %s -t %s -r %s ' \
          '--interpolation Linear' \
          % (in_image, out_deformed_image, in_bspline_transformation,
             in_affine_transformation, in_reference_image)
    os.system(cmd)

    return out_deformed_image


def ants_registration_syn_quick(fix_image, moving_image, prefix_output=None):
    """
    Small wrapper for antsRegistrationSyNQuick.sh.

    This function calls antsRegistrationSyNQuick.sh in order to register
    non-linearly moving_image towards fix_image.

    Args:
        fix_image (str): The target image.
        moving_image (str): The source
        prefix_output (str): Prefix for output files
            (format: <prefix_output>[Warped|0GenericAffine|1Warp|
            InverseWarped|1InverseWarp])

    Returns:
        The deformed image with the deformation parameters.
    """
    import os

    if prefix_output is None:
        prefix_output = 'SyN_Quick'

    image_warped = os.path.abspath(prefix_output + 'Warped.nii.gz')
    affine_matrix = os.path.abspath(prefix_output + '0GenericAffine.mat')
    warp = os.path.abspath(prefix_output + '1Warp.nii.gz')
    inverse_warped = os.path.abspath(prefix_output + 'InverseWarped.nii.gz')
    inverse_warp = os.path.abspath(prefix_output + '1InverseWarp.nii.gz')

    cmd = 'antsRegistrationSyNQuick.sh -t b -d 3 -f %s -m %s -o %s' \
          % (fix_image, moving_image, prefix_output)
    os.system(cmd)

    return image_warped, affine_matrix, warp, inverse_warped, inverse_warp


def ants_combine_transform(in_file, transforms_list, reference):
    """
    Apply a transformation obtained with antsRegistrationSyNQuick.sh.

    This function applies a rigid & deformable B-Spline syn transformation
    which has been estimated previously with antsRegistrationSyNQuick script.

    Args:
        in_file (str): File containing the input image to be transformed.
        reference (str): File defining the spacing, origin, size, and
            direction of the output warped image.
        transforms_list (str): File containing the transformations
            obtained by antsRegistrationSyNQuick.

    Returns:
        out_warp (str): File containing the deformed image according to
            transforms_list in_affine_transformation and
            in_bspline_transformation transformations.
    """
    import os

    out_warp = os.path.abspath('out_warp.nii.gz')

    transforms = ""
    for trans in transforms_list:
        transforms += " " + trans
    cmd = 'antsApplyTransforms -o [out_warp.nii.gz,1] -i %s -r %s -t %s' % \
          (in_file, reference, transforms)
    os.system(cmd)

    return out_warp
