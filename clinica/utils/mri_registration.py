"""This module contains functions used for the registration aspects."""


def convert_flirt_transformation_to_mrtrix_transformation(
    in_source_image: str,
    in_reference_image: str,
    in_flirt_matrix: str,
    name_output_matrix=None,
) -> str:
    """Convert flirt matrix to mrtrix matrix.

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
        name_output_matrix (str, optional): Name of the output matrix.
            Defaults to "mrtrix_matrix.mat".

    Returns:
        out_mrtrix_matrix (str): Transformation matrix in MRtrix format.
    """
    import os

    from clinica.utils.check_dependency import check_mrtrix

    check_mrtrix()

    assert os.path.isfile(in_source_image)
    assert os.path.isfile(in_reference_image)
    assert os.path.isfile(in_flirt_matrix)

    out_mrtrix_matrix = os.path.abspath(name_output_matrix or "mrtrix_matrix.mat")

    cmd = f"transformconvert {in_flirt_matrix} {in_source_image} {in_reference_image} flirt_import {out_mrtrix_matrix}"
    os.system(cmd)

    return out_mrtrix_matrix
