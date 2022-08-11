from os import PathLike
from typing import Union

from nibabel.nifti1 import Nifti1Header
from numpy import ndarray


def initialize_tissues_spm_segment(
    tissue_probability_map: Union[PathLike, None] = None
) -> tuple:

    """Prepare data structure for SPM segment interface

    Parameters
    ----------
    tissue_probability_map: Union(PathLike, None)
        Path to the niftii tissue probability map

    Returns
    -------
    tuple
        "tissues" data structure for SPMSegment
    """

    import clinica.pydra.t1_volume.t1_volume_utils as spm_utils

    parameters = {}
    parameters.setdefault("tissue_classes", [1, 2, 3])
    parameters.setdefault("dartel_tissues", [1, 2, 3])
    parameters.setdefault("save_warped_unmodulated", True)
    parameters.setdefault("save_warped_modulated", False)

    if tissue_probability_map:
        parameters.setdefault("tissue_probability_maps", tissue_probability_map)
    else:
        parameters.setdefault("tissue_probability_maps", spm_utils.get_tpm())

    tissue_tuples = get_tissue_tuples(
        parameters["tissue_probability_maps"],
        parameters["tissue_classes"],
        parameters["dartel_tissues"],
        parameters["save_warped_unmodulated"],
        parameters["save_warped_modulated"],
    )
    return tissue_tuples


def init_input_node(bids_name: str) -> None:
    """Extract "sub-<participant_id>_ses-<session_label>" from <bids_name> and print begin message.

    Parameters
    ----------
    bids_name : str
        The BIDS name of the t1w image

    Warns
    -----
    Outputs ID of subjet/session being processed
    """

    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    subject_id = get_subject_id(bids_name)
    print_begin_image(subject_id)

    return subject_id, bids_name


def print_end_pipeline(bids_name: str) -> None:
    """Prints end message.

    Parameters
    ----------
    bids_name : str
        The BIDS name of the t1w image

    Warns
    -----
    Prints out message for image <bids_name>
    """
    from clinica.utils.ux import print_end_image

    print_end_image(bids_name)


def zip_list_files(class_images: list, zip_files: bool = False) -> list:
    """Create list of (Optionally) zipped nifti images
    Parameters
    ----------
    class_images : list
    zip_files: bool

    Returns
    -------
    list
        list of (optionally) zipped nifti files
    """
    from clinica.utils.filemanip import zip_nii

    if zip_files:
        return [zip_nii(tissue, True) for tissue in class_images]

    return [tissue for tissue in class_images]


def get_tissue_tuples(
    tissue_map: PathLike,
    tissue_classes: list,
    dartel_tissues: list,
    save_warped_unmodulated: bool,
    save_warped_modulated: bool,
) -> list:
    """Extract list of tuples, one for each tissue clas.

    Parameters
    ---------
    tissue_map : PathLike
        Path to tissue maps
    tissue_classes: list
        Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    dartel_tissues: list
        Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'
    save_warped_unmodulated : bool
        Save warped unmodulated images for tissues specified in --tissue_classes
    save_warped_modulated: bool
        Save warped modulated images for tissues specified in --tissue_classes

    Returns
    -------
    list
        List of tuples according to NewSegment input por tissues

    Notes
    -----
     The returned list contains tissue probability map (4D), 1-based index to frame
     - number of gaussians
     - which maps to save [Native, DARTEL] - a tuple of two boolean values
     - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values

    """
    tissues = []

    for i in range(1, 7):
        n_gaussians = 2

        if i == 4 or i == 5:
            n_gaussians = i - 1

        native_space = False
        dartel_input = False
        warped_unmodulated = False
        warped_modulated = False

        if i in tissue_classes:
            native_space = True
            if save_warped_unmodulated:
                warped_unmodulated = True
            if save_warped_modulated:
                warped_modulated = True

        if i in dartel_tissues:
            dartel_input = True

        tissues.append(
            (
                (tissue_map, i),
                n_gaussians,
                (native_space, dartel_input),
                (warped_unmodulated, warped_modulated),
            )
        )

    return tissues


def is_centered(nii_volume: PathLike, threshold_l2: int = 50) -> bool:
    """Checks if a NIfTI volume is centered on the origin of the world coordinate system.

    Parameters
    ---------
    nii_volume : PathLike
        path to NIfTI volume threshold_l2: maximum distance between origin of the world coordinate system and the center of the volume to
        be considered centered. The threshold were SPM segmentation stops working is around 100 mm (it was determined empirically after several
        trials on a generated dataset), so default value is 50mm in order to have a security margin, even when dealing with co-registered files afterward.

    Returns
    -------
        bool

    Notes
    -----------

    SPM has troubles to segment files if the center of the volume is not close from the origin of the world coordinate
    system. A series of experiment have been conducted: we take a volume whose center is on the origin of the world
    coordinate system. We add an offset using coordinates of affine matrix [0, 3], [1, 3], [2, 3] (or by modifying the
    header['srow_x'][3], header['srow_y'][3], header['srow_z'][3], this is strictly equivalent).

    It has been determined that volumes were still segmented with SPM when the L2 distance between origin and center of
    the volume did not exceed 100 mm. Above this distance, either the volume is either not segmented (SPM error), or the
    produced segmentation is wrong (not the shape of a brain anymore)
    """

    import numpy as np

    center = get_world_coordinate_of_center(nii_volume)

    distance_from_origin = np.linalg.norm(center, ord=2)

    if distance_from_origin < threshold_l2:
        return True
    else:
        # If center is a np.nan,
        return False


def get_world_coordinate_of_center(nii_volume: PathLike) -> ndarray:
    """Extract the world coordinates of the center of the image.

    Parameters
    ---------
    nii_volume : PathLike
        path to nii volume

    Returns
    -------
    tuple
        coordinates in the world space

    References
    ------
    https://brainder.org/2012/09/23/the-nifti-file-format/

    """

    from os.path import isfile

    import nibabel as nib
    import numpy as np
    from nibabel.filebasedimages import ImageFileError

    # from clinica.utils.stream import cprint
    # assert isinstance(nii_volume, str), "input argument nii_volume must be a str"

    assert isfile(nii_volume), "input argument must be a path to a file"

    try:
        orig_nifti = nib.load(str(nii_volume))
    except ImageFileError:
        print(f"File {nii_volume} could not be read by nibabel. Is it a valid NIfTI")
        return np.nan

    head = orig_nifti.header

    if isinstance(head, nib.freesurfer.mghformat.MGHHeader):
        # If MGH volume
        center_coordinates_world = vox_to_world_space_method_3_bis(
            head["dims"][0:3] / 2, head
        )
    else:
        # Standard NIfTI volume
        center_coordinates = get_center_volume(head)

        if head["qform_code"] > 0:
            center_coordinates_world = vox_to_world_space_method_2(
                center_coordinates, head
            )
        elif head["sform_code"] > 0:
            center_coordinates_world = vox_to_world_space_method_3(
                center_coordinates, head
            )
        elif head["sform_code"] == 0:
            center_coordinates_world = vox_to_world_space_method_1(
                center_coordinates, head
            )
        else:
            center_coordinates_world = np.nan
    return center_coordinates_world


def get_center_volume(header: dict) -> ndarray:
    """Get the voxel coordinates of the center of the data, using header information.

    Parameters
    ----------
        header: Nifti1Header
            Contains image metadata

    Returns
    -------
    ndarray
        Voxel coordinates of the center of the volume
    """
    import numpy as np

    center_x = header["dim"][1] / 2
    center_y = header["dim"][2] / 2
    center_z = header["dim"][3] / 2
    return np.array([center_x, center_y, center_z])


def vox_to_world_space_method_1(
    coordinates_vol: ndarray, header: Nifti1Header
) -> ndarray:
    """Convert coordinates to world space

    Parameters
    ----------
    coordinates_vol: ndarray
        Coordinate in the volume (raw data)
    header: Nifti1Header
        Contains image metadata

    Returns
    -------
    ndarray
        Coordinates in the world space

    Notes
    -----
    This method is for compatibility with analyze and is not supposed to be used as the main orientation method. But it
    is used if sform_code = 0. The world coordinates are determined simply by scaling by the voxel size by their
    dimension stored in pixdim.

    References
    ----------
    https://brainder.org/2012/09/23/the-nifti-file-format/
    """
    import numpy as np

    return np.array(coordinates_vol) * np.array(
        header["pixdim"][1], header["pixdim"][2], header["pixdim"][3]
    )


def vox_to_world_space_method_2(
    coordinates_vol: ndarray, header: Nifti1Header
) -> ndarray:
    """Convert coordinates to world space (method 2)


    Parameters
    ----------
    coordinates_vol : ndarray
        Coordinate in the volume (raw data)
    header : Nifti1Header
        Image header containing metadata

    Returns
    -------
    ndarray
        Coordinates in the world space

    Notes
    -----
    This method is used when short qform_code is larger than zero. To get the coordinates, we multiply a rotation
    matrix (r_mat) by coordinates_vol, then perform Hadamard with pixel dimension pixdim (like in method 1). Then we add
    an offset (qoffset_x, qoffset_y, qoffset_z)
    """
    import numpy as np

    def get_r_matrix(h):
        """Get rotation matrix.

        More information here: https://brainder.org/2012/09/23/the-nifti-file-format/

        Parameters
        ----------
            h: header

        Returns
        -------
            Rotation matrix
        """
        b = h["quatern_b"]
        c = h["quatern_c"]
        d = h["quatern_d"]
        a = np.sqrt(1 - (b**2) - (c**2) - (d**2))
        r = np.zeros((3, 3))
        r[0, 0] = (a**2) + (b**2) - (c**2) - (d**2)
        r[0, 1] = 2 * ((b * c) - (a * d))
        r[0, 2] = 2 * ((b * d) + (a * c))
        r[1, 0] = 2 * ((b * c) + (a * d))
        r[1, 1] = (a**2) + (c**2) - (b**2) - (d**2)
        r[1, 2] = 2 * ((c * d) - (a * b))
        r[2, 0] = 2 * ((b * d) - (a * c))
        r[2, 1] = 2 * ((b * d) - (a * c))
        r[2, 2] = (a**2) + (d**2) - (b**2) - (c**2)
        return r

    i = coordinates_vol[0]
    j = coordinates_vol[1]
    k = coordinates_vol[2]
    if header["qform_code"] > 0:
        r_mat = get_r_matrix(header)
    else:
        # Should never be reached
        raise ValueError("qform_code must be greater than 0 to use this method")
    q = header["pixdim"][0]
    if q not in [-1, 1]:
        print("q was " + str(q), ", now is 1")
        q = 1
    return np.dot(r_mat, np.array([i, j, q * k])) * np.array(
        header["pixdim"][1:4]
    ) + np.array([header["qoffset_x"], header["qoffset_y"], header["qoffset_z"]])


def vox_to_world_space_method_3(coordinates_vol: ndarray, header: Nifti1Header):
    """Convert coordinates to world space (method 3)


    Parameters
    ----------
    coordinates_vol : ndarray
        Coordinate in the volume (raw data)
    header : Nifti1Header
        Image header containing metadata

    Returns
    -------
    ndarray
        Coordinates in the world space

    Notes
    -----
    This method is used when sform_code is larger than zero. It relies on a full affine matrix, stored in the header in
    the fields srow_[x,y,y], to map voxel to world coordinates.
    When a nifti file is created with raw data and affine=..., this is this method that is used to decipher the
    voxel-to-world correspondence.

    """
    import numpy as np

    def get_aff_matrix(h: Nifti1Header) -> ndarray:
        """Get affine transformation matrix.

        Parameters
        ----------
            h : Nifti1Header

        Returns
        -------
        ndarray
            Affine transformation matrix

        References
        ----------
        https://brainder.org/2012/09/23/the-nifti-file-format/

        """
        mat = np.zeros((4, 4))
        mat[0, 0] = h["srow_x"][0]
        mat[0, 1] = h["srow_x"][1]
        mat[0, 2] = h["srow_x"][2]
        mat[0, 3] = h["srow_x"][3]
        mat[1, 0] = h["srow_y"][0]
        mat[1, 1] = h["srow_y"][1]
        mat[1, 2] = h["srow_y"][2]
        mat[1, 3] = h["srow_y"][3]
        mat[2, 0] = h["srow_z"][0]
        mat[2, 1] = h["srow_z"][1]
        mat[2, 2] = h["srow_z"][2]
        mat[2, 3] = h["srow_z"][3]
        mat[3, 3] = 1
        return mat

    if header["sform_code"] > 0:
        aff = get_aff_matrix(header)
    else:
        # Should never be reached
        raise ValueError("sform_code has a value > 0, so method 3 cannot be used")

    homogeneous_coord = np.concatenate(
        (np.array(coordinates_vol), np.array([1])), axis=0
    )
    return np.dot(aff, homogeneous_coord)[0:3]


def vox_to_world_space_method_3_bis(coordinates_vol, header):
    """
    This method relies on the same technique as method 3, but for images created by FreeSurfer (MGHImage, MGHHeader).
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: nib.freesurfer.mghformat.MGHHeader object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    affine_transformation_matrix = header.get_affine()
    homogeneous_coord = np.concatenate(
        (np.array(coordinates_vol), np.array([1])), axis=0
    )
    return np.dot(affine_transformation_matrix, homogeneous_coord)[0:3]


def get_tpm() -> PathLike:
    """Get Tissue Probability Map (TPM) from SPM.

    Returns
    -------
    str
        TPM.nii path from SPM
    """
    import os
    from glob import glob
    from os.path import join

    spm_home = os.getenv("SPM_HOME")

    if not spm_home:
        spm_home = os.getenv("SPMSTANDALONE_HOME")

    if not spm_home:
        raise RuntimeError(
            "Could not determine location of your SPM installation. Neither $SPM_HOME "
            "or $SPMSTANDALONE_HOME are present in your environment"
        )

    tpm_file_glob = glob(join(spm_home, "**/TPM.nii"), recursive=True)

    if len(tpm_file_glob) == 0:
        raise RuntimeError(f"No file found for TPM.nii in your $SPM_HOME in {spm_home}")

    if len(tpm_file_glob) > 1:
        error_str = f"Multiple files found for TPM.nii in your SPM_HOME {spm_home}:"

        for file in tpm_file_glob:
            error_str += "\n\t" + file

        raise RuntimeError(error_str)

    return tpm_file_glob[0]
