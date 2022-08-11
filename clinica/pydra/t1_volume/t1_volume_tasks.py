from multiprocessing.dummy import Array
from os import PathLike
from pathlib import PosixPath
from typing import Union
from pydra.mark import annotate, task
from numpy import ndarray
from nibabel.nifti1 import Nifti1Header


@task
@annotate({"return": {"out_file": PathLike}})
def zip_nii(in_var: Union[PathLike, list], same_dir: bool = False) -> PathLike:
    """Zips a file or a list of files

    Parameters
    ---------
    in_var : Union[PathLike, list]
        file or list of files to zip

    same_dir : bool
        True if we want to zip in the origin directory, False if in cwd

    Returns
    ------
    PathLike
        path of the output

    """

    import gzip
    import shutil
    import os

    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_var is None) or isinstance(in_var, _Undefined):
        return None

    if not isinstance(in_var, str):  # type(in_file) is list:
        return [zip_nii(f, same_dir) for f in in_var]

    orig_dir, base, ext = split_filename(str(in_var))

    # Already compressed

    if os.path.splitext("ext") == ".gz":
        return in_var

    # Not compressed

    out_file = os.abspath(
        os.join(orig_dir if same_dir else os.getcwd(), base + ext + ".gz")
    )

    with open(in_var, "rb") as f_in, gzip.open(out_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


@task
@annotate({"return": {"out_status": None}})
def check_volume_location_in_world_coordinate_system(
    nifti_list: list,
    bids_dir: PathLike,
    modality: str = "t1w",
    skip_question: bool = False,
) -> None:
    """Check if images are centered around the origin of the world coordinate

    Parameters
    ----
    nifti_list: list
        list of path to nifti files

    bids_dir: str
        path to bids directory associated with this check

    modality: str
        the modality of the image

    skip_question: bool
        if True, assume answer is yes

    Returns
    -------
    None

    Warns
    ------
    If volume is not centered on origin of the world coordinate system

    Notes
    -----
    the NIfTI file list provided in argument are approximately centered around the origin of the
    world coordinates. Otherwise, issues may arise with further processing such as SPM segmentation. When not centered,
    we warn the user of the problem propose to exit clinica to run clinica iotools center-nifti or to continue with the execution
    of the pipeline

    """

    import sys
    from os.path import abspath, basename

    import click
    import numpy as np

    nifti_list = []
    if isinstance(nifti_list, PosixPath):
        if not is_centered(nifti_list):
            list_non_centered_files.append(nifti_list)
    elif isinstance(nifti_list, list):
        list_non_centered_files = [file for file in nifti_list if not is_centered(file)]

    if len(list_non_centered_files) > 0:
        centers = [
            get_world_coordinate_of_center(file) for file in list_non_centered_files
        ]
        l2_norm = [np.linalg.norm(center, ord=2) for center in centers]

        # File column width : 3 spaces more than the longest string to display
        file_width = 3 + max(len(basename(file)) for file in list_non_centered_files)

        # Center column width (with a fixed minimum size) : 3 spaces more than the longest string to display
        center_width = max(
            len("Coordinate of center") + 3,
            3 + max(len(str(center)) for center in centers),
        )

        warning_message = (
            f"It appears that {str(len(list_non_centered_files))} files "
            "have a center way out of the origin of the world coordinate system. SPM has a high "
            "probability to fail on these files (for co-registration or segmentation):\n\n"
        )
        warning_message += (
            "%-" + str(file_width) + "s%-" + str(center_width) + "s%-s"
        ) % ("File", "Coordinate of center", "Distance to origin")
        # 18 is the length of the string 'Distance to origin'
        warning_message += "\n" + "-" * (file_width + center_width + 18) + "\n"
        for file, center, l2 in zip(list_non_centered_files, centers, l2_norm):
            warning_message += (
                "%-" + str(file_width) + "s%-" + str(center_width) + "s%-25.2f\n"
            ) % (basename(file), str(center), l2)

        cmd_line = f"`clinica iotools center-nifti {abspath(bids_dir)} {abspath(bids_dir)}_centered --modality {modality}`"

        warning_message += (
            "\nIf you are trying to launch the t1-freesurfer pipeline, you can ignore this message "
            "if you do not want to run the pet-surface pipeline afterward."
        )

        warning_message += (
            "\nClinica provides a tool to counter this problem by replacing the center of the volume"
            " at the origin of the world coordinates.\nUse the following command line to correct the "
            f"header of the faulty NIFTI volumes in a new folder:\n{cmd_line}"
            "You will find more information on the command by typing "
            "clinica iotools center-nifti in the console."
        )

        click.echo(warning_message)

        if not skip_question:
            if not click.confirm("Do you still want to launch the pipeline?"):
                click.echo("Clinica will now exit...")
                sys.exit(0)
    return


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


def get_tpm():
    """Extracts Tissue Probability Map (TPM) from SPM.

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
        # Try MCR to get a hint on SPM location
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
