from os import PathLike
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import nibabel as nib
import numpy as np
import pandas as pd

__all__ = [
    "center_nifti_origin",
    "check_volume_location_in_world_coordinate_system",
    "check_relative_volume_location_in_world_coordinate_system",
    "center_all_nifti",
]


def center_nifti_origin(input_image: PathLike, output_image: PathLike) -> PathLike:
    """Put the origin of the coordinate system at the center of the image.

    Parameters
    ----------
    input_image : PathLike
        Path to the input image.

    output_image : PathLike
        Path to the output image (where the result will be stored).

    Returns
    -------
    PathLike :
        The path of the output image created.
    """
    input_image = nib.load(Path(input_image))
    output_image = Path(output_image)
    canonical_image = nib.as_closest_canonical(input_image)
    header = canonical_image.header
    new_image = nib.Nifti1Image(
        canonical_image.get_fdata(caching="unchanged"),
        affine=_compute_qform(header),
        header=header,
    )
    # Without deleting already-existing file, nib.save causes a severe bug on Linux system
    if output_image.is_file():
        output_image.unlink()
    nib.save(new_image, output_image)
    if not output_image.is_file():
        raise RuntimeError(
            f"NIfTI file created but Clinica could not save it to {output_image}. "
            "Please check that the output folder has the correct permissions."
        )

    return output_image


def _compute_qform(header: nib.Nifti1Header) -> np.ndarray:
    qform = np.zeros((4, 4))
    for i in range(1, 4):
        qform[i - 1, i - 1] = header["pixdim"][i]
        qform[i - 1, 3] = -1.0 * header["pixdim"][i] * header["dim"][i] / 2.0
    return qform


def check_volume_location_in_world_coordinate_system(
    nifti_list: List[PathLike],
    bids_dir: PathLike,
    modality: str = "t1w",
    skip_question: bool = False,
) -> bool:
    """Check if images are centered around the origin of the world coordinate.

    Parameters
    ----------
    nifti_list : list of PathLike
        List of path to nifti files.

    bids_dir : PathLike
        Path to bids directory associated with this check.

    modality : str, optional
        The modality of the image. Default='t1w'.

    skip_question : bool, optional
        If True, assume answer is yes. Default=False.

    Returns
    -------
    bool :
        True if they are centered, False otherwise

    Warns
    ------
    If volume is not centered on origin of the world coordinate system

    Notes
    -----
    The NIfTI file list provided in argument are approximately centered around the origin of the world coordinates.
    Otherwise, issues may arise with further processing such as SPM segmentation.
    When not centered, we warn the user of the problem propose to exit clinica to run clinica iotools center-nifti
    or to continue with the execution of the pipeline.
    """
    import click
    import numpy as np

    bids_dir = Path(bids_dir)
    list_non_centered_files = [
        Path(file) for file in nifti_list if not _is_centered(Path(file))
    ]
    if len(list_non_centered_files) == 0:
        return True
    centers = [
        _get_world_coordinate_of_center(file) for file in list_non_centered_files
    ]
    l2_norm = [np.linalg.norm(center, ord=2) for center in centers]
    click.echo(
        _build_warning_message(
            list_non_centered_files, centers, l2_norm, bids_dir, modality
        )
    )
    if not skip_question:
        _ask_for_confirmation()

    return False


def _ask_for_confirmation() -> None:
    import sys

    import click

    if not click.confirm("Do you still want to launch the pipeline?"):
        click.echo("Clinica will now exit...")
        sys.exit(0)


def _build_warning_message(
    non_centered_files: List[Path],
    centers: List[np.array],
    l2_norm: List[np.ndarray],
    bids_dir: Path,
    modality: str,
) -> str:
    warning_message = (
        f"It appears that {len(non_centered_files)} files "
        "have a center way out of the origin of the world coordinate system. SPM has a high "
        "probability to fail on these files (for co-registration or segmentation):\n\n"
    )
    df = pd.DataFrame(
        [
            (a, b, c)
            for a, b, c in zip([f.name for f in non_centered_files], centers, l2_norm)
        ],
        columns=["File", "Coordinate of center", "Distance to origin"],
    )
    warning_message += str(df)
    cmd_line = f"`clinica iotools center-nifti {bids_dir.resolve()} {bids_dir.resolve()}_centered --modality {modality}`"
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

    return warning_message


def check_relative_volume_location_in_world_coordinate_system(
    label_1: str,
    nifti_list1: List[PathLike],
    label_2: str,
    nifti_list2: List[PathLike],
    bids_dir: PathLike,
    modality: str,
    skip_question: bool = False,
):
    """Check if the NIfTI file list `nifti_list1` and `nifti_list2` provided in argument are not too far apart,
    otherwise coreg in SPM may fail. Norm between center of volumes of 2 files must be less than 80 mm.

    Parameters
    ----------
    label_1 : str
        Label of the first nifti_list1 files (used in potential warning message).

    nifti_list1 : list of PathLike
        First list of files.

    label_2 : str
        Label of the second nifti_list.

    nifti_list2 : list of PathLike
        Second list of files, must be same length as `nifti_list1`.

    bids_dir : PathLike
        BIDS directory (used in potential warning message).

    modality : str
        String that must be used in argument of:
        clinica iotools bids --modality <MODALITY>
        (used in potential warning message).

    skip_question : bool, optional
        Disable prompts for user input (default is False)
    """
    import warnings

    from clinica.utils.stream import cprint

    warning_message = _build_warning_message_relative_volume_location(
        label_1, nifti_list1, label_2, nifti_list2, bids_dir, modality
    )
    cprint(msg=warning_message, lvl="warning")
    warnings.warn(warning_message)
    if not skip_question:
        _ask_for_confirmation()


def _build_warning_message_relative_volume_location(
    label_1: str,
    nifti_list1: List[PathLike],
    label_2: str,
    nifti_list2: List[PathLike],
    bids_dir: PathLike,
    modality: str,
) -> str:
    bids_dir = Path(bids_dir)
    file_couples = [(Path(f1), Path(f2)) for f1, f2 in zip(nifti_list1, nifti_list2)]
    df = pd.DataFrame(
        _get_problematic_pairs_with_l2_norm(file_couples),
        columns=[label_1, label_2, "Relative distance"],
    )
    if len(df) == 0:
        return
    warning_message = (
        f"It appears that {len(df)} pairs of files have an important relative offset. "
        "SPM co-registration has a high probability to fail on these files:\n\n"
    )
    warning_message += str(df)
    warning_message += (
        "\nClinica provides a tool to counter this problem by replacing the center "
        "of the volume at the origin of the world coordinates.\nUse the following "
        "command line to correct the header of the faulty NIFTI volumes in a new folder:\n\n"
        f"`clinica iotools center-nifti {bids_dir.resolve()} {bids_dir.resolve()}_centered --modality {modality}`\n\n"
        "You will find more information on the command by typing `clinica iotools center-nifti` in the console."
    )
    return warning_message


def _get_problematic_pairs_with_l2_norm(
    file_couples: List[Tuple[Path, Path]],
    threshold: float = 80.0,
) -> List[Tuple[str, str, float]]:
    l2_norm = _compute_l2_norm(file_couples)
    pairs_with_problems = [i for i, norm in enumerate(l2_norm) if norm > threshold]

    return [
        (file_couples[k][0].name, file_couples[k][1].name, l2_norm[k])
        for k in pairs_with_problems
    ]


def _compute_l2_norm(file_couples: List[Tuple[Path, Path]]) -> List[float]:
    center_coordinates = [
        (_get_world_coordinate_of_center(f[0]), _get_world_coordinate_of_center(f[1]))
        for f in file_couples
    ]

    return [np.linalg.norm(center[0] - center[1]) for center in center_coordinates]


def center_all_nifti(
    bids_dir: PathLike,
    output_dir: PathLike,
    modalities: Optional[Iterable[str]] = None,
    center_all_files: bool = False,
) -> List[Path]:
    """Center all the NIfTI images of the input BIDS folder into the empty output_dir specified in argument.

    All the files from bids_dir are copied into output_dir, then all the NIfTI images found are replaced by their
    centered version if their center is off the origin by more than 50 mm.

    Parameters
    ----------
    bids_dir : PathLike
        Path to the BIDS directory.

    output_dir : PathLike
        Path to the output directory where the centered files will be written to.

    modalities : iterable of str, optional
        Process these modalities only. Process all modalities otherwise.

    center_all_files:  bool, default=False
        Center files that may cause problem for SPM if set to False, all files otherwise.

    Returns
    -------
    list of Path
        Centered NIfTI files.
    """
    from shutil import copy, copytree

    from clinica.utils.stream import cprint

    bids_dir, output_dir = _validate_bids_and_output_dir(bids_dir, output_dir)
    for f in bids_dir.iterdir():
        if f.is_dir() and not (output_dir / f.name).is_dir():
            copytree(f, output_dir / f.name, copy_function=copy)
        elif f.is_file() and not (output_dir / f.name).is_file():
            copy(f, output_dir / f.name)
    nifti_files_filtered: List[Path] = []
    for f in output_dir.glob("**/*.nii*"):
        if modalities is None or any(
            elem.lower() in f.name.lower() for elem in modalities
        ):
            nifti_files_filtered.append(f)
    if not center_all_files:
        nifti_files_filtered = [
            file for file in nifti_files_filtered if not _is_centered(file)
        ]
    errors: List[str] = []
    for f in nifti_files_filtered:
        cprint(msg=f"Handling file {f}", lvl="debug")
        try:
            center_nifti_origin(f, f)
        except Exception as e:
            errors.append(str(e))
    if errors:
        raise RuntimeError(
            f"Clinica encountered {len(errors)} error(s) while trying to center all NIfTI images.\n"
            + "\n".join(errors)
        )
    return nifti_files_filtered


def _validate_bids_and_output_dir(
    bids_dir: PathLike, output_dir: PathLike
) -> Tuple[Path, Path]:
    from clinica.utils.exceptions import ClinicaBIDSError
    from clinica.utils.inputs import check_bids_folder

    bids_dir = Path(bids_dir)
    output_dir = Path(output_dir)
    if bids_dir == output_dir:
        raise ClinicaBIDSError(
            f"Input BIDS ({bids_dir}) and output ({output_dir}) directories must be different."
        )
    check_bids_folder(bids_dir)

    return bids_dir, output_dir


def _is_centered(nii_volume: Path, threshold_l2: int = 50) -> bool:
    """Checks if a NIfTI volume is centered on the origin of the world coordinate system.

    Parameters
    ---------
    nii_volume : Path
        The path to the NIfTI volume to check.

    threshold_l2: int, optional
        Maximum distance between the origin of the world coordinate system and the center of the volume to
        be considered centered. The threshold where SPM segmentation stops working is around 100 mm
        (it was determined empirically after several trials on a generated dataset), so default value
        is 50mm in order to have a security margin, even when dealing with co-registered files afterward.

    Returns
    -------
    bool :
        True if the volume is centered, False otherwise.

    Notes
    ------
    SPM has troubles to segment files if the center of the volume is not close from the origin of the world coordinate
    system. A series of experiment have been conducted: we take a volume whose center is on the origin of the world
    coordinate system. We add an offset using coordinates of affine matrix [0, 3], [1, 3], [2, 3] (or by modifying the
    header['srow_x'][3], header['srow_y'][3], header['srow_z'][3], this is strictly equivalent).
    It has been determined that volumes were still segmented with SPM when the L2 distance between origin and center of
    the volume did not exceed 100 mm. Above this distance, either the volume is either not segmented (SPM error), or the
    produced segmentation is wrong (not the shape of a brain anymore)
    """
    center = _get_world_coordinate_of_center(nii_volume)
    if center is None:
        raise ValueError(
            f"Unable to compute the world coordinates of center for image {nii_volume}."
            "Please verify the image data and header."
        )
    distance_from_origin = np.linalg.norm(center, ord=2)

    return distance_from_origin < threshold_l2


def _get_world_coordinate_of_center(nii_volume: Path) -> Optional[np.ndarray]:
    """Extract the world coordinates of the center of the image.

    Parameters
    ---------
    nii_volume : PathLike
        The path to the nii volume.

    Returns
    -------
    np.ndarray :
        The coordinates in the world space.

    References
    ------
    https://brainder.org/2012/09/23/the-nifti-file-format/
    """
    from nibabel.filebasedimages import ImageFileError

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.stream import cprint

    if not nii_volume.is_file():
        raise ClinicaException(
            f"The input {nii_volume} does not appear to be a path to a file."
        )
    try:
        orig_nifti = nib.load(nii_volume)
    except ImageFileError:
        cprint(
            msg=f"File {nii_volume} could not be read by nibabel. Is it a valid NIfTI file ?",
            lvl="warning",
        )
        return None

    header = orig_nifti.header
    if isinstance(header, nib.freesurfer.mghformat.MGHHeader):
        return _scale_with_affine_transformation_matrix(
            header["dims"][0:3] / 2, header, freesurfer_image=True
        )
    try:
        center_coordinates = _get_center_volume(header)
        if header["qform_code"] > 0:
            return _scale_with_rotation_matrix(center_coordinates, header)
        if header["sform_code"] > 0:
            return _scale_with_affine_transformation_matrix(center_coordinates, header)
        if header["sform_code"] == 0:
            return _scale_coordinates_by_pixdim(center_coordinates, header)
    except KeyError as e:
        raise ValueError(
            f"Cannot get the world coordinates of the center of image {nii_volume}. "
            f"This is most likely due to missing data in the header : {e}."
        )
    return None


def _get_center_volume(header: nib.Nifti1Header) -> np.ndarray:
    """Get the voxel coordinates of the center of the data, using header information.

    Parameters
    ----------
    header: Nifti1Header
        Image header which contains image metadata.

    Returns
    -------
    ndarray :
        Voxel coordinates of the center of the volume
    """
    return np.array([x / 2 for x in header["dim"][1:4]])


def _scale_with_affine_transformation_matrix(
    coordinates_vol: np.ndarray,
    header: nib.Nifti1Header,
    freesurfer_image: bool = False,
) -> np.ndarray:
    """Convert coordinates to world space using the affine transformation matrix.

    Parameters
    ----------
    coordinates_vol : ndarray
        Coordinate in the volume (raw data).

    header : Nifti1Header
        Image header containing metadata.
        The header must have sform_code > 0.

    freesurfer_image : bool, optional
        Whether the image for which the scaling is computed was obtained with
        Freesurfer (MGHImage) or not.
        Default=False.

    Returns
    -------
    ndarray
        Coordinates in the world space

    Notes
    -----
    This method is used when sform_code is larger than zero.
    It relies on a full affine matrix, stored in the header in the fields srow_[x,y,y],
    to map voxel to world coordinates. When a nifti file is created with raw data and affine=...,
    this is this method that is used to decipher the voxel-to-world correspondence.
    """
    homogeneous_coord = np.concatenate(
        (np.array(coordinates_vol), np.array([1])), axis=0
    )
    affine = (
        header.get_affine()
        if freesurfer_image
        else _get_affine_transformation_matrix(header)
    )

    return np.dot(affine, homogeneous_coord)[0:3]


def _get_affine_transformation_matrix(header: nib.Nifti1Header) -> np.ndarray:
    """Get affine transformation matrix.

    Parameters
    ----------
    header : Nifti1Header
        The image header for which to compute the affine transformation matrix.

    Returns
    -------
    ndarray :
        The computed affine transformation matrix

    References
    ----------
    https://brainder.org/2012/09/23/the-nifti-file-format/
    """
    matrix = np.zeros((4, 4))
    matrix[0, 0] = header["srow_x"][0]
    matrix[0, 1] = header["srow_x"][1]
    matrix[0, 2] = header["srow_x"][2]
    matrix[0, 3] = header["srow_x"][3]
    matrix[1, 0] = header["srow_y"][0]
    matrix[1, 1] = header["srow_y"][1]
    matrix[1, 2] = header["srow_y"][2]
    matrix[1, 3] = header["srow_y"][3]
    matrix[2, 0] = header["srow_z"][0]
    matrix[2, 1] = header["srow_z"][1]
    matrix[2, 2] = header["srow_z"][2]
    matrix[2, 3] = header["srow_z"][3]
    matrix[3, 3] = 1

    return matrix


def _scale_with_rotation_matrix(
    coordinates_vol: np.ndarray,
    header: nib.Nifti1Header,
) -> np.ndarray:
    """Convert coordinates to world space using the rotation matrix.

    Parameters
    ----------
    coordinates_vol : ndarray
        Coordinates in the volume (raw data).

    header : Nifti1Header
        Image header containing metadata.
        The header must have a qform_code > 0.

    Returns
    -------
    ndarray :
        Coordinates in the world space.

    Notes
    -----
    This method is used when short qform_code is larger than zero.
    To get the coordinates, we multiply a rotation matrix (r_mat) by coordinates_vol,
    then perform Hadamard with pixel dimension pixdim (like in method 1).
    Then we add an offset (qoffset_x, qoffset_y, qoffset_z)
    """
    from clinica.utils.stream import cprint

    q = header["pixdim"][0]
    offset = np.array([header[f"qoffset_{x}"] for x in ("x", "y", "z")])
    if q not in (-1, 1):
        cprint(
            f"Pixdim of provided header was {q} while either -1 or 1 was expected. Using 1.",
            lvl="warning",
        )
        q = 1
    return (
        np.dot(
            _get_rotation_matrix(header),
            np.array([coordinates_vol[0], coordinates_vol[1], q * coordinates_vol[2]]),
        )
        * np.array(header["pixdim"][1:4])
        + offset
    )


def _get_rotation_matrix(header: nib.Nifti1Header) -> np.ndarray:
    """Get the rotation matrix from the provided image header.

    More information here: https://brainder.org/2012/09/23/the-nifti-file-format/

    Parameters
    ----------
    header : Nifti1Header
        The header

    Returns
    -------
    np.ndarray :
        The rotation matrix.
    """
    b = header["quatern_b"]
    c = header["quatern_c"]
    d = header["quatern_d"]
    a = np.sqrt(1 - (b**2) - (c**2) - (d**2))
    rotation = np.zeros((3, 3))
    rotation[0, 0] = (a**2) + (b**2) - (c**2) - (d**2)
    rotation[0, 1] = 2 * ((b * c) - (a * d))
    rotation[0, 2] = 2 * ((b * d) + (a * c))
    rotation[1, 0] = 2 * ((b * c) + (a * d))
    rotation[1, 1] = (a**2) + (c**2) - (b**2) - (d**2)
    rotation[1, 2] = 2 * ((c * d) - (a * b))
    rotation[2, 0] = 2 * ((b * d) - (a * c))
    rotation[2, 1] = 2 * ((b * d) - (a * c))
    rotation[2, 2] = (a**2) + (d**2) - (b**2) - (c**2)

    return rotation


def _scale_coordinates_by_pixdim(
    coordinates_vol: np.ndarray,
    header: nib.Nifti1Header,
) -> np.ndarray:
    """Convert coordinates to world space by scaling with pixel dimension.

    Parameters
    ----------
    coordinates_vol: ndarray
        Coordinates in the volume (raw data).

    header: Nifti1Header
        Contains image metadata

    Returns
    -------
    ndarray
        Coordinates in the world space

    Notes
    -----
    This method is for compatibility with analyze and is not supposed to be used as the main orientation method.
    But it is used if sform_code = 0. The world coordinates are determined simply by scaling by the voxel size
    by their dimension stored in pixdim.

    References
    ----------
    https://brainder.org/2012/09/23/the-nifti-file-format/
    """
    return np.array(coordinates_vol) * np.array(
        [header["pixdim"][1], header["pixdim"][2], header["pixdim"][3]]
    )
