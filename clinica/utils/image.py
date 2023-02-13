import functools
import warnings
from os import PathLike

import nibabel as nib
import numpy as np
from nibabel.nifti1 import Nifti1Image


def compute_aggregated_volume(
    image_filename: PathLike, aggregator, volumes_to_keep=None
) -> np.ndarray:
    """Computes the aggregated 3D volumes from a 4D image and an aggregator function.

    The aggregation is computes on the last (fourth) dimension.
    It is possible to compute the aggregation on a subset of the volumes through the
    parameter 'volumes_to_keep'.

    Parameters
    ----------
    image_filename : str
        The path to the input image.

    aggregator : Callable
        The aggregator function. Example: np.average, np.median.

    volumes_to_keep : np.array, optional
        The volumes to be kept during aggregation. This is a 1D index array.
        If None, all volumes of the input image will be kept during aggregation.
        Default = None.

    Returns
    -------
    np.ndarray:
        The 3D volume data array obtained from the aggregation.
    """
    volumes = np.array(nib.four_to_three(nib.load(image_filename)))
    if volumes_to_keep is None:
        volumes_to_keep = volumes
    else:
        volumes_to_keep = volumes[volumes_to_keep]

    return aggregator(
        [
            volume.get_fdata(dtype="float32").astype(np.float32)
            for volume in volumes_to_keep
        ],
        axis=0,
    )


def get_new_image_like(
    old_image: Nifti1Image, new_image_data: np.ndarray
) -> Nifti1Image:
    """Build a new Nifti1Image from the provided image and new data.

    Parameters
    ----------
    old_image : Nifti1Image
        The old image from which to get the header and affine.

    new_image_data : np.ndarray
        The data for the new image to build.

    Returns
    -------
    Nifti1Image :
        The new image.
    """
    old_img = nib.load(old_image)
    hdr = old_img.header.copy()
    hdr.set_data_shape(new_image_data.shape)
    hdr.set_xyzt_units("mm")
    hdr.set_data_dtype(np.float32)

    return nib.Nifti1Image(new_image_data, old_img.affine, hdr)


def merge_volumes(
    volume1: PathLike, volume2: PathLike, axis, out_file: str = None
) -> str:
    """Merge 'volume1' and 'volume2' in the provided 'axis' dimension.

    Parameters
    ----------
    volume1 : PathLike
        Path to the image containing the first set of volumes.

    volume2 : PathLike
        Path to the image containing the second set of volumes.

    axis : int
        The axis along which the concatenation should be done.

    out_file : str, optional
        The file on which the concatenated image should be written.
        If None, the file name will be 'merge_files.nii.gz', and the
        base path of volume1 will be used. In this case, if volume1
        and volume2 have different paths, the base path of volume1
        will still be used but a warning will be given.

    Returns
    -------
    out_file : str
        Path to the image containing the two sets of volumes merged.
    """
    from pathlib import Path

    if axis not in range(-1, 4):
        raise ValueError(
            f"Axis should be an integer value in [-1, 3]. You provided {axis}."
        )

    volume1 = Path(volume1)
    volume2 = Path(volume2)
    if out_file is None:
        out_file = volume1.parent / "merged_files.nii.gz"
        if volume1.parent != volume2.parent:
            warnings.warn(
                f"Merging volumes {volume1} and {volume2} in {out_file}."
                "Please provide a custom location through the out_file "
                "argument if needed."
            )
    img1 = nib.load(volume1)
    img2 = nib.load(volume2)
    vol = np.concatenate((img1.get_fdata(), img2.get_fdata()), axis=axis)
    new = get_new_image_like(volume1, vol)
    nib.save(new, out_file)

    return out_file


merge_volumes_time_dimension = functools.partial(merge_volumes, axis=-1)
