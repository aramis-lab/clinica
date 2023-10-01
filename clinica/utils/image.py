from os import PathLike
from pathlib import Path
from typing import Callable, Optional, Tuple, Union

import nibabel as nib
import numpy as np
from nibabel.nifti1 import Nifti1Image


def compute_aggregated_volume(
    image_filename: PathLike,
    aggregator: Optional[Callable] = None,
    volumes_to_keep: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Computes the aggregated 3D volumes from a 4D image and an aggregator function.

    The aggregation is computed on the last (fourth) dimension.
    It is possible to compute the aggregation on a subset of the volumes through the
    parameter 'volumes_to_keep'.

    Parameters
    ----------
    image_filename : str
        The path to the input image.

    aggregator : Callable, optional
        The aggregator function. Example: np.average, np.median.
        If None, no aggregation will be performed.

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
    data = [
        volume.get_fdata(dtype="float32").astype(np.float32)
        for volume in volumes_to_keep
    ]
    if aggregator is None:
        return np.stack(data, axis=-1)
    return aggregator(data, axis=0)


def get_new_image_like(old_image: PathLike, new_image_data: np.ndarray) -> Nifti1Image:
    """Builds a new Nifti1Image from the provided image and new data.

    Parameters
    ----------
    old_image : PathLike
        The path to the old image file from which to get the header and affine.

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


def merge_nifti_images_in_time_dimension(
    images: Tuple[Union[str, PathLike], ...], out_file: Optional[PathLike] = None
) -> PathLike:
    """Concatenates the provided images in the 4th dimension.

    The provided images must all be 3D or 4D. For 3D images, a dummy
    4th dimension will be added before concatenation.

    Parameters
    ----------
    images : tuple of str or tuple of Pathlike
        The paths to the images that should get merged.

    out_file : PathLike, optional
        The path to the file in which to write the merged
        volumes. If None, the volumes will be written to a
        file named 'merged_files.nii.gz' in the current folder.

    Returns
    -------
    out_file : PathLike
        Path to merged volumes.
    """
    import os

    images = _check_existence(images)
    out_file = out_file or os.path.abspath("merged_files.nii.gz")
    volumes = _check_volumes_from_images(images)
    merged_volume = np.concatenate(volumes, axis=-1)
    merged_image = get_new_image_like(images[0], merged_volume)
    nib.save(merged_image, out_file)

    return out_file


def _check_existence(filenames: Tuple[PathLike, ...]) -> Tuple[Path, ...]:
    """Converts all element in provided tuple to Path objects and check for existence."""
    from pathlib import Path

    if len(filenames) < 2:
        raise ValueError("At least 2 files are required.")
    filenames = tuple(Path(f) for f in filenames)
    missing_files = tuple(f for f in filenames if not f.exists())
    if missing_files:
        raise FileNotFoundError(f"the following file(s) are missing: {missing_files}")
    return filenames


def _check_volumes_from_images(images: Tuple[Path, ...]) -> Tuple[np.ndarray, ...]:
    """Loads the images and check the dimensions."""
    images = tuple(nib.load(i) for i in images)
    volumes = tuple(i.get_fdata() for i in images)
    four_dimensional_volumes = []
    for volume in volumes:
        if volume.ndim == 3:
            four_dimensional_volumes.append(volume[..., np.newaxis])
        elif volume.ndim == 4:
            four_dimensional_volumes.append(volume)
        else:
            raise ValueError(
                f"Only 3D or 4D images can be concatenated. A {volume.ndim}D image was found."
            )

    return tuple(four_dimensional_volumes)


def merge_nifti_images_in_time_dimension_task(image1: str, image2: str) -> str:
    """Merges the two provided volumes in the time (4th) dimension."""
    # This is needed because Nipype needs to have self-contained functions
    from clinica.utils.image import merge_nifti_images_in_time_dimension  # noqa

    return str(merge_nifti_images_in_time_dimension((image1, image2)))


def remove_dummy_dimension_from_image(image: str, output: str) -> str:
    """Remove all dummy dimensions (i.e. equal to 1) from an image.

    Parameters
    ----------
    image : str
        Path to the input image.

    output : str
        Path to the desired output image.

    Returns
    -------
    str :
        The path to the output image.
    """
    import nibabel as nib
    from nilearn.image import new_img_like

    img = new_img_like(image, nib.load(image).get_fdata().squeeze())
    nib.save(img, output)

    return output
