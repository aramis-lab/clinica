"""This module contains utilities used by the DWIPreprocessingUsingT1 pipeline."""

import functools
import shutil
from os import PathLike
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from clinica.pipelines.dwi.utils import DWIDataset

__all__ = [
    "rename_into_caps",
    "change_itk_transform_type",
    "broadcast_matrix_filename_to_match_b_vector_length",
    "ants_warp_image_multi_transform",
    "rotate_b_vectors",
    "init_input_node",
    "print_end_pipeline",
    "prepare_reference_b0",
    "insert_b0_into_dwi",
    "compute_reference_b0",
    "configure_working_directory",
    "register_b0",
    "extract_sub_ses_folder_name",
    "split_dwi_dataset_with_b_values",
]


def rename_into_caps(
    dwi_filename: Path,
    dwi_preproc_filename: Path,
    b_values_preproc_filename: Path,
    b_vectors_preproc_filename: Path,
    b0_brain_mask_filename: Path,
) -> Tuple[str, ...]:
    """Rename the outputs of the pipelines into CAPS.

    Parameters
    ----------
    dwi_filename : str
        The path to the input BIDS DWI to extract the <source_file>.

    dwi_preproc_filename : str
        The path to the preprocessed DWI file.

    dwi_preproc_filename : str
        The path to the preprocessed DWI file.

    b_values_preproc_filename : str
        The path to the preprocessed b-values file.

    b_vectors_preproc_filename : str
        The path to the preprocessed b-vectors file.

    b0_brain_mask_filename : str
        The path to the B0 mask file.

    Returns
    -------
    Tuple[str, str, str, str]
        The different outputs in CAPS format.
    """
    from clinica.pipelines.dwi.utils import rename_files

    return rename_files(
        dwi_filename,
        {
            dwi_preproc_filename: "_space-T1w_desc-preproc_dwi.nii.gz",
            b_values_preproc_filename: "_space-T1w_desc-preproc_dwi.bval",
            b_vectors_preproc_filename: "_space-T1w_desc-preproc_dwi.bvec",
            b0_brain_mask_filename: "_space-T1w_brainmask.nii.gz",
        },
    )


def change_itk_transform_type(input_affine_file: Path) -> Path:
    """Change ITK transform type.

    This function takes in the affine.txt produced by the c3d_affine_tool (which converted
    an FSL FLIRT affine.mat into the affine.txt) it then modifies the 'Transform Type' of
    this affine.txt so that it is compatible with the antsApplyTransforms tool and
    produces a new affine file titled 'updated_affine.txt'.

    Parameters
    ----------
    input_affine_file : Path
        The path to the input affine that should be changed.

    Returns
    -------
    update_affine_file : Path
        The path to the updated affine.
    """
    original_content = input_affine_file.read_text()
    updated_affine_file = input_affine_file.with_name("updated_affine.txt")
    updated_affine_file.write_text(
        original_content.replace(
            "Transform: MatrixOffsetTransformBase_double_3_3",
            "Transform: AffineTransform_double_3_3",
        )
    )

    return updated_affine_file


def broadcast_matrix_filename_to_match_b_vector_length(
    matrix_filename: Path, b_vectors_filename: Path
) -> List[Path]:
    """Return a list of the matrix filename repeated as many times as there are B-vectors."""
    b_vectors = np.loadtxt(b_vectors_filename).T

    return _broadcast_filename_into_list(matrix_filename, len(b_vectors))


def _broadcast_filename_into_list(filename: Path, desired_length: int) -> List[Path]:
    """Return a list of provided filename repeated to match the desired length."""
    return [filename] * desired_length


def ants_warp_image_multi_transform(fix_image, moving_image, ants_warp_affine) -> str:
    import os
    from pathlib import Path

    out_warp = Path("warped_epi.nii.gz").resolve()

    cmd = (
        f"WarpImageMultiTransform 3 {moving_image} {out_warp} -R {fix_image} "
        f"{ants_warp_affine[0]} {ants_warp_affine[1]} {ants_warp_affine[2]}"
    )
    os.system(cmd)

    return str(out_warp)


def rotate_b_vectors(
    b_vectors_filename: Path,
    matrix_filenames: List[str],
    output_dir: Optional[Path] = None,
) -> Path:
    """Rotate the B-vectors contained in the input b_vectors_filename file
    according to the provided list of matrices.

    Parameters
    ----------
    b_vectors_filename : str
        Path to the B-vectors file to rotate.

    matrix_filenames : list of str
        List of paths to rotation matrices to apply to the B-vectors.

    output_dir : str, optional
        If specified, the rotated B-vectors will be written in this folder
        rather than in the same folder as the provided B-vectors.

    Returns
    -------
    rotated_b_vectors_filename : str
        Path to the rotated B-vectors.

    Notes
    -----
    The input affine matrix transforms points in the destination image to their corresponding
    coordinates in the original image. Therefore, this matrix should be inverted first, as
    we want to know the target position of :math:`\\vec{r}`.
    """
    from pathlib import Path

    stem = (
        Path(b_vectors_filename.stem).stem
        if b_vectors_filename.suffix == ".gz"
        else b_vectors_filename.stem
    )
    rotated_b_vectors_filename = b_vectors_filename.with_name(f"{stem}_rotated.bvec")
    if output_dir:
        rotated_b_vectors_filename = output_dir / rotated_b_vectors_filename.name
    else:
        rotated_b_vectors_filename = Path(rotated_b_vectors_filename.name).resolve()

    b_vectors = np.loadtxt(b_vectors_filename).T

    if len(b_vectors) != len(matrix_filenames):
        raise RuntimeError(
            f"Number of b-vectors ({len(b_vectors)}) and rotation "
            f"matrices ({len(matrix_filenames)}) should match."
        )
    rotated_b_vectors = []
    for b_vector, matrix_filename in zip(b_vectors, matrix_filenames):
        if np.all(b_vector == 0.0):
            rotated_b_vectors.append(b_vector)
        else:
            inv_rot = np.linalg.inv(np.loadtxt(matrix_filename))[:3, :3]
            new_b_vector = inv_rot.dot(b_vector)
            rotated_b_vectors.append((new_b_vector / np.linalg.norm(new_b_vector)))

    np.savetxt(rotated_b_vectors_filename, np.array(rotated_b_vectors).T, fmt="%0.15f")

    return rotated_b_vectors_filename


def init_input_node(
    t1w_filename: str,
    dwi_filename: str,
    dwi_json_filename: str,
    b_vectors_filename: str,
    b_values_filename: str,
) -> tuple:
    """Initialize the pipeline.

    Parameters
    ----------
    t1w_filename : str
        The path to the T1w image in BIDS format.

    dwi_filename : str
        The path to the diffusion weighted image in BIDS format.

    dwi_json_filename : str
        The path to the DWI JSON file in BIDS format and containing
        'TotalReadoutTime' and 'PhaseEncodingDirection' metadata
        (see BIDS specifications).

    b_vectors_filename : str
        The path to the b-vectors file in BIDS format.

    b_values_filename : str
        The path of the b-values file in BIDS format.

    Returns
    -------
    image_id : str
        The subject ID extracted from the t1w image path.

    t1w_filename : str
        The path to the T1w image in BIDS format.

    dwi_filename : str
        The path to the diffusion weighted image in BIDS format.

    b_vectors_filename : str
        The path to the b-vectors file in BIDS format.

    b_values_filename : str
        The path of the b-values file in BIDS format.

    total_readout_time : str
        The total readout time extracted from the dwi JSON file.

    phase_encoding_direction : str
        The phase encoding direction for the dwi image, extracted
        from the dwi JSON file.
    """
    from clinica.pipelines.dwi.utils import DWIDataset
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    from ..utils import check_dwi_volume, get_readout_time_and_phase_encoding_direction

    image_id = get_subject_id(t1w_filename)
    check_dwi_volume(
        DWIDataset(
            dwi=dwi_filename,
            b_values=b_values_filename,
            b_vectors=b_vectors_filename,
        )
    )
    (
        total_readout_time,
        phase_encoding_direction,
    ) = get_readout_time_and_phase_encoding_direction(dwi_json_filename)
    print_begin_image(
        image_id,
        ["TotalReadoutTime", "PhaseEncodingDirection"],
        [str(total_readout_time), phase_encoding_direction],
    )

    return (
        image_id,
        t1w_filename,
        dwi_filename,
        b_vectors_filename,
        b_values_filename,
        total_readout_time,
        phase_encoding_direction,
    )


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.ux import print_end_image

    print_end_image(image_id)


def prepare_reference_b0(
    dwi_dataset: DWIDataset,
    b_value_threshold: float = 5.0,
    working_directory: Optional[Path] = None,
) -> Tuple[Path, DWIDataset]:
    """Prepare reference b=0 image.

    This function prepares the data for further corrections.
    It co-registers the B0 images and then average them in order to
    obtain only one average B0 image.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        DWI dataset for which to prepare the reference B0.

    b_value_threshold : float, optional
        Threshold for B0 volumes. Volumes in the DWI image for which the
        corresponding b_value is <= b_value_threshold will be considered
        as b0 volumes.
        Default=5.0.

    working_directory : Path, optional
        Path to temporary folder where results are stored.
        Defaults to None.

    Returns
    -------
    reference_b0 : Path
        Path to the average of the B0 images or the only B0 image.

    reference_dataset : DWIDataset
        DWI dataset containing the B0 images merged and inserted at index 0.
    """
    from ..utils import check_dwi_dataset

    dwi_dataset = check_dwi_dataset(dwi_dataset)
    working_directory = configure_working_directory(dwi_dataset.dwi, working_directory)
    small_b_dataset, large_b_dataset = split_dwi_dataset_with_b_values(
        dwi_dataset,
        b_value_threshold=b_value_threshold,
        working_directory=working_directory,
    )
    reference_b0 = compute_reference_b0(
        small_b_dataset.dwi, dwi_dataset.b_values, b_value_threshold, working_directory
    )
    reference_dataset = insert_b0_into_dwi(reference_b0, large_b_dataset)

    return reference_b0, reference_dataset


def insert_b0_into_dwi(b0_filename: PathLike, dwi_dataset: DWIDataset) -> DWIDataset:
    """Inserts a b0 volume into the DWI dataset as the first volume.

    Also updates the b-values and b-vectors files accordingly.

    Parameters
    ----------
    b0_filename : PathLike
        Path to image file containing one b=0 volume (could be the average of a b0 dataset).

    dwi_dataset : DWIDataset
        DWI dataset in which to insert the b0 volume.

    Returns
    -------
    DWIDataset :
        The diffusion dataset : b0 volume + dwi volumes.
    """
    from ..utils import check_dwi_dataset, check_file

    dwi_dataset = check_dwi_dataset(dwi_dataset)
    b0_filename = check_file(b0_filename)

    return DWIDataset(
        dwi=_insert_b0_into_dwi_image(dwi_dataset.dwi, b0_filename),
        b_values=_insert_b0_into_b_values(dwi_dataset.b_values),
        b_vectors=_insert_b0_into_b_vectors(dwi_dataset.b_vectors),
    )


def _insert_b0_into_dwi_image(dwi_filename: Path, b0_filename: Path) -> Path:
    """Insert the provided B0 volume into the DWI image.

    This insertion is done at index 0 along the 4th / time dimension.
    """
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    from ..utils import add_suffix_to_filename

    output_dwi_filename = merge_nifti_images_in_time_dimension(
        (b0_filename, dwi_filename),
        out_file=add_suffix_to_filename(
            _remove_entity_from_filename(dwi_filename, "large_b"),
            "merged",
        ),
    )

    return output_dwi_filename


def _remove_entity_from_filename(filename: Path, entity: str) -> Path:
    """Removes the provided entity from the given file name.

    Parameters
    ----------
    filename : Path
        The path object on which to operate removal of entity.

    entity : str
        The entity which should be removed.

    Returns
    -------
    Path :
        The new Path object without the entity.
    """
    return Path(filename.parent / filename.name.replace(f"_{entity}", ""))


def _insert_b0_into_b_values(b_values_filename: Path) -> Path:
    """Insert a 0 value at index 0 into the b-values file."""
    from ..utils import add_suffix_to_filename

    b_values = np.loadtxt(b_values_filename)
    b_values = np.insert(b_values, 0, 0)
    output_b_values_filename = add_suffix_to_filename(
        _remove_entity_from_filename(b_values_filename, "large_b"), "merged"
    )
    _write_b_values(output_b_values_filename, b_values)

    return output_b_values_filename


def _write_numpy(filename: Path, data: np.ndarray, fmt: str, delimiter: str) -> None:
    """Writes the provided array to the provided filename using the provided formatting."""
    np.savetxt(filename, data, fmt=fmt, delimiter=delimiter)


_write_b_vectors = functools.partial(_write_numpy, fmt="%10.5f", delimiter=" ")
_write_b_values = functools.partial(_write_numpy, fmt="%d", delimiter=" ")


def _insert_b0_into_b_vectors(b_vectors_filename: Path) -> Path:
    """Insert a 0 vector into the b-vectors file."""
    from ..utils import add_suffix_to_filename

    b_vectors = np.loadtxt(b_vectors_filename)
    b_vectors = np.insert(b_vectors, 0, 0.0, axis=1)
    output_b_vectors_filename = add_suffix_to_filename(
        _remove_entity_from_filename(b_vectors_filename, "large_b"), "merged"
    )
    _write_b_vectors(output_b_vectors_filename, b_vectors)

    return output_b_vectors_filename


def compute_reference_b0(
    extracted_b0: Path,
    b_values_filename: Path,
    b_value_threshold: float,
    working_directory: Path,
    clean_working_dir: bool = False,
) -> Path:
    """Compute the reference B0.

    This function calls the b0_flirt FSL pipeline under the hood.
    The pipeline will write files to the provided working directory.

    .. warning::
        If the option clean_working_dir is set to True, this directory
        will be cleaned after execution of the pipeline.

    Parameters
    ----------
    extracted_b0 : Path
        The path to the image of the extracted B0 volume.

    b_values_filename : Path
        The path to the b-values file.

    b_value_threshold : float
        Threshold on b-values to decide whether a DWI volume is
        B0 or not.

    working_directory : Path
        Path to the directory where output files should be written.

    clean_working_dir : bool, optional
        If True, the working directory will be cleaned at the end of
        the execution.
        Default=False.

    Returns
    -------
    registered_b0_filename : Path
        Path to the output image holding the registered B0.

    Raises
    ------
    ValueError :
        If the number of B0 volumes is <= 0.
    """
    from ..utils import compute_average_b0

    if (
        nb_b0s := _count_b0s(
            b_value_filename=b_values_filename, b_value_threshold=b_value_threshold
        )
    ) <= 0:
        raise ValueError(
            f"The number of b0s should be strictly positive (b-val file: {b_values_filename})."
        )
    if nb_b0s == 1:
        return extracted_b0
    registered_b0s = register_b0(
        nb_b0s,
        extracted_b0,
        working_directory=working_directory,
    )
    out_reference_b0 = compute_average_b0(registered_b0s, squeeze=False)
    registered_b0_file_name = extracted_b0.with_name("reference_b0_volume.nii.gz")
    shutil.copy(out_reference_b0, registered_b0_file_name)
    if clean_working_dir:
        shutil.rmtree(working_directory)

    return registered_b0_file_name


def configure_working_directory(
    dwi_filename: Path,
    working_directory: Optional[Path] = None,
) -> Path:
    """Configures a temporary working directory for writing the output files of
    the b0 co-registration.

    Parameters
    ----------
    dwi_filename : Path
        The path to the DWI image file. This is used to create the name of the folder.

    working_directory : Path, optional
        The folder to be used if provided by the user.
        If None, then create a temporary folder to work in.

    Returns
    -------
    Path:
        The path to the configured working directory.
    """
    import hashlib
    import tempfile
    from pathlib import Path

    working_directory = working_directory or Path(tempfile.mkdtemp())
    working_directory = (
        working_directory / hashlib.md5(str(dwi_filename).encode()).hexdigest()
    )
    working_directory.mkdir(parents=True, exist_ok=True)

    return working_directory


def register_b0(
    nb_b0s: int,
    extracted_b0_filename: Path,
    working_directory: Path,
) -> Path:
    """Run the FSL pipeline 'b0_flirt_pipeline' in order to co-register the b0 images.

    This function is a simple wrapper around the b0_flirt_pipeline which configures it,
    runs it, and returns the path to the merged file of interest.

    Parameters
    ----------
    nb_b0s : int
        The number of B0 volumes in the dataset.

    extracted_b0_filename : Path
        The extracted B0 volumes to co-register.

    working_directory : Path
        The working directory in which the pipeline will write intermediary files.

    Returns
    -------
    Path:
        The path to the nifti image containing the co-registered B0 volumes.
        If the pipeline ran successfully, this file should be located in :
        working_directory / b0_coregistration / concat_ref_moving / merged_files.nii.gz
    """
    from .workflows import b0_flirt_pipeline

    b0_flirt = b0_flirt_pipeline(num_b0s=nb_b0s)
    b0_flirt.inputs.inputnode.in_file = str(extracted_b0_filename)
    b0_flirt.base_dir = str(working_directory)
    b0_flirt.run()

    return (
        working_directory
        / "b0_coregistration"
        / "concat_ref_moving"
        / "merged_files.nii.gz"
    )


def extract_sub_ses_folder_name(file_path: str) -> str:
    """This function extracts the name of the folder corresponding to a subject and a session.

    Parameters
    ----------
    file_path: str
        Path to a temporary file for the subject and session of interest.

    Returns
    -------
    str:
    Name of the folder corresponding to a subject and a session.

    Examples
    --------
    >>> extract_sub_ses_folder_name("/localdrive10TB/users/matthieu.joulot/wd/dwi-preprocessing-using-t1/epi_pipeline/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/MergeDWIs/Jacobian_image_maths_thresh_merged.nii.gz")
    4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea
    """
    from pathlib import Path

    return (Path(Path(file_path).parent).parent).name


def split_dwi_dataset_with_b_values(
    dwi_dataset: DWIDataset,
    b_value_threshold: float = 5.0,
    working_directory: Optional[Path] = None,
) -> Tuple[DWIDataset, DWIDataset]:
    """Splits the DWI dataset in two through the B-values.

    Splits the DWI volumes into two datasets :
     - the first dataset is relative to volumes having a b-value <= b_value_threshold.
     - the second dataset is relative to volumes having a b-value > b_value_threshold.

    The function writes 6 files (3 files for each dataset), and returns a
    length 2 tuple containing the two DWI datasets.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        DWI image dataset to split.

    b_value_threshold : float, optional
        Defines the b0 volumes as all volumes b-value <= b_value_threshold.
        Defaults to 5.0.

    working_directory : Path, optional
        An optional working directory in which to write the small and large
        DWI datasets. If not provided, they will be written in the same
        directory as the provided DWI dataset.

    Returns
    -------
    small_b_dataset : DWIDataset
        DWI dataset for volumes for which b-value <= b_value_threshold.

    large_b_dataset : DWIDataset
        DWI dataset for volumes for which b-value > b_value_threshold.

    Raises
    ------
    ValueError:
        If b_value_threshold < 0.
    """
    from ..utils import check_b_value_threshold, check_dwi_dataset, get_b0_filter

    check_b_value_threshold(b_value_threshold)
    dwi_dataset = check_dwi_dataset(dwi_dataset)
    b_values, _ = _check_b_values_and_b_vectors(dwi_dataset)
    small_b_filter = get_b0_filter(
        dwi_dataset.b_values, b_value_threshold=b_value_threshold
    )
    large_b_filter = np.array(
        [i for i in range(len(b_values)) if i not in small_b_filter]
    )

    return (
        _build_dwi_dataset_from_filter(
            dwi_dataset,
            "small_b",
            small_b_filter,
            working_directory=working_directory,
        ),
        _build_dwi_dataset_from_filter(
            dwi_dataset,
            "large_b",
            large_b_filter,
            working_directory=working_directory,
        ),
    )


def _check_b_values_and_b_vectors(
    dwi_dataset: DWIDataset,
) -> Tuple[np.ndarray, np.ndarray]:
    """Opens the b-values and b-vectors files and transpose the b-vectors if needed."""
    import warnings

    b_values = np.loadtxt(dwi_dataset.b_values)
    b_vectors = np.loadtxt(dwi_dataset.b_vectors)
    if b_values.shape[0] == b_vectors.shape[0]:
        warnings.warn(
            "Warning: The b-vectors file should be column-wise. The b-vectors will be transposed",
            UserWarning,
        )
        b_vectors = b_vectors.T

    return b_values, b_vectors


def _build_dwi_dataset_from_filter(
    dwi_dataset: DWIDataset,
    filter_name: str,
    filter_array: np.ndarray,
    working_directory: Optional[Path] = None,
) -> DWIDataset:
    """Builds a new DWI dataset from a given DWI dataset and a filter.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        The DWI dataset to filter.

    filter_name : str
        The name of the filter. This will be used to build the
        file names associated with the new dataset.

    filter_array : np.ndarray
        1D array of indices to filter the DWI dataset.

    working_directory : Path, optional
        An optional working directory in which to write the filtered
        DWI dataset. If not provided, it will be written in the same
        directory as the provided DWI dataset.

    Returns
    -------
    DWIDataset :
        The new filtered DWI dataset.
    """
    from ..utils import check_dwi_dataset

    return check_dwi_dataset(
        DWIDataset(
            dwi=_filter_dwi(
                dwi_dataset,
                filter_name,
                filter_array,
                working_directory=working_directory,
            ),
            b_values=_filter_b_values(
                dwi_dataset,
                filter_name,
                filter_array,
                working_directory=working_directory,
            ),
            b_vectors=_filter_b_vectors(
                dwi_dataset,
                filter_name,
                filter_array,
                working_directory=working_directory,
            ),
        )
    )


def _filter_dwi(
    dwi_dataset: DWIDataset,
    filter_name: str,
    filter_array: np.ndarray,
    working_directory: Optional[Path] = None,
) -> Path:
    """Filters the dwi component of the provided DWI dataset."""
    from clinica.utils.image import compute_aggregated_volume, get_new_image_like

    from ..utils import add_suffix_to_filename

    dwi_filename = add_suffix_to_filename(dwi_dataset.dwi, filter_name)
    if working_directory:
        dwi_filename = working_directory / dwi_filename.name
    data = compute_aggregated_volume(
        dwi_dataset.dwi, aggregator=None, volumes_to_keep=filter_array
    )
    img = get_new_image_like(dwi_dataset.dwi, data)
    img.to_filename(dwi_filename)

    return dwi_filename


def _filter_b_values(
    dwi_dataset: DWIDataset,
    filter_name: str,
    filter_array: np.ndarray,
    working_directory: Optional[Path] = None,
) -> Path:
    """Filters the b-values component of the provided DWI dataset."""
    from ..utils import add_suffix_to_filename

    b_values, _ = _check_b_values_and_b_vectors(dwi_dataset)
    b_values_filename = add_suffix_to_filename(dwi_dataset.b_values, filter_name)
    if working_directory:
        b_values_filename = working_directory / b_values_filename.name
    _write_b_values(b_values_filename, b_values[filter_array])

    return b_values_filename


def _filter_b_vectors(
    dwi_dataset: DWIDataset,
    filter_name: str,
    filter_array: np.ndarray,
    working_directory: Optional[Path] = None,
) -> Path:
    """Filters the b-vectors component of the provided DWI dataset."""
    from ..utils import add_suffix_to_filename

    _, b_vectors = _check_b_values_and_b_vectors(dwi_dataset)
    b_vectors_filename = add_suffix_to_filename(dwi_dataset.b_vectors, filter_name)
    if working_directory:
        b_vectors_filename = working_directory / b_vectors_filename.name
    _write_b_vectors(b_vectors_filename, np.array([b[filter_array] for b in b_vectors]))

    return b_vectors_filename


def _count_b0s(b_value_filename: PathLike, b_value_threshold: float = 5.0) -> int:
    """Counts the number of volumes where b<=low_bval.

    Parameters
    ----------
    b_value_filename : PathLike
        Path to the b-value file.

    b_value_threshold : float, optional
        Defines the threshold for the b0 volumes as all volumes which
        satisfy b-value <= b_value_threshold.
        Defaults to 5.0.

    Returns
    -------
    num_b0s : int
        Number of b0 volumes.
    """
    from ..utils import get_b0_filter

    return len(get_b0_filter(b_value_filename, b_value_threshold=b_value_threshold))
