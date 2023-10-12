import shutil
from pathlib import Path
from typing import Optional, Tuple

from clinica.utils.dwi import DWIDataset


def rename_into_caps(
    dwi_filename: str,
    dwi_preproc_filename: str,
    b_values_preproc_filename: str,
    b_vectors_preproc_filename: str,
    b0_brain_mask_filename: str,
) -> tuple:
    """Rename the outputs of the pipelines into CAPS.

    Parameters
    ----------
    dwi_filename : str
        The path to the input BIDS DWI to extract the <source_file>.

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
    from clinica.utils.dwi import rename_files

    return rename_files(
        dwi_filename,
        {
            dwi_preproc_filename: "_space-T1w_preproc.nii.gz",
            b_values_preproc_filename: "_space-T1w_preproc.bval",
            b_vectors_preproc_filename: "_space-T1w_preproc.bvec",
            b0_brain_mask_filename: "_space-T1w_brainmask.nii.gz",
        },
    )


def change_itk_transform_type(input_affine_file: str) -> str:
    """Change ITK transform type.

    This function takes in the affine.txt produced by the c3d_affine_tool (which converted
    an FSL FLIRT affine.mat into the affine.txt) it then modifies the 'Transform Type' of
    this affine.txt so that it is compatible with the antsApplyTransforms tool and
    produces a new affine file titled 'updated_affine.txt'.

    Parameters
    ----------
    input_affine_file : str
        Path to the input affine that should be changed.

    Returns
    -------
    update_affine_file : str
        Path to the updated affine.
    """
    from pathlib import Path

    input_affine_file = Path(input_affine_file)
    original_content = input_affine_file.read_text()
    updated_affine_file = input_affine_file.with_name("updated_affine.txt")
    updated_affine_file.write_text(
        original_content.replace(
            "Transform: MatrixOffsetTransformBase_double_3_3",
            "Transform: AffineTransform_double_3_3",
        )
    )

    return str(updated_affine_file)


def broadcast_matrix_filename_to_match_b_vector_length(
    matrix_filename: str, b_vectors_filename: str
) -> list:
    """Return a list of the matrix filename repeated as many times as there are B-vectors."""
    import numpy as np

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import (  # noqa
        broadcast_filename_into_list,
    )

    b_vectors = np.loadtxt(b_vectors_filename).T

    return broadcast_filename_into_list(matrix_filename, len(b_vectors))


def broadcast_filename_into_list(filename: str, desired_length: int) -> list:
    """Return a list of provided filename repeated to match the desired length."""
    return [filename] * desired_length


def ants_warp_image_multi_transform(fix_image, moving_image, ants_warp_affine):
    import os

    out_warp = os.path.abspath("warped_epi.nii.gz")

    cmd = (
        f"WarpImageMultiTransform 3 {moving_image} {out_warp} -R {fix_image} "
        f"{ants_warp_affine[0]} {ants_warp_affine[1]} {ants_warp_affine[2]}"
    )
    os.system(cmd)

    return out_warp


def rotate_b_vectors(
    b_vectors_filename: str, matrix_filenames: list, output_dir: str = None
) -> str:
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
    import os
    from pathlib import Path

    import numpy as np

    b_vectors_filename = Path(b_vectors_filename)
    stem = (
        Path(b_vectors_filename.stem).stem
        if b_vectors_filename.suffix == ".gz"
        else b_vectors_filename.stem
    )
    rotated_b_vectors_filename = b_vectors_filename.with_name(f"{stem}_rotated.bvec")
    if output_dir:
        rotated_b_vectors_filename = Path(output_dir) / rotated_b_vectors_filename.name
    else:
        rotated_b_vectors_filename = os.path.abspath(rotated_b_vectors_filename.name)

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

    return str(rotated_b_vectors_filename)


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
    from clinica.utils.dwi import (
        DWIDataset,
        check_dwi_volume,
        get_readout_time_and_phase_encoding_direction,
    )
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(t1w_filename)
    check_dwi_volume(
        DWIDataset(
            dwi=dwi_filename, b_values=b_values_filename, b_vectors=b_vectors_filename
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


def prepare_reference_b0_task(
    dwi_filename: str,
    b_values_filename: str,
    b_vectors_filename: str,
    b_value_threshold: float = 5.0,
    working_directory=None,
) -> tuple:
    """Task called be Nipype to execute prepare_reference_b0."""
    from pathlib import Path

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import (  # noqa
        prepare_reference_b0,
    )
    from clinica.utils.dwi import DWIDataset

    if working_directory:
        working_directory = Path(working_directory)

    b0_reference_filename, reference_dwi_dataset = prepare_reference_b0(
        DWIDataset(
            dwi=dwi_filename,
            b_values=b_values_filename,
            b_vectors=b_vectors_filename,
        ),
        b_value_threshold,
        working_directory,
    )

    return str(b0_reference_filename), *(str(_) for _ in reference_dwi_dataset)


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
    from clinica.utils.dwi import (
        check_dwi_dataset,
        insert_b0_into_dwi,
        split_dwi_dataset_with_b_values,
    )

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
        The path to the directory where output files should be written.

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
    from clinica.utils.dwi import compute_average_b0, count_b0s

    nb_b0s = count_b0s(
        b_value_filename=b_values_filename, b_value_threshold=b_value_threshold
    )
    if nb_b0s <= 0:
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
    Path :
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
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        b0_flirt_pipeline,
    )

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


def generate_index_file_task(
    b_values_filename: str,
    image_id=None,
    output_dir=None,
) -> str:
    """Wrapper for Nipype."""
    from pathlib import Path

    from clinica.utils.dwi import generate_index_file

    return str(
        generate_index_file(
            Path(b_values_filename),
            image_id,
            Path(output_dir) if output_dir else None,
        )
    )


def generate_acq_file_task(
    dwi_filename: str,
    fsl_phase_encoding_direction: str,
    total_readout_time: str,
    image_id=None,
    output_dir=None,
) -> str:
    """Wrapper for Nipype."""
    from pathlib import Path

    from clinica.utils.dwi import generate_acq_file

    return str(
        generate_acq_file(
            Path(dwi_filename),
            fsl_phase_encoding_direction,
            total_readout_time,
            image_id,
            Path(output_dir) if output_dir else None,
        )
    )
