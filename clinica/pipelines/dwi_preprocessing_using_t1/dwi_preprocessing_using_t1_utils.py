import os
import shutil
from os import PathLike
from typing import Tuple

from clinica.utils.dwi import DWIDataset


def rename_into_caps(
    in_bids_dwi: str,
    fname_dwi: str,
    fname_bval: str,
    fname_bvec: str,
    fname_brainmask: str,
):
    """Rename the outputs of the pipelines into CAPS.

    Parameters
    ----------
    in_bids_dwi : str
        Path to input BIDS DWI to extract the <source_file>

    fname_dwi : str
        Name of preprocessed DWI file.

    fname_bval : str
        Name of preprocessed bval file.

    fname_bvec : str
        Name of preprocessed bvec file.

    fname_brainmask : str
        Name of B0 mask file.

    Returns
    -------
    Tuple[str, str, str, str]
        The different outputs in CAPS format.
    """
    from clinica.utils.dwi import rename_files

    return rename_files(
        in_bids_dwi,
        {
            fname_dwi: "_space-T1w_preproc.nii.gz",
            fname_bval: "_space-T1w_preproc.bval",
            fname_bvec: "_space-T1w_preproc.bval",
            fname_brainmask: "_space-T1w_brainmask.nii.gz",
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


def rotate_b_vectors(b_vectors_filename: str, matrix_filenames: list) -> str:
    """Rotate the B-vectors contained in the input b_vectors_filename file
    according to the provided list of matrices.

    Parameters
    ----------
    b_vectors_filename : str
        Path to the B-vectors file to rotate.

    matrix_filenames : list of str
        List of paths to rotation matrices to apply to the B-vectors.

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

    import numpy as np

    b_vectors_filename = Path(b_vectors_filename)
    stem = (
        Path(b_vectors_filename.stem).stem
        if b_vectors_filename.suffix == ".gz"
        else b_vectors_filename.stem
    )
    rotated_b_vectors_filename = b_vectors_filename.with_name(f"{stem}_rotated.bvec")
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


def ants_apply_transforms(
    fixed_image, moving_image, transforms, warped_image, output_warped_image=True
) -> None:
    import os
    import subprocess

    # Convert to absolute path.
    warped_image = os.path.abspath(warped_image)

    # Whether we want the warped image or transformation field as output.
    output = f"{warped_image}" if output_warped_image else f"[{warped_image}, 1]"

    # Common template for the command.
    cmd = (
        f"antsApplyTransforms -o {output} -i {moving_image} -r {fixed_image} "
        f"-t {transforms[0]} -t {transforms[1]} -t {transforms[2]}"
    )

    # Perform the actual call.
    subprocess.run(cmd, shell=True)

    return warped_image


def init_input_node(t1w, dwi, bvec, bval, dwi_json):
    """Initialize the pipeline."""
    from clinica.utils.dwi import bids_dir_to_fsl_dir, check_dwi_volume
    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        get_subject_id,
        handle_missing_keys_dwi,
    )
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(t1w)
    check_dwi_volume(DWIDataset(dwi=dwi, b_values=bval, b_vectors=bvec))

    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        dwi_json,
        [
            "TotalReadoutTime",
            "PhaseEncodingDirection",
        ],
        handle_missing_keys=handle_missing_keys_dwi,
    )
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    # Print begin message
    print_begin_image(
        image_id,
        ["TotalReadoutTime", "PhaseEncodingDirection"],
        [str(total_readout_time), phase_encoding_direction],
    )

    return (
        image_id,
        t1w,
        dwi,
        bvec,
        bval,
        total_readout_time,
        phase_encoding_direction,
    )


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.ux import print_end_image

    print_end_image(image_id)


def prepare_reference_b0(
    in_dwi: str,
    in_bval: str,
    in_bvec: str,
    low_bval: int = 5,
    working_directory=None,
):
    """Prepare reference b=0 image.

    This function prepares the data for further corrections.
    It co-registers the B0 images and then average them in order to
    obtain only one average B0 image.

    Parameters
    ----------
    in_dwi : str
        Path to the input DWI file.

    in_bvec : str
        Path to the vector file of the diffusion directions of the DWI dataset.

    in_bval : str
        Path to the B-values file.

    low_bval : int, optional
        Threshold for B0 images: images with b<=low_bval will be considered as b0 images.
        Default=5.

    working_directory : str, optional
        Temporary folder results where the results are stored.
        Defaults to None.

    Returns
    -------
    out_reference_b0 : str
        Path to the average of the B0 images or the only B0 image.

    out_b0_dwi_merge : str
        Path to the average of B0 images merged to the DWIs.

    out_updated_bval : str
        Path to the updated gradient values table.

    out_updated_bvec : str
        Path to the updated gradient vectors table.
    """
    from clinica.utils.dwi import b0_dwi_split, check_dwi_dataset, insert_b0_into_dwi

    dwi_dataset = check_dwi_dataset(
        DWIDataset(dwi=in_dwi, b_values=in_bval, b_vectors=in_bvec)
    )
    working_directory = configure_working_directory(dwi_dataset.dwi, working_directory)
    small_b_dataset, large_b_dataset = b0_dwi_split(dwi_dataset, low_bval=low_bval)
    reference_b0 = compute_reference_b0(
        small_b_dataset.dwi, dwi_dataset.b_values, low_bval, working_directory
    )
    reference_dataset = insert_b0_into_dwi(reference_b0, large_b_dataset)

    return str(reference_b0), *tuple(str(_) for _ in reference_dataset)


def compute_reference_b0(
    extracted_b0: PathLike,
    in_bval: PathLike,
    low_bval: int,
    working_directory,
    clean_working_dir: bool = True,
) -> PathLike:
    """Compute the reference B0.

    This function calls the b0_flirt FSL pipeline under the hood.
    The pipeline will write files to the provided working directory.

    .. warning::
        If the option clean_working_dir is set to True, this directory
        will be cleaned after execution of the pipeline.

    Raises
    ------
    ValueError :
        If the number of B0 volumes is <= 0.
    """
    from clinica.utils.dwi import compute_average_b0, count_b0s

    nb_b0s = count_b0s(in_bval=in_bval, low_bval=low_bval)
    if nb_b0s <= 0:
        raise ValueError(
            f"The number of b0s should be strictly positive (b-val file: {in_bval})."
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
    in_dwi: PathLike, working_directory: PathLike = None
) -> PathLike:
    """Configures a temporary working directory for writing the output files of
    the b0 co-registration.

    Parameters
    ----------
    in_dwi : str
        Path to DWI file. This is used to create the name of the folder.

    working_directory : str, optional
        The folder to be used if provided by the user.
        If None, then create a temporary folder to work in.

    Returns
    -------
    str:
        The configured working directory.
    """
    import hashlib
    import tempfile
    from pathlib import Path

    working_directory = working_directory or tempfile.mkdtemp()
    working_directory = (
        Path(working_directory) / hashlib.md5(str(in_dwi).encode()).hexdigest()
    )
    working_directory.mkdir(parents=True)

    return working_directory


def register_b0(
    nb_b0s: int, extracted_b0: PathLike, working_directory: PathLike
) -> PathLike:
    """Run the FSL pipeline 'b0_flirt_pipeline' in order to co-register the b0 images.

    This function is a simple wrapper around the b0_flirt_pipeline which configures it,
    runs it, and returns the path to the merged file of interest.

    Parameters
    ----------
    nb_b0s : int
        The number of B0 volumes in the dataset.

    extracted_b0 : str
        The extracted B0 volumes to co-register.

    working_directory : str
        The working directory in which the pipeline will write intermediary files.

    Returns
    -------
    str:
        The path to the nifti image containing the co-registered B0 volumes.
        If the pipeline ran successfully, this file should be located in :
        working_directory / b0_coregistration / concat_ref_moving / merged_files.nii.gz
    """
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        b0_flirt_pipeline,
    )

    b0_flirt = b0_flirt_pipeline(num_b0s=nb_b0s)
    b0_flirt.inputs.inputnode.in_file = str(extracted_b0)
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


def delete_temp_dirs(checkpoint: str, dir_to_del: list, base_dir: str) -> None:
    """This function deletes the directories of the given list".

    Parameters
    ----------
    checkpoint: str
    Path to a file. Used to ensure, that the temporary directories we want to delete are not useful anymore, and to verify that the subject and session are right.

    dir_to_del: list
    Names of the directories we want to delete.

    base_dir: str
    Path to the working directory.
    """
    import shutil
    from pathlib import Path

    from clinica.utils.stream import cprint

    subject_session_folder_name = extract_sub_ses_folder_name(checkpoint)
    for a in dir_to_del:
        for z in Path(base_dir).rglob(f"*{a}*"):
            if (Path(z).parent).name == subject_session_folder_name:
                shutil.rmtree(z)
                cprint(msg=f"Temporary folder {z} deleted", lvl="info")
