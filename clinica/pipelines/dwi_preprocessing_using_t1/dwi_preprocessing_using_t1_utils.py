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

    # Extract image ID
    image_id = get_subject_id(t1w)

    # Check that the number of DWI, bvec & bval are the same
    check_dwi_volume(dwi, bvec, bval)

    # Read metadata from DWI JSON file:
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


def prepare_reference_b0(in_dwi, in_bval, in_bvec, low_bval=5, working_directory=None):
    """Prepare reference b=0 image.

    This function prepare the data for further corrections. It co-registers the B0 images
    and then average it in order to obtain only one average B0 images.

    Args:
        in_dwi (str): Input DWI file.
        in_bvec (str): Vector file of the diffusion directions of the DWI dataset.
        in_bval (str): B-values file.
        low_bval (optional, int): Set b<=low_bval such that images are considered b0. Defaults to 5.
        working_directory (str): Temporary folder results where the results are stored. Defaults to None.

    Returns:
        out_reference_b0 (str): Average of the B0 images or the only B0 image.
        out_b0_dwi_merge (str): Average of B0 images merged to the DWIs.
        out_updated_bval (str): Updated gradient values table.
        out_updated_bvec (str): Updated gradient vectors table.
    """
    import hashlib
    import os
    import tempfile

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        b0_flirt_pipeline,
    )
    from clinica.utils.dwi import (
        b0_average,
        b0_dwi_split,
        count_b0s,
        insert_b0_into_dwi,
    )

    # Count the number of b0s
    nb_b0s = count_b0s(in_bval=in_bval, low_bval=low_bval)

    # Split dataset into two datasets: the b0 and the b>low_bval datasets
    [extracted_b0, out_split_dwi, out_split_bval, out_split_bvec] = b0_dwi_split(
        in_dwi=in_dwi, in_bval=in_bval, in_bvec=in_bvec, low_bval=low_bval
    )

    if nb_b0s == 1:
        # The reference b0 is the extracted b0
        # cprint('Only one b0 for %s' % in_dwi)
        out_reference_b0 = extracted_b0
    elif nb_b0s > 1:
        # Register the b0 onto the first b0
        b0_flirt = b0_flirt_pipeline(num_b0s=nb_b0s)
        b0_flirt.inputs.inputnode.in_file = extracted_b0
        if working_directory is None:
            working_directory = tempfile.mkdtemp()
        tmp_dir = os.path.join(
            working_directory, hashlib.md5(in_dwi.encode()).hexdigest()
        )
        b0_flirt.base_dir = tmp_dir
        b0_flirt.run()
        # BUG: Nipype does allow to extract the output after running the
        # workflow: we need to 'guess' where the output will be generated
        # out_node = b0_flirt.get_node('outputnode')
        registered_b0s = os.path.abspath(
            os.path.join(
                tmp_dir, "b0_coregistration", "concat_ref_moving", "merged_files.nii.gz"
            )
        )
        # cprint('B0 s will be averaged (file = ' + registered_b0s + ')')
        # Average the b0s to obtain the reference b0
        out_reference_b0 = b0_average(in_file=registered_b0s)
    else:
        raise ValueError(
            f"The number of b0s should be strictly positive (b-val file: {in_bval})."
        )

    # Merge datasets such that bval(DWI) = (0 b1 ... bn)
    [out_b0_dwi_merge, out_updated_bval, out_updated_bvec] = insert_b0_into_dwi(
        in_b0=out_reference_b0,
        in_dwi=out_split_dwi,
        in_bval=out_split_bval,
        in_bvec=out_split_bvec,
    )

    return out_reference_b0, out_b0_dwi_merge, out_updated_bval, out_updated_bvec


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
