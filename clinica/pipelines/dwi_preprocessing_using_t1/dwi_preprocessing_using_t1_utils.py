def rename_into_caps(in_bids_dwi, fname_dwi, fname_bval, fname_bvec, fname_brainmask):
    """Rename the outputs of the pipelines into CAPS.


    Args:
        in_bids_dwi (str): Input BIDS DWI to extract the <source_file>
        fname_dwi (str): Preprocessed DWI file.
        fname_bval (str): Preprocessed bval.
        fname_bvec (str): Preprocessed bvec.
        fname_brainmask (str): B0 mask.

    Returns:
        Tuple[str, str, str, str]: The different outputs in CAPS format.
    """
    import os

    from nipype.interfaces.utility import Rename
    from nipype.utils.filemanip import split_filename

    # Extract <source_file> in format sub-CLNC01_ses-M00[_acq-label]_dwi
    _, source_file_dwi, _ = split_filename(in_bids_dwi)

    # Extract base path from fname:
    base_dir_dwi, _, _ = split_filename(fname_dwi)
    base_dir_bval, _, _ = split_filename(fname_bval)
    base_dir_bvec, _, _ = split_filename(fname_bvec)
    base_dir_brainmask, _, _ = split_filename(fname_brainmask)

    # Rename into CAPS DWI:
    rename_dwi = Rename()
    rename_dwi.inputs.in_file = fname_dwi
    rename_dwi.inputs.format_string = os.path.join(
        base_dir_dwi, f"{source_file_dwi}_space-T1w_preproc.nii.gz"
    )
    out_caps_dwi = rename_dwi.run()

    # Rename into CAPS bval:
    rename_bval = Rename()
    rename_bval.inputs.in_file = fname_bval
    rename_bval.inputs.format_string = os.path.join(
        base_dir_bval, f"{source_file_dwi}_space-T1w_preproc.bval"
    )
    out_caps_bval = rename_bval.run()

    # Rename into CAPS bvec:
    rename_bvec = Rename()
    rename_bvec.inputs.in_file = fname_bvec
    rename_bvec.inputs.format_string = os.path.join(
        base_dir_bvec, f"{source_file_dwi}_space-T1w_preproc.bvec"
    )
    out_caps_bvec = rename_bvec.run()

    # Rename into CAPS DWI:
    rename_brainmask = Rename()
    rename_brainmask.inputs.in_file = fname_brainmask
    rename_brainmask.inputs.format_string = os.path.join(
        base_dir_brainmask, f"{source_file_dwi}_space-T1w_brainmask.nii.gz"
    )
    out_caps_brainmask = rename_brainmask.run()

    return (
        out_caps_dwi.outputs.out_file,
        out_caps_bval.outputs.out_file,
        out_caps_bvec.outputs.out_file,
        out_caps_brainmask.outputs.out_file,
    )


def change_itk_transform_type(input_affine_file):
    """Change ITK transform type.

    This function takes in the affine.txt produced by the c3d_affine_tool (which converted
    an FSL FLIRT affine.mat into the affine.txt) it then modifies the 'Transform Type' of
    this affine.txt so that it is compatible with the antsApplyTransforms tool and
    produces a new affine file titled 'updated_affine.txt'.

    """
    import os

    new_file_lines = []

    with open(input_affine_file) as f:
        for line in f:
            if "Transform:" in line:
                if "MatrixOffsetTransformBase_double_3_3" in line:
                    transform_line = "Transform: AffineTransform_double_3_3\n"
                    new_file_lines.append(transform_line)
            else:
                new_file_lines.append(line)

    updated_affine_file = os.path.join(os.getcwd(), "updated_affine.txt")

    with open(updated_affine_file, "wt") as f:
        for line in new_file_lines:
            f.write(line)

    return updated_affine_file


def expend_matrix_list(in_matrix, in_bvec):
    import numpy as np

    bvecs = np.loadtxt(in_bvec).T
    out_matrix_list = [in_matrix]

    out_matrix_list = out_matrix_list * len(bvecs)

    return out_matrix_list


def ants_warp_image_multi_transform(fix_image, moving_image, ants_warp_affine):
    import os

    out_warp = os.path.abspath("warped_epi.nii.gz")

    cmd = (
        f"WarpImageMultiTransform 3 {moving_image} {out_warp} -R {fix_image} "
        f"{ants_warp_affine[0]} {ants_warp_affine[1]} {ants_warp_affine[2]}"
    )
    os.system(cmd)

    return out_warp


def rotate_bvecs(in_bvec, in_matrix):
    """Rotate the input bvec file accordingly with a list of matrices.

    Notes:
        The input affine matrix transforms points in the destination image to their corresponding
        coordinates in the original image. Therefore, this matrix should be inverted first, as
        we want to know the target position of :math:`\\vec{r}`.

    """
    import os

    import numpy as np

    name, fext = os.path.splitext(os.path.basename(in_bvec))
    if fext == ".gz":
        name, _ = os.path.splitext(name)
    out_file = os.path.abspath(f"{name}_rotated.bvec")
    # Warning, bvecs.txt are not in the good configuration, need to put '.T'
    bvecs = np.loadtxt(in_bvec).T
    new_bvecs = []

    if len(bvecs) != len(in_matrix):
        raise RuntimeError(
            f"Number of b-vectors ({len(bvecs)}) and rotation matrices ({len(in_matrix)}) should match."
        )

    for bvec, mat in zip(bvecs, in_matrix):
        if np.all(bvec == 0.0):
            new_bvecs.append(bvec)
        else:
            invrot = np.linalg.inv(np.loadtxt(mat))[:3, :3]
            newbvec = invrot.dot(bvec)
            new_bvecs.append((newbvec / np.linalg.norm(newbvec)))

    np.savetxt(out_file, np.array(new_bvecs).T, fmt="%0.15f")
    return out_file


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
    from clinica.utils.filemanip import extract_metadata_from_json, get_subject_id
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(t1w)

    # Check that the number of DWI, bvec & bval are the same
    check_dwi_volume(dwi, bvec, bval)

    # Read metadata from DWI JSON file:
    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        dwi_json, ["TotalReadoutTime", "PhaseEncodingDirection"]
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
    Path to a file. Used to ensure, that the tempory directories we want to delete are not useful anymore, and to verify that the subject and session are right.

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
