from pathlib import Path

import pydra
import toimportfunctions
from pydra.mark import annotate, task

working_dir = "/localdrive10TB/users/matthieu.joulot/werk"


@task
@annotate(
    {
        "return": {
            "image_id": str,
            "t1w": str,
            "dwi": str,
            "bvec": str,
            "bval": str,
            "total_readout_time": str,
            "phase_encoding_direction": str,
        }
    }
)
def init_input_node(t1w, dwi, bvec, bval, dwi_json):
    """Initialize the pipeline."""

    # Extract image ID
    image_id = toimportfunctions.get_subject_id(t1w)

    # Check that the number of DWI, bvec & bval are the same
    toimportfunctions.check_dwi_volume(dwi, bvec, bval)

    # Read metadata from DWI JSON file:
    [
        total_readout_time,
        phase_encoding_direction,
    ] = toimportfunctions.extract_metadata_from_json(
        dwi_json, ["TotalReadoutTime", "PhaseEncodingDirection"]
    )
    phase_encoding_direction = toimportfunctions.bids_dir_to_fsl_dir(
        phase_encoding_direction
    )

    # Print begin message
    toimportfunctions.print_begin_image(
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


@task
@annotate(
    {
        "return": {
            "out_reference_b0": str,
            "out_b0_dwi_merge": str,
            "out_updated_bval": str,
            "out_updated_bvec": str,
        }
    }
)
def prepare_reference_b0_1(
    in_dwi, in_bval, in_bvec, low_bval=5, working_directory=working_dir
):
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

    # Count the number of b0s
    nb_b0s = toimportfunctions.count_b0s(in_bval=in_bval, low_bval=low_bval)

    # Split dataset into two datasets: the b0 and the b>low_bval datasets
    [
        extracted_b0,
        out_split_dwi,
        out_split_bval,
        out_split_bvec,
    ] = toimportfunctions.b0_dwi_split(
        in_dwi=in_dwi, in_bval=in_bval, in_bvec=in_bvec, low_bval=low_bval
    )

    if nb_b0s == 1:
        # The reference b0 is the extracted b0
        # cprint('Only one b0 for %s' % in_dwi)
        out_reference_b0 = extracted_b0
    elif nb_b0s > 1:
        if working_directory is None:
            working_directory = tempfile.mkdtemp()
        tmp_dir = os.path.join(
            working_directory, hashlib.md5(in_dwi.encode()).hexdigest()
        )
        b0_flirt_results = toimportfunctions.build_n_run_flirt(
            nb_b0s, extracted_b0, tmp_dir
        )
        out_reference_b0 = toimportfunctions.b0_average(
            in_file=b0_flirt_results[1].get_output_field("out_file")
        )
    else:
        raise ValueError(
            f"The number of b0s should be strictly positive (b-val file: {in_bval})."
        )
    # Merge datasets such that bval(DWI) = (0 b1 ... bn)
    [
        out_b0_dwi_merge,
        out_updated_bval,
        out_updated_bvec,
    ] = toimportfunctions.insert_b0_into_dwi(
        in_b0=out_reference_b0,
        in_dwi=out_split_dwi,
        in_bval=out_split_bval,
        in_bvec=out_split_bvec,
    )
    return out_reference_b0, out_b0_dwi_merge, out_updated_bval, out_updated_bvec


@task
@annotate({"return": {"out_acq": str}})
def generate_acq_file(
    in_dwi, fsl_phase_encoding_direction, total_readout_time, image_id=None
):
    """Generate [`image_id`]_acq.txt file for FSL eddy command.

    Args:
        in_dwi (str): DWI file.
        fsl_phase_encoding_direction (str): PhaseEncodingDirection from BIDS specifications in FSL format (i.e. x/y/z instead of i/j/k).
        total_readout_time (str): TotalReadoutTime from BIDS specifications.
        image_id (str, optional): Optional prefix. Defaults to None.

    Returns:
        out_acq: [`image_id`]_acq.txt or acq.txt file.
    """
    import os

    import nibabel as nb
    import numpy as np

    if image_id:
        out_acq = os.path.abspath(f"{image_id}_acq.txt")
    else:
        out_acq = os.path.abspath("acq.txt")
    vols = nb.load(in_dwi).get_fdata().shape[-1]
    arr = np.ones([vols, 4])
    for i in range(vols):
        if fsl_phase_encoding_direction == "y-":
            arr[i, :] = np.array((0, -1, 0, total_readout_time))
        elif fsl_phase_encoding_direction == "y":
            arr[i, :] = np.array((0, 1, 0, total_readout_time))
        elif fsl_phase_encoding_direction == "x":
            arr[i, :] = np.array((1, 0, 0, total_readout_time))
        elif fsl_phase_encoding_direction == "x-":
            arr[i, :] = np.array((-1, 0, 0, total_readout_time))
        elif fsl_phase_encoding_direction == "z":
            arr[i, :] = np.array((0, 1, 0, total_readout_time))
        elif fsl_phase_encoding_direction == "z-":
            arr[i, :] = np.array((0, 0, -1, total_readout_time))
        else:
            raise RuntimeError(
                f"FSL PhaseEncodingDirection (found value: {fsl_phase_encoding_direction}) "
                f"is unknown, it should be a value in (x, y, z, x-, y-, z-)"
            )

    np.savetxt(out_acq, arr, fmt="%d " * 3 + "%f")

    return out_acq


@task
@annotate({"return": {"out_index": str}})
def generate_index_file(in_bval, low_bval=5.0, image_id=None):
    """Generate [`image_id`]_index.txt file for FSL eddy command.

    Args:
        in_bval (str): Bval file.
        low_bval (float): Define the b0 volumes as all volume bval <= low_bval. Default to 5.0.
        image_id (str, optional): Optional prefix. Defaults to None.

    Returns:
        out_index: [`image_id`]_index.txt or index.txt file.
    """
    import os

    import numpy as np

    assert os.path.isfile(in_bval)
    bvals = np.loadtxt(in_bval)
    idx_low_bvals = np.where(bvals <= low_bval)
    b0_index = idx_low_bvals[0].tolist()

    if not b0_index:
        raise ValueError(
            f"Could not find b-value <= {low_bval} in bval file ({in_bval}). Found values: {bvals}"
        )

    if image_id:
        out_index = os.path.abspath(f"{image_id}_index.txt")
    else:
        out_index = os.path.abspath("index.txt")

    vols = len(bvals)
    index_list = []
    for i in range(0, len(b0_index)):
        if i == (len(b0_index) - 1):
            index_list.extend([i + 1] * (vols - b0_index[i]))
        else:
            index_list.extend([i + 1] * (b0_index[i + 1] - b0_index[i]))
    index_array = np.asarray(index_list)
    try:
        len(index_list) == vols
    except ValueError:
        raise ValueError(
            "It seems that you do not define the index file for FSL eddy correctly!"
        )
    np.savetxt(out_index, index_array.T)

    return out_index


@task
@annotate({"return": {"out_matrix_list": str}})
def expend_matrix_list(in_matrix, in_bvec):
    import numpy as np

    bvecs = np.loadtxt(in_bvec).T
    out_matrix_list = [in_matrix]

    out_matrix_list = out_matrix_list * len(bvecs)

    return out_matrix_list


@task
@annotate({"return": {"out_file": str}})
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


@task
@annotate({"return": {"updated_affine_file": str}})
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


@task
@annotate({"return": {"warped_image": str}})
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


@task
@annotate({"return": {"out_b0_average": str}})
def compute_average_b0(in_dwi: str, in_bval: str, low_bval: float = 5.0):
    """Compute average b0 volume from DWI dataset."""
    import os

    import nibabel
    import numpy as np

    assert os.path.isfile(in_dwi)
    assert os.path.isfile(in_bval)
    assert low_bval >= 0

    imgs = np.array(nibabel.four_to_three(nibabel.load(in_dwi)))
    bvals = np.loadtxt(in_bval)
    low_bvals = np.where(bvals <= low_bval)[0]

    fname_b0, ext_b0 = os.path.splitext(os.path.basename(in_dwi))
    if ext_b0 == ".gz":
        fname_b0, ext2 = os.path.splitext(fname_b0)
        ext_b0 = ext2 + ext_b0
    out_b0_average = os.path.abspath(f"{fname_b0}_avg_b0{ext_b0}")

    b0s_data = [imgs[i].get_data() for i in low_bvals]
    avg_b0_data = np.average(np.array(b0s_data), axis=0)

    hdr = imgs[0].get_header().copy()
    hdr.set_data_shape(avg_b0_data.shape)
    nibabel.Nifti1Image(avg_b0_data, imgs[0].get_affine(), hdr).to_filename(
        out_b0_average
    )

    return out_b0_average


@task
@annotate({"return": {"out_file": str}})
def rotate_bvecs_alt(in_bvec, in_matrix):
    """Rotate the input bvec file accordingly with a list of matrices.

    Notes:
        The input affine matrix transforms points in the destination image to their corresponding
        coordinates in the original image. Therefore, this matrix should be inverted first, as
        we want to know the target position of :math:`\\vec{r}`.

    """
    import os

    import numpy as np

    print("in_bvec type: ", type(in_bvec[0]))
    print("in_bvec: ", in_bvec)
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
