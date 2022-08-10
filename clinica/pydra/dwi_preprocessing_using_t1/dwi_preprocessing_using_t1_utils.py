def merge_volumes_tdim(in_file1, in_file2):
    """Merge 'in_file1' and 'in_file2' in the t dimension.

    Args:
        in_file1 (str): First set of volumes.
        in_file2 (str): Second set of volumes.

    Returns:
        out_file (str): The two sets of volumes merged.
    """
    import os

    out_file = os.path.abspath("merged_files.nii.gz")
    cmd = f"fslmerge -t {out_file} {in_file1} {in_file2}"
    os.system(cmd)
    return out_file


####################################################
def bids_dir_to_fsl_dir(bids_dir):
    """Converts BIDS PhaseEncodingDirection parameters (i,j,k,i-,j-,k-) to FSL direction (x,y,z,x-,y-,z-)."""
    fsl_dir = bids_dir.lower()
    if "i" not in fsl_dir and "j" not in fsl_dir and "k" not in fsl_dir:
        raise ValueError(
            f"Unknown PhaseEncodingDirection {fsl_dir}: it should be a value in (i, j, k, i-, j-, k-)"
        )

    return fsl_dir.replace("i", "x").replace("j", "y").replace("k", "z")


def check_dwi_volume(in_dwi, in_bvec, in_bval):
    """Check that # DWI = # B-val = # B-vec.

    Raises
        ValueError if # DWI, # B-val and # B-vec mismatch.
    """
    import nibabel as nib
    import numpy as np

    bvals = np.loadtxt(in_bval)
    num_b_vals = len(bvals)

    bvecs = np.loadtxt(in_bvec)
    _, num_b_vecs = bvecs.shape

    img = nib.load(in_dwi)
    _, _, _, num_dwis = img.shape

    if not (num_b_vals == num_b_vecs == num_dwis):
        raise IOError(
            f"Number of DWIs, b-vals and b-vecs mismatch "
            f"(# DWI = {num_dwis}, # B-vec = {num_b_vecs}, #B-val = {num_b_vals}) "
        )


####################################################
def extract_metadata_from_json(json_file, list_keys):
    """Extract fields from JSON file."""
    import datetime
    import json

    list_values = []
    try:
        with open(json_file, "r") as file:
            data = json.load(file)
            for key in list_keys:
                list_values.append(data[key])
    except EnvironmentError:
        raise EnvironmentError(
            f"[Error] Clinica could not open the following JSON file: {json_file}"
        )
    except KeyError as e:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_message = f"[{now}] Error: Clinica could not find the e key in the following JSON file: {json_file}"
        raise EnvironmentError(error_message)
    finally:
        file.close()

    return list_values


####################################################
def get_subject_id(bids_or_caps_file: str) -> str:
    """Extract "sub-<participant_id>_ses-<session_label>" from BIDS or CAPS file."""
    import re

    m = re.search(r"(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)", bids_or_caps_file)

    if not m:
        raise ValueError(
            f"Input filename {bids_or_caps_file} is not in a BIDS or CAPS compliant format."
            " It does not contain the subject and session information."
        )

    subject_id = m.group(1) + "_" + m.group(2)

    return subject_id


def print_begin_image(image_id, list_keys=None, list_values=None):
    """Print begin run pipeline message for a given image `image_id`."""

    if list_keys is not None:
        assert len(list_keys) == len(list_values)

    begin_message = f"Running pipeline for {image_id.replace('_', ' | ')}"
    if list_keys and list_values:
        begin_message += " ("
        begin_message += ", ".join(
            f"{key} = {key_value}" for key, key_value in zip(list_keys, list_values)
        )
        begin_message += ")"
    print(f"{begin_message}")


def b0_average(in_file, out_file=None):
    """
    Average the b0 volumes.

    Args:
        in_file (str): The b0 volumes already registered.
        out_file (str, optional): Name of the file. Defaults to None.

    Returns:
        The mean of the b0 volumes.

    Warnings:
        The b0 volumes must be registered.
    """
    import os

    import nibabel as nb
    import numpy as np

    if not out_file:
        fname, ext = os.path.splitext(os.path.basename(in_file))
        if ext == ".gz":
            fname, ext2 = os.path.splitext(fname)
            ext = ext2 + ext
        out_file = os.path.abspath(f"{fname}_avg_b0{ext}")

    imgs = np.array(nb.four_to_three(nb.load(in_file)))
    b0s = [im.get_fdata(dtype="float32").astype(np.float32) for im in imgs]
    b0 = np.average(np.array(b0s), axis=0)

    hdr = imgs[0].get_header().copy()
    hdr.set_data_shape(b0.shape)
    hdr.set_xyzt_units("mm")
    hdr.set_data_dtype(np.float32)
    nb.Nifti1Image(b0, imgs[0].get_affine(), hdr).to_filename(out_file)

    return out_file


def b0_dwi_split(in_dwi, in_bval, in_bvec, low_bval=5.0):
    """Split DWI dataset.

    Split the DWI volumes into two datasets :
     - the first dataset contains the set of b<=low_bval volumes.
     - the second dataset contains the set of DWI volumes.

    Args:
        in_dwi (str): DWI dataset.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.
        low_bval (int, optional): Define the b0 volumes as all volume bval <= lowbval. Defaults to 5.0.

    Returns:
        out_b0 (str): The set of b<=low_bval volumes.
        out_dwi (str): Output. The set of b>low_bval volumes.
        out_bvals (str): The b-values corresponding to the out_dwi.
        out_bvecs (str): The b-vecs corresponding to the out_dwi.
    """
    import os
    import warnings

    import nibabel as nib
    import numpy as np

    assert os.path.isfile(in_dwi)
    assert os.path.isfile(in_bval)
    assert os.path.isfile(in_bvec)
    assert low_bval >= 0

    im = nib.load(in_dwi)
    data = im.get_fdata(dtype="float32")
    hdr = im.get_header().copy()
    bvals = np.loadtxt(in_bval)
    bvecs = np.loadtxt(in_bvec)

    if bvals.shape[0] == bvecs.shape[0]:
        warnings.warn(
            "Warning: The b-vectors file should be column-wise. The b-vectors will be transposed",
            UserWarning,
        )
        bvecs = bvecs.T

    lowbs = np.where(bvals <= low_bval)[0]

    fname_b0, ext_b0 = os.path.splitext(os.path.basename(in_dwi))
    if ext_b0 == ".gz":
        fname_b0, ext2 = os.path.splitext(fname_b0)
        ext_b0 = ext2 + ext_b0
    out_b0 = os.path.abspath(f"{fname_b0}_b0{ext_b0}")
    # out_b0 = op.abspath('b0.nii.gz')
    b0data = data[..., lowbs]
    hdr.set_data_shape(b0data.shape)
    nib.Nifti1Image(b0data, im.get_affine(), hdr).to_filename(out_b0)

    dwi_bvals = np.where(bvals > low_bval)[0]
    out_dwi = os.path.abspath("dwi.nii.gz")
    dwi_data = data[..., dwi_bvals]
    hdr.set_data_shape(dwi_data.shape)
    nib.Nifti1Image(dwi_data, im.get_affine(), hdr).to_filename(out_dwi)

    bvals_dwi = bvals[dwi_bvals]
    out_bvals = os.path.abspath("bvals")
    np.savetxt(out_bvals, bvals_dwi, fmt="%d", delimiter=" ")

    bvecs_dwi = np.array(
        [
            bvecs[0][dwi_bvals].tolist(),
            bvecs[1][dwi_bvals].tolist(),
            bvecs[2][dwi_bvals].tolist(),
        ]
    )
    out_bvecs = os.path.abspath("bvecs")
    np.savetxt(out_bvecs, bvecs_dwi, fmt="%10.5f", delimiter=" ")

    return out_b0, out_dwi, out_bvals, out_bvecs


####################################################
def count_b0s(in_bval, low_bval=5.0):
    """Count the number of volumes where b<=low_bval.

    Args:
        in_bval (str): bval file.
        low_bval (int, optional): Define the b0 volumes as all volume bval <= lowbval. Defaults to 5.0.

    Returns:
        num_b0s: Number of b0s.
    """
    import numpy as np

    bvals = np.loadtxt(in_bval)
    print("bvals: ", bvals)
    num_b0s = len(np.where(bvals <= low_bval)[0])

    return num_b0s


def insert_b0_into_dwi(in_b0, in_dwi, in_bval, in_bvec):
    """Insert a b0 volume into the dwi dataset as the first volume and update the bvals and bvecs files.

    Args:
        in_b0 (str): One b=0 volume (could be the average of a b0 dataset).
        in_dwi (str): The set of DWI volumes.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.

    Returns:
        out_dwi (str): Diffusion dataset : b0 volume + dwi volumes.
        out_bval (str): B-values update.
        out_bvec (str): Directions of diffusion update.
    """
    import os

    import numpy as np

    assert os.path.isfile(in_b0)
    assert os.path.isfile(in_dwi)
    assert os.path.isfile(in_bval)
    assert os.path.isfile(in_bvec)

    out_dwi = merge_volumes_tdim(in_b0, in_dwi)

    lst = np.loadtxt(in_bval).tolist()
    lst.insert(0, 0)
    out_bvals = os.path.abspath("bvals")
    np.savetxt(out_bvals, np.matrix(lst), fmt="%d", delimiter=" ")

    bvecs = np.loadtxt(in_bvec)
    bvecs_0 = bvecs[0].tolist()
    bvecs_0.insert(0, 0.0)
    bvecs_1 = bvecs[1].tolist()
    bvecs_1.insert(0, 0.0)
    bvecs_2 = bvecs[2].tolist()
    bvecs_2.insert(0, 0.0)
    bvecs_dwi = np.array([bvecs_0, bvecs_1, bvecs_2])
    out_bvecs = os.path.abspath("bvecs")
    np.savetxt(out_bvecs, bvecs_dwi, fmt="%10.5f", delimiter=" ")

    return out_dwi, out_bvals, out_bvecs


def build_flirt_workflow(nb_b0s, extracted_b0, working_dir):
    # imports
    from pydra import Workflow

    # Main workflow: workflows will be chained inside or that is the plan
    workflow = Workflow(
        name="complete_wf",
        input_spec=["num_b0s", "in_file", "base_dir"],
        num_b0s=nb_b0s,
        in_file=extracted_b0,
        base_dir=working_dir,
    )

    ref = ref_wf()
    ref.inputs.num_b0s = workflow.lzin.num_b0s
    ref.inputs.in_file = workflow.lzin.in_file
    ref.inputs.base_dir = workflow.lzin.base_dir
    workflow.add(ref)

    moving = moving_wf(nb_b0s, extracted_b0, working_dir)
    workflow.add(moving)

    b0_flirt_workflow = flirt_wf()
    b0_flirt_workflow.inputs.split_moving = moving.lzout.split_moving
    b0_flirt_workflow.inputs.dilate = ref.lzout.ref_weight
    b0_flirt_workflow.inputs.fslroi_ref = ref.lzout.reference
    workflow.add(b0_flirt_workflow)

    workflow.set_output(
        [
            ("out_xmf", b0_flirt_workflow.lzout.out_matrix_file),
            ("out_file", b0_flirt_workflow.lzout.out_file),
        ]
    )
    return workflow


def run_b0_flirt_workflow(workflow):
    from pydra import Submitter, Workflow
    from pydra.tasks.nipype1.utils import Nipype1Task

    with Submitter(plugin="cf") as submitter:
        submitter(workflow)

    return workflow.result(return_inputs=True)


def ref_wf():
    from nipype.interfaces.fsl.maths import MathsCommand
    from nipype.interfaces.fsl.preprocess import BET
    from nipype.interfaces.fsl.utils import ExtractROI
    from pydra import Workflow
    from pydra.tasks.nipype1.utils import Nipype1Task

    fslroi_ref = ExtractROI(args="0 1")
    bet_ref = BET(frac=0.3, robust=True)
    dilate = MathsCommand(nan2zeros=True, args="-kernel sphere 5 -dilM")

    ref = Workflow(
        name="ref",
        input_spec=["num_b0s", "in_file", "base_dir"],
    )
    ref.add(
        Nipype1Task(
            name="fslroi_ref",
            interface=fslroi_ref,
            in_file=ref.lzin.in_file,
        )
    )
    ref.add(
        Nipype1Task(
            name="bet_ref",
            interface=bet_ref,
            in_file=ref.fslroi_ref.lzout.roi_file,
            mask=True,
        )
    )
    ref.add(
        Nipype1Task(
            name="dilate",
            interface=dilate,
            in_file=ref.bet_ref.lzout.mask_file,
        )
    )
    ref.set_output(
        [
            ("ref_weight", ref.dilate.lzout.out_file),
            ("reference", ref.fslroi_ref.lzout.roi_file),
            ("bet_ref", ref.bet_ref.lzout.mask_file),
        ]
    )
    return ref


def run_ref_wf(ref_wf):
    from pydra import Submitter

    with Submitter(plugin="cf") as submitter:
        submitter(ref_wf)

    results_ref = ref_wf.result(return_inputs=True)

    return results_ref


def moving_wf(nb_b0s, extracted_b0, working_dir):
    from nipype.interfaces.fsl.utils import ExtractROI, Split
    from pydra import Workflow
    from pydra.tasks.nipype1.utils import Nipype1Task

    tsize = nb_b0s - 1
    fslroi_moving = ExtractROI(args=("1 " + str(tsize)))
    split_moving = Split(dimension="t")
    moving = Workflow(
        name="moving",
        input_spec=["num_b0s", "in_file", "base_dir"],
        num_b0s=nb_b0s,
        in_file=extracted_b0,
        base_dir=working_dir,
        # dilate=results_ref[1].get_output_field('ref_weight'),
        # fslroi_ref=results_ref[1].get_output_field('reference'),
    )
    moving.add(
        Nipype1Task(
            name="fslroi_moving",
            interface=fslroi_moving,
            in_file=moving.lzin.in_file,
        )
    )
    moving.add(
        Nipype1Task(
            name="split_moving",
            interface=split_moving,
            in_file=moving.fslroi_moving.lzout.roi_file,
        )
    )
    moving.set_output(
        [
            ("split_moving", moving.split_moving.lzout.out_files),
            # ("out_matrix", b0_moving_workflow.flirt.lzout.out_matrix_file),
        ]
    )
    return moving


def run_moving_wf(moving_wf):
    from pydra import Submitter

    with Submitter(plugin="cf") as submitter:
        submitter(moving_wf)

    results_moving = moving_wf.result(return_inputs=True)

    return results_moving


def flirt_wf():
    from nipype.interfaces.fsl.maths import Threshold
    from nipype.interfaces.fsl.preprocess import FLIRT
    from nipype.interfaces.fsl.utils import Merge
    from pydra import Workflow
    from pydra.mark import annotate, task
    from pydra.tasks.nipype1.utils import Nipype1Task

    flirt = FLIRT(
        interp="spline",
        dof=6,
        bins=50,
        save_log=True,
        cost="corratio",
        cost_func="corratio",
        padding_size=10,
        searchr_x=[-4, 4],
        searchr_y=[-4, 4],
        searchr_z=[-4, 4],
        fine_search=1,
        coarse_search=10,
    )
    merge = Merge(dimension="t")
    thres = Threshold(thresh=0.0)

    @task
    @annotate({"return": {"out_file": str}})
    def merge_volumes_tdim(in_file1, in_file2):
        """Merge 'in_file1' and 'in_file2' in the t dimension.

        Args:
            in_file1 (str): First set of volumes.
            in_file2 (str): Second set of volumes.

        Returns:
            out_file (str): The two sets of volumes merged.
        """
        import os

        out_file = os.path.abspath("merged_files.nii.gz")
        cmd = f"fslmerge -t {out_file} {in_file1} {in_file2}"
        os.system(cmd)
        return out_file

    b0_flirt_workflow = Workflow(
        name="b0_flirt_workflow",
        input_spec=["split_moving", "dilate", "fslroi_ref"],
    )

    b0_flirt_workflow.add(
        Nipype1Task(
            name="flirt",
            interface=flirt,
            in_file=b0_flirt_workflow.lzin.split_moving,
            ref_weight=b0_flirt_workflow.lzin.dilate,
            in_weight=b0_flirt_workflow.lzin.dilate,
            reference=b0_flirt_workflow.lzin.fslroi_ref,
        ).split("in_file")
    )
    # b0_flirt_workflow.flirt.split("split_moving")
    b0_flirt_workflow.add(
        Nipype1Task(
            name="thres",
            interface=thres,
            in_file=b0_flirt_workflow.flirt.lzout.out_file,
        )
    )
    b0_flirt_workflow.thres.combine("flirt.in_file")
    b0_flirt_workflow.add(
        Nipype1Task(
            name="merge",
            interface=merge,
            in_files=b0_flirt_workflow.thres.lzout.out_file,
        )
    )
    b0_flirt_workflow.add(
        merge_volumes_tdim(
            name="merge_volumes_tdim",
            interface=merge_volumes_tdim,
            in_file1=b0_flirt_workflow.merge.lzout.merged_file,
            in_file2=b0_flirt_workflow.lzin.fslroi_ref,
        )
    )
    b0_flirt_workflow.set_output(
        [
            ("out_matrix_file", b0_flirt_workflow.flirt.lzout.out_matrix_file),
            ("out_file", b0_flirt_workflow.merge_volumes_tdim.lzout.out_file),
        ]
    )

    return b0_flirt_workflow


def run_flirt_wf(flirt_wf):
    from pydra import Submitter

    with Submitter(plugin="cf") as submitter:
        submitter(flirt_wf)

    results_flirt = flirt_wf.result(return_inputs=True)

    return results_flirt


def build_n_run_flirt(nb_b0s, extracted_b0, working_dir):

    ref = ref_wf()
    ref.inputs.num_b0s = nb_b0s
    ref.inputs.in_file = extracted_b0
    ref.inputs.base_dir = working_dir
    ref_results = run_ref_wf(ref)

    moving = moving_wf(nb_b0s, extracted_b0, working_dir)
    moving_results = run_moving_wf(moving)

    b0_flirt_workflow = flirt_wf()
    b0_flirt_workflow.inputs.split_moving = moving_results[1].get_output_field(
        "split_moving"
    )
    b0_flirt_workflow.inputs.dilate = ref_results[1].get_output_field("ref_weight")
    b0_flirt_workflow.inputs.fslroi_ref = ref_results[1].get_output_field("reference")
    flirt_results = run_flirt_wf(b0_flirt_workflow)

    return flirt_results


def build_epi_wf(t1w, out_updated_bvec, out_corrected):
    import dwi_preprocessing_using_t1_tasks as dt
    import nipype.interfaces.utility as niu
    from nipype.interfaces.ants import CreateJacobianDeterminantImage
    from nipype.interfaces.ants.registration import RegistrationSynQuick
    from nipype.interfaces.c3 import C3dAffineTool
    from nipype.interfaces.fsl import MultiImageMaths
    from nipype.interfaces.fsl.epi import Eddy
    from nipype.interfaces.fsl.maths import Threshold
    from nipype.interfaces.fsl.preprocess import BET, FLIRT
    from nipype.interfaces.fsl.utils import Merge, Split
    from nipype.interfaces.io import BIDSDataGrabber
    from nipype.interfaces.mrtrix3 import DWIBiasCorrect
    from nipype.interfaces.mrtrix3.preprocess import DWIBiasCorrect
    from pydra import Workflow
    from pydra.tasks.nipype1.utils import Nipype1Task

    bet_2 = BET(frac=0.3, robust=True)
    use_cuda = False
    initrand = False
    eddy = Eddy(repol=True, use_cuda=use_cuda, initrand=initrand)
    split_epi = Split(dimension="t")
    select_epi = niu.Select()
    flirt_b0_2_t1 = FLIRT(dof=6, interp="spline", cost="normmi", cost_func="normmi")
    ants_reg = RegistrationSynQuick(transform_type="br", dimension=3)
    c3d_affine_tool = C3dAffineTool(itk_transform=True, fsl2ras=True)
    merge_transform = niu.Merge(3)
    jacobian = CreateJacobianDeterminantImage(
        imageDimension=3, outputImage="Jacobian_image.nii.gz"
    )
    jacmult = MultiImageMaths(op_string="-mul %s")
    thres_epi = Threshold(thresh=0.0)
    merge_epi = Merge(dimension="t")

    working_dir = "/localdrive10TB/users/matthieu.joulot/werk"
    workflow = Workflow(
        name="epi",
        input_spec=["t1w", "out_updated_bvec", "out_corrected"],
        t1w=t1w,
        out_updated_bvec=out_updated_bvec,
        out_corrected=out_corrected,
    )
    workflow.add(
        Nipype1Task(
            name="split_epi",
            interface=split_epi,
            # in_file=eddy_out_corrected,
            # in_file=workflow.lzin.out_corrected,
            in_file=workflow.lzin.t1w,
        )
    )
    # workflow2.split("split_epi.lzout.out_files")
    workflow.add(
        Nipype1Task(
            name="select_epi",
            interface=select_epi,
            inlist=workflow.split_epi.lzout.out_files,
            index=[0],
        )
    )
    workflow.add(
        Nipype1Task(
            name="flirt_b0_2_t1",
            interface=flirt_b0_2_t1,
            # in_file=workflow.init_input_node.lzout.t1w,
            in_file=workflow.select_epi.lzout.out,
            reference=workflow.lzin.t1w,
        )
    )

    workflow.add(
        dt.expend_matrix_list(
            name="expend_matrix_list",
            # interface=expend_matrix_list,
            in_matrix=workflow.flirt_b0_2_t1.lzout.out_matrix_file,
            # in_bvec=workflow.lzin.out_rotated_bvecs,
            in_bvec=workflow.lzin.out_updated_bvec,
        )
    )
    workflow.add(
        dt.rotate_bvecs(
            name="rotate_bvecs",
            # interface=rotate_bvecs,
            # in_bvec=workflow.lzin.out_rotated_bvecs,
            in_bvec=workflow.lzin.out_updated_bvec,
            in_matrix=workflow.expend_matrix_list.lzout.out_matrix_list,
        )
    )
    workflow.add(
        Nipype1Task(
            name="ants_reg",
            interface=ants_reg,
            moving_image=workflow.flirt_b0_2_t1.lzout.out_file,
            # moving_image = workflow.lzin.out_corrected,
            # fixed_image=workflow.init_input_node.lzout.t1w,
            fixed_image=workflow.lzin.t1w,
            # fixed_image=workflow.printer.lzout.var_out,
        )
    )
    workflow.add(
        Nipype1Task(
            name="c3d_affine_tool",
            interface=c3d_affine_tool,
            source_file=workflow.select_epi.lzout.out,
            reference_file=workflow.lzin.t1w,
            transform_file=workflow.flirt_b0_2_t1.lzout.out_matrix_file,
        )
    )
    workflow.add(
        dt.change_itk_transform_type(
            name="change_itk_transform_type",
            input_affine_file=workflow.c3d_affine_tool.lzout.itk_transform,
        )
    )
    workflow.add(
        Nipype1Task(
            name="merge_transform",
            interface=merge_transform,
            in1=workflow.change_itk_transform_type.lzout.updated_affine_file,
            in2=workflow.ants_reg.lzout.out_matrix,
            in3=workflow.ants_reg.lzout.forward_warp_field,
        )
    )

    workflow.add(
        dt.ants_apply_transforms(
            name="apply_transform_image",
            fixed_image=workflow.lzin.t1w,
            moving_image=workflow.split_epi.lzout.out_files,
            transforms=workflow.merge_transform.lzout.out,
            warped_image="out_warped_field.nii.gz",
            output_warped_image=True,
        ).split("moving_image")
    )
    workflow.add(
        dt.ants_apply_transforms(
            name="apply_transform_field",
            # interface=ants_apply_transforms,
            fixed_image=workflow.lzin.t1w,
            moving_image=workflow.split_epi.lzout.out_files,
            transforms=workflow.merge_transform.lzout.out,
            warped_image="out_warped.nii.gz",
            output_warped_image=False,
        ).split(splitter=("moving_image"))
    )
    workflow.add(
        Nipype1Task(
            name="jacobian",
            interface=jacobian,
            deformationField=workflow.apply_transform_field.lzout.warped_image,
        )
    )
    workflow.jacobian.combine(["apply_transform_field.moving_image"])
    workflow.apply_transform_image.combine(["moving_image"])
    workflow.add(
        Nipype1Task(
            name="jacmult",
            interface=jacmult,
            operand_files=workflow.apply_transform_image.lzout.warped_image,
            in_file=workflow.jacobian.lzout.jacobian_image,
        ).split(("operand_files", "in_file"))
    )
    workflow.add(
        Nipype1Task(
            name="thres_epi",
            interface=thres_epi,
            in_file=workflow.jacmult.lzout.out_file,
        )
    )
    workflow.thres_epi.combine(["jacmult.operand_files", "jacmult.in_file"])

    workflow.add(
        Nipype1Task(
            name="merge_epi",
            interface=merge_epi,
            in_files=workflow.thres_epi.lzout.out_file,
        )
    )
    # ## epi pipeline end
    workflow.set_output(
        [
            ("trans", workflow.apply_transform_image.lzout.warped_image),
            ("merge_epi", workflow.merge_epi.lzout.merged_file),
            ("rotated_bvecs", workflow.rotate_bvecs.lzout.out_file),
        ]
    )
    return workflow


def build_eddy_wf(
    out_b0_dwi_merge,
    phase_encoding_direction,
    total_readout_time,
    out_updated_bval,
    mask_file,
    out_updated_bvec,
):
    import dwi_preprocessing_using_t1_tasks as dt
    from nipype.interfaces.fsl.epi import Eddy
    from pydra import Workflow
    from pydra.tasks.nipype1.utils import Nipype1Task

    use_cuda = False
    initrand = False
    eddy = Eddy(repol=True, use_cuda=use_cuda, initrand=initrand)
    eddy_wf = Workflow(
        name="eddyy",
        input_spec=[
            "out_b0_dwi_merge",
            "phase_encoding_direction",
            "total_readout_time",
            "out_updated_bval",
            "mask_file",
            "out_updated_bvec",
        ],
        out_b0_dwi_merge=out_b0_dwi_merge,
        phase_encoding_direction=phase_encoding_direction,
        total_readout_time=total_readout_time,
        out_updated_bval=out_updated_bval,
        mask_file=mask_file,
        out_updated_bvec=out_updated_bvec,
    )
    ## Eddy pipeline parts: generate_acq, generate_index, eddy
    eddy_wf.add(
        dt.generate_acq_file(
            name="generate_acq_file",
            # interface=generate_acq_file,
            in_dwi=eddy_wf.lzin.out_b0_dwi_merge,
            fsl_phase_encoding_direction=eddy_wf.lzin.phase_encoding_direction,
            total_readout_time=eddy_wf.lzin.total_readout_time,
            image_id=None,
        )
    )
    eddy_wf.add(
        dt.generate_index_file(
            name="generate_index_file",
            # interface=generate_index_file,
            in_bval=eddy_wf.lzin.out_updated_bval,
            low_bval=5.0,
            image_id=None,
        )
    )
    # eddy_wf.add(
    #     Nipype1Task(
    #         name="eddy",
    #         interface=eddy,
    #         in_file=eddy_wf.lzin.out_b0_dwi_merge,
    #         in_mask=eddy_wf.lzin.mask_file,
    #         in_bval=eddy_wf.lzin.out_updated_bval,
    #         in_bvec=eddy_wf.lzin.out_updated_bvec,
    #         in_acqp=eddy_wf.generate_acq_file.lzout.out_acq,
    #         in_index=eddy_wf.generate_index_file.lzout.out_index,
    #     )
    # )
    eddy_wf.set_output(
        [
            ("out_acq", eddy_wf.generate_acq_file.lzout.out_acq),
            # ("out_rotated", eddy_wf.eddy.lzout.out_rotated_bvecs),
            # ("out_corrected", ),
        ]
    )
    return eddy_wf
