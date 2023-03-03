from nipype.pipeline.engine import Workflow


def eddy_fsl_pipeline(
    use_cuda: bool,
    initrand: bool,
    compute_mask: bool = True,
    image_id: bool = False,
    field: bool = False,
    output_dir=None,
    name: str = "eddy_fsl",
) -> Workflow:
    """Build a FSL-based pipeline for head motion correction and eddy current distortion correction.

    This pipeline needs FSL as it relies on the Eddy executable.

    Parameters
    ----------
    use_cuda : bool
        Boolean to indicate whether cuda should be used or not.

    initrand : bool
        ????

    compute_mask : bool, optional
        If True, the eddy_fsl pipeline computes the brain mask
        of the reference B0 volume using FSL BET. And provide the
        result to Eddy.
        If False, a mask should be provided to the eddy_fsl
        pipeline through the "in_mask" connection.
        Default=True.

    image_id : bool, optional
        Boolean to indicate whether the Eddy node should expect
        a value for its 'out_base' parameter. This is used for
        building the names of the output files.
        Default=False.

    field : bool, optional
        Boolean to indicate whether the Eddy node should expect
        a value for its 'field' input.
        Default=False.

    output_dir: str, optional
        Path to output directory.
        If provided, the pipeline will write its output in this folder.
        Default to None.

    name : str, optional
        Name of the pipeline. Default='eddy_fsl'.

    Returns
    -------
    Workflow :
        The Nipype workflow.
        This workflow has the following inputs:
            - "dwi_filename": The path to the DWI image
            - "b_vectors_filename": The path to the associated B-vectors file
            - "b_values_filename": The path to the associated B-values file
            - "in_mask": The path to the mask image to be provided to Eddy
              This will be used only if compute_mask is False.
            - "image_id": Prefix to be used for output files
            - "field": The path to the field image to be used by Eddy
            - "reference_b0": The path to the reference B0 image
              This will be used to compute the mask only if compute_mask is True.
            - "total_readout_time": The total readout time extracted from JSON metadata
            - "phase_encoding_direction": The phase encoding direction extracted from JSON metadata

        And the following outputs:
            - "out_parameter": Path to the file storing the output parameters
            - "out_corrected": Path to the corrected image
            - "out_rotated_bvecs": Path to the file holding the rotated B-vectors
    """
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces.fsl import BET
    from nipype.interfaces.fsl.epi import Eddy

    from clinica.utils.dwi import generate_acq_file, generate_index_file

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "dwi_filename",
                "b_vectors_filename",
                "b_values_filename",
                "in_mask",
                "image_id",
                "field",
                "reference_b0",
                "total_readout_time",
                "phase_encoding_direction",
            ]
        ),
        name="inputnode",
    )

    mask_reference_b0 = pe.Node(
        BET(frac=0.3, mask=True, robust=True), name="mask_reference_b0"
    )

    generate_acq = pe.Node(
        niu.Function(
            input_names=[
                "dwi_filename",
                "fsl_phase_encoding_direction",
                "total_readout_time",
                "image_id",
            ],
            output_names=["out_file"],
            function=generate_acq_file,
        ),
        name="generate_acq",
    )

    generate_index = pe.Node(
        niu.Function(
            input_names=["b_values_filename"],
            output_names=["out_file"],
            function=generate_index_file,
        ),
        name="generate_index",
    )

    eddy = pe.Node(interface=Eddy(), name="eddy_fsl")
    eddy.inputs.repol = True
    eddy.inputs.use_cuda = use_cuda
    eddy.inputs.initrand = initrand

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["out_parameter", "out_corrected", "out_rotated_bvecs"]
        ),
        name="outputnode",
    )

    if output_dir:
        write_results = pe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = output_dir
        write_results.inputs.parameterization = False

    wf = pe.Workflow(name=name)

    connections = [
        (inputnode, generate_acq, [("dwi_filename", "dwi_filename")]),
        (inputnode, generate_acq, [("total_readout_time", "total_readout_time")]),
        (
            inputnode,
            generate_acq,
            [("phase_encoding_direction", "fsl_phase_encoding_direction")],
        ),
        (inputnode, generate_index, [("b_values_filename", "b_values_filename")]),
        (inputnode, generate_index, [("image_id", "image_id")]),
        (inputnode, eddy, [("b_vectors_filename", "in_bvec")]),
        (inputnode, eddy, [("b_values_filename", "in_bval")]),
        (inputnode, eddy, [("dwi_filename", "in_file")]),
        (generate_acq, eddy, [("out_file", "in_acqp")]),
        (generate_index, eddy, [("out_file", "in_index")]),
        (eddy, outputnode, [("out_parameter", "out_parameter")]),
        (eddy, outputnode, [("out_corrected", "out_corrected")]),
        (eddy, outputnode, [("out_rotated_bvecs", "out_rotated_bvecs")]),
    ]

    if image_id:
        connections += [
            (inputnode, generate_acq, [("image_id", "image_id")]),
            (inputnode, eddy, [("image_id", "out_base")]),
        ]

    if field:
        connections += [(inputnode, eddy, [("field", "field")])]

    if compute_mask:
        connections += [
            (inputnode, mask_reference_b0, [("reference_b0", "in_file")]),
            (mask_reference_b0, eddy, [("mask_file", "in_mask")]),
        ]
    else:
        connections += [(inputnode, eddy, [("in_mask", "in_mask")])]

    if output_dir:
        connections += [
            (outputnode, write_results, [("out_parameter", "out_parameter")]),
            (outputnode, write_results, [("out_corrected", "out_corrected")]),
            (outputnode, write_results, [("out_rotated_bvecs", "out_rotated_bvecs")]),
        ]

    wf.connect(connections)

    return wf


def epi_pipeline(
    base_dir: str,
    delete_cache: bool = False,
    name="susceptibility_distortion_correction_using_t1",
):
    """Perform EPI correction.

    This workflow allows to correct for echo-planar induced susceptibility artifacts without fieldmap
    (e.g. ADNI Database) by elastically register DWIs to their respective baseline T1-weighted
    structural scans using an inverse consistent registration algorithm with a mutual information cost
    function (SyN algorithm). This workflow allows also a coregistration of DWIs with their respective
    baseline T1-weighted structural scans in order to latter combine tracks and cortex parcellation.

    Parameters
    ----------
    base_dir: str
        Working directory, which contains all of the intermediary data generated.

    delete_cache: bool
        If True, part of the temporary data is automatically deleted after usage.

    name: str, optional
        Name of the pipeline.

    Warnings
    --------
        This workflow rotates the b-vectors.

    Notes
    -----
        Nir et al. (2015): Connectivity network measures predict volumetric atrophy in mild cognitive impairment
        Leow et al. (2007): Statistical Properties of Jacobian Maps and the Realization of
        Unbiased Large Deformation Nonlinear Image Registration
    """
    import nipype.interfaces.ants as ants
    import nipype.interfaces.c3 as c3
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from .dwi_preprocessing_using_t1_utils import (
        ants_apply_transforms,
        broadcast_matrix_filename_to_match_b_vector_length,
        change_itk_transform_type,
        delete_temp_dirs,
        rotate_b_vectors,
    )

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["T1", "DWI", "bvec"]), name="inputnode"
    )

    split = pe.Node(fsl.Split(dimension="t"), name="SplitDWIs")
    pick_ref = pe.Node(niu.Select(), name="Pick_b0")
    pick_ref.inputs.index = [0]

    flirt_b0_2_t1 = pe.Node(interface=fsl.FLIRT(dof=6), name="flirt_B0_2_T1")
    flirt_b0_2_t1.inputs.interp = "spline"
    flirt_b0_2_t1.inputs.cost = "normmi"
    flirt_b0_2_t1.inputs.cost_func = "normmi"

    expend_matrix = pe.Node(
        interface=niu.Function(
            input_names=["matrix_filename", "b_vectors_filename"],
            output_names=["out_matrix_list"],
            function=broadcast_matrix_filename_to_match_b_vector_length,
        ),
        name="expend_matrix",
    )

    rot_bvec = pe.Node(
        niu.Function(
            input_names=["matrix_filenames", "b_vectors_filename"],
            output_names=["out_file"],
            function=rotate_b_vectors,
        ),
        name="Rotate_Bvec",
    )

    ants_registration = pe.Node(
        interface=ants.registration.RegistrationSynQuick(
            transform_type="br", dimension=3
        ),
        name="antsRegistrationSyNQuick",
    )

    c3d_flirt2ants = pe.Node(c3.C3dAffineTool(), name="fsl_reg_2_itk")
    c3d_flirt2ants.inputs.itk_transform = True
    c3d_flirt2ants.inputs.fsl2ras = True

    change_transform = pe.Node(
        niu.Function(
            input_names=["input_affine_file"],
            output_names=["updated_affine_file"],
            function=change_itk_transform_type,
        ),
        name="change_transform_type",
    )

    merge_transform = pe.Node(niu.Merge(3), name="MergeTransforms")

    # the nipype function is not used since the results it gives are not
    # as good as the ones we get by using the command directly.
    apply_transform_image = pe.MapNode(
        interface=niu.Function(
            input_names=[
                "fixed_image",
                "moving_image",
                "transforms",
                "warped_image",
                "output_warped_image",
            ],
            output_names=["warped_image"],
            function=ants_apply_transforms,
        ),
        iterfield=["moving_image"],
        name="warp_image",
    )
    apply_transform_image.inputs.warped_image = "out_warped.nii.gz"
    apply_transform_image.inputs.output_warped_image = True

    apply_transform_field = pe.MapNode(
        interface=niu.Function(
            input_names=[
                "fixed_image",
                "moving_image",
                "transforms",
                "warped_image",
                "output_warped_image",
            ],
            output_names=["warped_image"],
            function=ants_apply_transforms,
        ),
        iterfield=["moving_image"],
        name="warp_field",
    )
    apply_transform_field.inputs.warped_image = "out_warped_field.nii.gz"
    apply_transform_field.inputs.output_warped_image = False

    jacobian = pe.MapNode(
        interface=ants.CreateJacobianDeterminantImage(),
        iterfield=["deformationField"],
        name="jacobian",
    )

    jacobian.inputs.imageDimension = 3
    jacobian.inputs.outputImage = "Jacobian_image.nii.gz"

    jacmult = pe.MapNode(
        fsl.MultiImageMaths(op_string="-mul %s"),
        iterfield=["in_file", "operand_files"],
        name="ModulateDWIs",
    )

    thres = pe.MapNode(
        fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="RemoveNegative"
    )

    merge = pe.Node(fsl.Merge(dimension="t"), name="MergeDWIs")
    # Delete the temporary directory that takes too much place

    delete_warp_field_tmp = pe.Node(
        name="deletewarpfieldtmp",
        interface=niu.Function(
            inputs=["checkpoint", "dir_to_del", "base_dir"],
            function=delete_temp_dirs,
        ),
    )
    delete_warp_field_tmp.inputs.base_dir = base_dir
    delete_warp_field_tmp.inputs.dir_to_del = [
        apply_transform_field.name,
        jacobian.name,
        jacmult.name,
        thres.name,
        apply_transform_image.name,
    ]

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "DWI_2_T1_Coregistration_matrix",
                "epi_correction_deformation_field",
                "epi_correction_affine_transform",
                "epi_correction_image_warped",
                "DWIs_epicorrected",
                "warp_epi",
                "out_bvec",
            ]
        ),
        name="outputnode",
    )

    wf = pe.Workflow(name="epi_pipeline")
    # fmt: off
    wf.connect(
        [
            (inputnode, split, [("DWI", "in_file")]),
            (split, pick_ref, [("out_files", "inlist")]),
            (pick_ref, flirt_b0_2_t1, [("out", "in_file")]),
            (inputnode, flirt_b0_2_t1, [("T1", "reference")]),
            (inputnode, rot_bvec, [("bvec", "b_vectors_filename")]),
            (flirt_b0_2_t1, expend_matrix, [("out_matrix_file", "matrix_filename")]),
            (inputnode, expend_matrix, [("bvec", "b_vectors_filename")]),
            (expend_matrix, rot_bvec, [("out_matrix_list", "matrix_filenames")]),
            (inputnode, ants_registration, [("T1", "fixed_image")]),
            (flirt_b0_2_t1, ants_registration, [("out_file", "moving_image")]),
            (inputnode, c3d_flirt2ants, [("T1", "reference_file")]),
            (pick_ref, c3d_flirt2ants, [("out", "source_file")]),
            (flirt_b0_2_t1, c3d_flirt2ants, [("out_matrix_file", "transform_file")]),
            (c3d_flirt2ants, change_transform, [("itk_transform", "input_affine_file")]),
            (change_transform, merge_transform, [("updated_affine_file", "in1")]),
            (ants_registration, merge_transform, [("out_matrix", "in2")]),
            (ants_registration, merge_transform, [("forward_warp_field", "in3")]),
            (inputnode, apply_transform_image, [("T1", "fixed_image")]),
            (split, apply_transform_image, [("out_files", "moving_image")]),
            (merge_transform, apply_transform_image, [("out", "transforms")]),

            (inputnode, apply_transform_field, [("T1", "fixed_image")]),
            (split, apply_transform_field, [("out_files", "moving_image")]),
            (merge_transform, apply_transform_field, [("out", "transforms")]),

            (apply_transform_field, jacobian, [("warped_image", "deformationField")]),
            (apply_transform_image, jacmult, [("warped_image", "operand_files")]),
            (jacobian, jacmult, [("jacobian_image", "in_file")]),
            (jacmult, thres, [("out_file", "in_file")]),
            (thres, merge, [("out_file", "in_files")]),

            (merge, outputnode, [("merged_file", "DWIs_epicorrected")]),
            (flirt_b0_2_t1, outputnode, [("out_matrix_file", "DWI_2_T1_Coregistration_matrix")]),
            (ants_registration, outputnode, [("forward_warp_field", "epi_correction_deformation_field"),
                                             ("out_matrix", "epi_correction_affine_transform"),
                                             ("warped_image", "epi_correction_image_warped")]),
            (merge_transform, outputnode, [("out", "warp_epi")]),
            (rot_bvec, outputnode, [("out_file", "out_bvec")]),



        ]
    )
    if delete_cache:
        wf.connect(
            [
                (merge, delete_warp_field_tmp, [("merged_file", "checkpoint")])
            ]
        )
    # fmt: on
    return wf


def b0_flirt_pipeline(num_b0s, name="b0_coregistration"):
    """Rigid registration of the B0 dataset onto the first volume.

    Rigid registration is achieved using FLIRT and the normalized correlation.

    Parameters
    ----------
    num_b0s : int
        Number of B0 volumes in the dataset.

    name : str, optional
        Name of the workflow.
        Default="b0_coregistration".

    Inputnode:
        in_file(str): B0 dataset.

    Outputnode
        out_b0_reg(str): The set of B0 volumes registered to the first volume.
    """
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces import fsl

    from clinica.utils.image import merge_volumes_time_dimension_task

    inputnode = pe.Node(niu.IdentityInterface(fields=["in_file"]), name="inputnode")
    fslroi_ref = pe.Node(fsl.ExtractROI(args="0 1"), name="b0_reference")
    tsize = num_b0s - 1
    fslroi_moving = pe.Node(fsl.ExtractROI(args="1 " + str(tsize)), name="b0_moving")
    split_moving = pe.Node(fsl.Split(dimension="t"), name="split_b0_moving")

    bet_ref = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True), name="bet_ref")

    dilate = pe.Node(
        fsl.maths.MathsCommand(nan2zeros=True, args="-kernel sphere 5 -dilM"),
        name="mask_dilate",
    )

    flirt = pe.MapNode(
        fsl.FLIRT(
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
        ),
        name="b0_co_registration",
        iterfield=["in_file"],
    )

    merge = pe.Node(fsl.Merge(dimension="t"), name="merge_registered_b0s")
    thres = pe.MapNode(
        fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="remove_negative"
    )
    insert_ref = pe.Node(
        niu.Function(
            input_names=["volume1", "volume2"],
            output_names=["out_file"],
            function=merge_volumes_time_dimension_task,
        ),
        name="concat_ref_moving",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_file", "out_xfms"]), name="outputnode"
    )

    wf = pe.Workflow(name=name)
    # fmt: off
    wf.connect(
        [
            (inputnode, fslroi_ref, [("in_file", "in_file")]),
            (inputnode, fslroi_moving, [("in_file", "in_file")]),
            (fslroi_moving, split_moving, [("roi_file", "in_file")]),
            (fslroi_ref, bet_ref, [("roi_file", "in_file")]),
            (bet_ref, dilate, [("mask_file", "in_file")]),
            (dilate, flirt, [("out_file", "ref_weight"),
                             ("out_file", "in_weight")]),
            (fslroi_ref, flirt, [("roi_file", "reference")]),
            (split_moving, flirt, [("out_files", "in_file")]),
            (flirt, thres, [("out_file", "in_file")]),
            (thres, merge, [("out_file", "in_files")]),
            (merge, insert_ref, [("merged_file", "volume2")]),
            (fslroi_ref, insert_ref, [("roi_file", "volume1")]),
            (insert_ref, outputnode, [("out_file", "out_file")]),
            (flirt, outputnode, [("out_matrix_file", "out_xfms")])
        ]
    )
    # fmt: off

    return wf
