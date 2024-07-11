def get_wf(
    subject_id: str,
    session_id: str,
    pvc_psf_tsv: str,
    caps_dir: str,
    pet: str,
    orig_nu: str,
    white_surface_left: str,
    white_surface_right: str,
    working_directory_subjects: str,
    acq_label: str,
    csv_segmentation: str,
    suvr_reference_region: str,
    matscript_folder_inverse_deformation: str,
    destrieux_left: str,
    destrieux_right: str,
    desikan_left: str,
    desikan_right: str,
    is_longitudinal: bool,
    output_dir=None,
):
    """get_wf create a full workflow for only one subject, and then executes it.

    Parameters
    ----------
    subject_id : str
        The subject ID.

    session_id : str
        The session ID.

    pvc_psf_tsv : str
        The path the TSV file containing information on the point spread function (PSF).

    caps_dir : str
        The path to the CAPS directory.

    pet : str
        The path to the PET image in the BIDS directory.

    orig_nu : str
        The path to the 'orig_nu' file (must be in the CAPS directory, in 'mri').

    white_surface_left : str
        The path to the left white surface in native space of subject.

    white_surface_right : str
        The path to the right white surface in native space of subject.

    working_directory_subjects : str
        ???

    acq_label : str
        The PET tracer to consider.

    csv_segmentation : str
        The path to the CSV for the segmentation (problems encountered while using __file__).

    suvr_reference_region : str
        The label of the SUVR reference region.

    matscript_folder_inverse_deformation : str
        The path to the current folder containing the matlab script
        used to call the spm function for the inverse deformation.

    destrieux_left : str
        The path to the destrieux parcellation for the left hemisphere.

    destrieux_right : str
        The path to the destrieux parcellation for the right hemisphere.

    desikan_left : str
        The path to the desikan parcellation for the left hemisphere.

    desikan_right : str
        The path to the desikan parcellation for the right hemisphere.

    is_longitudinal : bool
        Whether the pipeline is PETSurface or PETSurfaceLongitudinal.

    output_dir : str, optional
        The path to the output directory.
    """
    import os
    from pathlib import Path

    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces.freesurfer import ApplyVolTransform, MRIConvert, Tkregister2
    from nipype.interfaces.fsl import Merge
    from nipype.interfaces.petpvc import PETPVC
    from nipype.interfaces.spm import Coregister, Normalize12

    from clinica.pipelines.pet.surface.tasks import (
        compute_average_pet_signal_based_on_annotations_task,
        compute_weighted_mean_surface_task,
        get_mid_surface_task,
        make_label_conversion_task,
        normalize_suvr_task,
        perform_gtmseg_task,
        project_onto_fsaverage_task,
        reformat_surfname_task,
        remove_nan_from_image_task,
        run_apply_inverse_deformation_field_spm_standalone_task,
        run_apply_inverse_deformation_field_task,
        run_mri_surf2surf_task,
        run_mri_vol2surf_task,
        run_mris_expand_task,
    )
    from clinica.pipelines.pet.surface.utils import (
        get_output_dir,
        get_regexp_substitutions,
    )
    from clinica.pipelines.pet.utils import get_suvr_mask, read_psf_information
    from clinica.pipelines.utils import FreeSurferAnnotationImage
    from clinica.utils.filemanip import get_subject_id, load_volume, unzip_nii
    from clinica.utils.pet import SUVRReferenceRegion, Tracer
    from clinica.utils.spm import get_tpm, use_spm_standalone_if_available
    from clinica.utils.ux import print_begin_image

    spm_standalone_is_available = use_spm_standalone_if_available()
    acq_label = Tracer(acq_label)
    suvr_reference_region = SUVRReferenceRegion(suvr_reference_region)
    image_id = get_subject_id(pet)
    try:
        load_volume(pet)
    except ValueError as e:
        raise ValueError(
            f"Clinica could not load volumes for {image_id.replace('_', ' | ')}. {str(e)}"
        )
    print_begin_image(image_id)

    # Creation of workflow
    # 1 Creation of node
    unzip_pet = pe.Node(
        niu.Function(
            input_names=["in_file"], output_names=["out_file"], function=unzip_nii
        ),
        name="unzip_pet",
    )

    unzip_orig_nu = unzip_pet.clone(name="unzip_orig_nu")

    unzip_mask = unzip_pet.clone(name="unzip_mask")
    unzip_mask.inputs.in_file = str(get_suvr_mask(suvr_reference_region))

    coreg = pe.Node(Coregister(), name="coreg")

    convert_mgh = pe.Node(MRIConvert(), name="convert_mgh")

    removenan = pe.Node(
        niu.Function(
            input_names=["image_path"],
            output_names=["vol_wo_nan"],
            function=remove_nan_from_image_task,
        ),
        name="removenan",
    )

    gtmsegmentation = pe.Node(
        niu.Function(
            input_names=["caps_dir", "subject_id", "session_id", "is_longitudinal"],
            output_names=["gtmseg_file"],
            function=perform_gtmseg_task,
        ),
        name="gtmseg",
    )
    gtmsegmentation.inputs.caps_dir = caps_dir
    gtmsegmentation.inputs.subject_id = subject_id
    gtmsegmentation.inputs.session_id = session_id
    gtmsegmentation.inputs.is_longitudinal = is_longitudinal

    tkregister = pe.Node(Tkregister2(reg_header=True), name="tkreg")

    convert_gtmseg = convert_mgh.clone(name="convert_gtmseg")

    labelconversion = pe.Node(
        niu.Function(
            input_names=["gtmsegfile", "csv"],
            output_names=["list_of_regions"],
            function=make_label_conversion_task,
        ),
        name="conversion_of_labels",
    )
    labelconversion.inputs.csv = csv_segmentation

    if not os.path.exists(labelconversion.inputs.csv):
        raise FileNotFoundError(
            f"CSV file : {labelconversion.inputs.csv} does not exist."
        )

    merge_volume = pe.Node(
        Merge(output_type="NIFTI_GZ", dimension="t"), name="merge_volume"
    )
    vol2vol = pe.Node(
        ApplyVolTransform(reg_header=True, interp="trilin"), name="vol2vol"
    )
    vol2vol_mask = pe.Node(
        ApplyVolTransform(reg_header=True, interp="nearest"), name="vol2vol_mask"
    )
    normalize12 = pe.Node(
        Normalize12(
            tpm=get_tpm(),
            affine_regularization_type="mni",
            jobtype="est",
            bias_fwhm=60,
            bias_regularization=0.0001,
            warping_regularization=[0, 0.001, 0.5, 0.05, 0.2],
        ),
        name="normalize_to_MNI",
    )
    apply_inverse_deformation = pe.Node(
        niu.Function(
            input_names=["target", "deformation_field", "img", "matscript_folder"],
            output_names=["freesurfer_space_eroded_mask"],
            function=(
                run_apply_inverse_deformation_field_spm_standalone_task
                if spm_standalone_is_available
                else run_apply_inverse_deformation_field_task
            ),
        ),
        name="applyInverseDeformation",
    )
    apply_inverse_deformation.inputs.matscript_folder = (
        matscript_folder_inverse_deformation
    )
    pons_normalization = pe.Node(
        niu.Function(
            input_names=["pet_image", "mask"],
            output_names=["suvr_image"],
            function=normalize_suvr_task,
        ),
        name="pons_normalization",
    )
    # read_psf_information expects a list of subjects/sessions and returns a list of PSF
    psf_info = read_psf_information(
        Path(pvc_psf_tsv),
        [subject_id],
        [session_id],
        acq_label,
    )[0]
    pvc = pe.Node(PETPVC(pvc="IY"), name="petpvc")
    pvc.inputs.fwhm_x = psf_info[0]
    pvc.inputs.fwhm_y = psf_info[1]
    pvc.inputs.fwhm_z = psf_info[2]

    reformat_surface_name = pe.Node(
        niu.Function(
            input_names=["hemisphere", "left_surface", "right_surface"],
            output_names=["out"],
            function=reformat_surfname_task,
        ),
        name="reformat_surface_name",
    )
    reformat_surface_name.inputs.left_surface = white_surface_left
    reformat_surface_name.inputs.right_surface = white_surface_right
    reformat_surface_name.iterables = ("hemisphere", ["lh", "rh"])

    mris_exp = pe.Node(
        niu.Function(
            input_names=["surface", "output_dir"],
            output_names=["out_surface"],
            function=run_mris_expand_task,
        ),
        name="mris_expand_white",
    )
    mris_exp.inputs.output_dir = output_dir

    surf_conversion = pe.MapNode(
        niu.Function(
            input_names=[
                "surface",
                "registration",
                "gtmseg_file",
                "subject_id",
                "session_id",
                "caps_dir",
                "is_longitudinal",
                "output_dir",
            ],
            output_names=["tval"],
            function=run_mri_surf2surf_task,
        ),
        name="surf_conversion",
        iterfield=["in_surface"],
    )
    surf_conversion.inputs.subject_id = subject_id
    surf_conversion.inputs.session_id = session_id
    surf_conversion.inputs.caps_dir = caps_dir
    surf_conversion.inputs.is_longitudinal = is_longitudinal
    surf_conversion.inputs.output_dir = output_dir

    vol_on_surf = pe.MapNode(
        niu.Function(
            input_names=[
                "pet_volume",
                "surface",
                "subject_id",
                "session_id",
                "caps_dir",
                "gtmsegfile",
                "is_longitudinal",
                "output_dir",
            ],
            output_names=["output"],
            function=run_mri_vol2surf_task,
        ),
        name="vol_on_surf",
        iterfield=["surface"],
    )
    vol_on_surf.inputs.subject_id = subject_id
    vol_on_surf.inputs.session_id = session_id
    vol_on_surf.inputs.caps_dir = caps_dir
    vol_on_surf.inputs.is_longitudinal = is_longitudinal
    vol_on_surf.inputs.output_dir = output_dir

    normal_average = pe.Node(
        niu.Function(
            input_names=["surfaces", "output_dir"],
            output_names=["out_surface"],
            function=compute_weighted_mean_surface_task,
        ),
        name="normal_average",
    )
    normal_average.inputs.output_dir = output_dir

    project_on_fsaverage = pe.Node(
        niu.Function(
            input_names=[
                "projection",
                "subject_id",
                "caps_dir",
                "session_id",
                "fwhm",
                "is_longitudinal",
                "output_dir",
            ],
            output_names=["out_fsaverage"],
            function=project_onto_fsaverage_task,
        ),
        name="project_on_fsaverage",
    )
    project_on_fsaverage.iterables = ("fwhm", [0, 5, 10, 15, 20, 25])
    project_on_fsaverage.inputs.subject_id = subject_id
    project_on_fsaverage.inputs.session_id = session_id
    project_on_fsaverage.inputs.caps_dir = caps_dir
    project_on_fsaverage.inputs.is_longitudinal = is_longitudinal
    project_on_fsaverage.inputs.output_dir = output_dir

    extract_mid_surface = pe.Node(
        niu.Function(
            input_names=["surfaces"],
            output_names=["mid_surface"],
            function=get_mid_surface_task,
        ),
        name="extract_mid_surface",
    )

    gather_pet_projection = pe.JoinNode(
        niu.IdentityInterface(fields=["pet_projection_lh_rh"]),
        name="gather_pet_projection_hemisphere",
        joinsource="reformat_surface_name",
        joinfield=["pet_projection_lh_rh"],
    )

    atlas_tsv = pe.Node(
        niu.Function(
            input_names=["pet_projections", "atlas_files"],
            output_names=["destrieux_tsv", "desikan_tsv"],
            function=compute_average_pet_signal_based_on_annotations_task,
        ),
        name="atlas_tsv",
    )
    atlas_tsv.inputs.atlas_files = {
        "destrieux": FreeSurferAnnotationImage.from_raw(
            destrieux_left, destrieux_right
        ),
        "desikan": FreeSurferAnnotationImage.from_raw(desikan_left, desikan_right),
    }

    # 2 creation of workflow : working dir, inputnode, outputnode and datasink
    name_workflow = subject_id.replace("-", "_") + "_" + session_id.replace("-", "_")
    if is_longitudinal:
        name_workflow += "_long"

    wf = pe.Workflow(name=name_workflow)
    wf.base_dir = working_directory_subjects

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "orig_nu",
                "pet",
                "psf",
                "white_surface_left",
                "white_surface_right",
            ]
        ),
        name="inputnode",
        mandatory_inputs=True,
    )

    inputnode.inputs.orig_nu = orig_nu
    inputnode.inputs.pet = pet
    inputnode.inputs.psf = pvc_psf_tsv
    inputnode.inputs.white_surface_right = white_surface_right
    inputnode.inputs.white_surface_left = white_surface_left

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "mid_surf",
                "projection_native_subject",
                "projection_fsaverage_smoothed",
                "destrieux_tsv",
                "desikan_tsv",
            ]
        ),
        name="outputnode",
        mandatory_inputs=True,
    )

    datasink = pe.Node(nio.DataSink(), name="sinker")
    datasink.inputs.base_directory = str(
        get_output_dir(is_longitudinal, caps_dir, subject_id, session_id)
    )
    datasink.inputs.parameterization = True
    datasink.inputs.regexp_substitutions = get_regexp_substitutions(
        acq_label, suvr_reference_region, is_longitudinal
    )

    # 3 Connecting the nodes

    # TODO(@arnaud.marcoux): Add titles for sections of connections.
    #   Could be useful to add sections title to group similar connections
    #   together.

    connections = [
        (inputnode, unzip_pet, [("pet", "in_file")]),
        (unzip_pet, coreg, [("out_file", "source")]),
        (inputnode, convert_mgh, [("orig_nu", "in_file")]),
        (convert_mgh, unzip_orig_nu, [("out_file", "in_file")]),
        (unzip_orig_nu, coreg, [("out_file", "target")]),
        (coreg, removenan, [("coregistered_source", "image_path")]),
        (removenan, vol2vol, [("vol_wo_nan", "source_file")]),
        (inputnode, tkregister, [("orig_nu", "target_image")]),
        (unzip_orig_nu, normalize12, [("out_file", "image_to_align")]),
        (unzip_mask, apply_inverse_deformation, [("out_file", "img")]),
        (
            normalize12,
            apply_inverse_deformation,
            [("deformation_field", "deformation_field")],
        ),
        (unzip_orig_nu, apply_inverse_deformation, [("out_file", "target")]),
        (
            apply_inverse_deformation,
            vol2vol_mask,
            [("freesurfer_space_eroded_mask", "source_file")],
        ),
        (gtmsegmentation, vol2vol_mask, [("gtmseg_file", "target_file")]),
        (gtmsegmentation, tkregister, [("gtmseg_file", "moving_image")]),
        (gtmsegmentation, convert_gtmseg, [("gtmseg_file", "in_file")]),
        (gtmsegmentation, vol2vol, [("gtmseg_file", "target_file")]),
        (vol2vol, pons_normalization, [("transformed_file", "pet_image")]),
        (vol2vol_mask, pons_normalization, [("transformed_file", "mask")]),
        (convert_gtmseg, labelconversion, [("out_file", "gtmsegfile")]),
        (labelconversion, merge_volume, [("list_of_regions", "in_files")]),
        (merge_volume, pvc, [("merged_file", "mask_file")]),
        (pons_normalization, pvc, [("suvr_image", "in_file")]),
        (reformat_surface_name, mris_exp, [("out", "surface")]),
        (mris_exp, extract_mid_surface, [("out_surface", "surfaces")]),
        (mris_exp, surf_conversion, [("out_surface", "surface")]),
        (tkregister, surf_conversion, [("reg_file", "registration")]),
        (gtmsegmentation, surf_conversion, [("gtmseg_file", "gtmseg_file")]),
        (pvc, vol_on_surf, [("out_file", "pet_volume")]),
        (surf_conversion, vol_on_surf, [("tval", "surface")]),
        (gtmsegmentation, vol_on_surf, [("gtmseg_file", "gtmsegfile")]),
        (vol_on_surf, normal_average, [("output", "surfaces")]),
        (normal_average, project_on_fsaverage, [("out_surface", "projection")]),
        (
            normal_average,
            gather_pet_projection,
            [("out_surface", "pet_projection_lh_rh")],
        ),
        (
            gather_pet_projection,
            atlas_tsv,
            [("pet_projection_lh_rh", "pet_projections")],
        ),
        (atlas_tsv, outputnode, [("destrieux_tsv", "destrieux_tsv")]),
        (atlas_tsv, outputnode, [("desikan_tsv", "desikan_tsv")]),
        (
            project_on_fsaverage,
            outputnode,
            [("out_fsaverage", "projection_fsaverage_smoothed")],
        ),
        (extract_mid_surface, outputnode, [("mid_surface", "mid_surf")]),
        (normal_average, outputnode, [("out_surface", "projection_native_subject")]),
        (
            outputnode,
            datasink,
            [("projection_fsaverage_smoothed", "projection_fsaverage")],
        ),
        (outputnode, datasink, [("mid_surf", "midsurface")]),
        (outputnode, datasink, [("projection_native_subject", "projection_native")]),
        (outputnode, datasink, [("destrieux_tsv", "destrieux_tsv")]),
        (outputnode, datasink, [("desikan_tsv", "desikan_tsv")]),
    ]
    wf.connect(connections)
    # wf.write_graph(graph2use='flat')
    wf.run()
