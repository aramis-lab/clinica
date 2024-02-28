def get_wf(
    subject_id,
    session_id,
    pvc_psf_tsv,
    caps_dir,
    pet,
    orig_nu,
    white_surface_left,
    white_surface_right,
    working_directory_subjects,
    acq_label,
    csv_segmentation,
    suvr_reference_region,
    matscript_folder_inverse_deformation,
    destrieux_left,
    destrieux_right,
    desikan_left,
    desikan_right,
    spm_standalone_is_available,
    is_longitudinal: bool,
):
    """get_wf create a full workflow for only one subject, and then executes it

    Args:
        subject_id (string): The subject ID
        session_id (string): The session ID
        pvc_psf_tsv (string): Path the TSV file containing information on the point spread function (PSF)
        caps_dir (string): Path to the CAPS directory
        pet (string): Path to the PET image in the bids directory
        orig_nu (string): Path to the orig_nu file (must be in the CAPS directory, in mri)
        white_surface_left (string): Path to the left white surface in native space of subject
        white_surface_right (string): Path to the right white surface in native space of subject
        working_directory_subjects (string):
        acq_label (string):
        csv_segmentation (string): Path to the CSV for the segmentation (problems encountered while using __file__)
        suvr_reference_region (string): Label of the SUVR reference region
        matscript_folder_inverse_deformation (string): Path to the current folder containing the matlab script used to call the spm function for the inverse deformation
        destrieux_left (string):
        destrieux_right (string):
        desikan_left (string):
        desikan_right (string):
        spm_standalone_is_available (string):
        is_longitudinal : bool

    Returns:
        Void
    """
    import os

    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces.freesurfer import ApplyVolTransform, MRIConvert, Tkregister2
    from nipype.interfaces.fsl import Merge
    from nipype.interfaces.petpvc import PETPVC
    from nipype.interfaces.spm import Coregister, Normalize12

    import clinica.pipelines.pet_surface.pet_surface_utils as utils
    from clinica.pipelines.pet.utils import get_suvr_mask, read_psf_information
    from clinica.utils.filemanip import get_subject_id, load_volume, unzip_nii
    from clinica.utils.spm import get_tpm
    from clinica.utils.ux import print_begin_image

    from .tasks import perform_gtmseg_task, remove_nan_from_image_task
    from .utils import get_output_dir

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
            function=utils.make_label_conversion,
        ),
        name="conversion_of_labels",
    )

    labelconversion.inputs.csv = csv_segmentation

    if not os.path.exists(labelconversion.inputs.csv):
        raise Exception("CSV file : " + labelconversion.inputs.csv + " does not exist.")

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

    if spm_standalone_is_available:
        fun_apply_inverse_deformation = (
            utils.runApplyInverseDeformationField_SPM_standalone
        )
    else:
        fun_apply_inverse_deformation = utils.runApplyInverseDeformationField
    apply_inverse_deformation = pe.Node(
        niu.Function(
            input_names=["target", "deformation_field", "img", "matscript_folder"],
            output_names=["freesurfer_space_eroded_mask"],
            function=fun_apply_inverse_deformation,
        ),
        name="applyInverseDeformation",
    )
    apply_inverse_deformation.inputs.matscript_folder = (
        matscript_folder_inverse_deformation
    )

    pons_normalization = pe.Node(
        niu.Function(
            input_names=["pet_path", "mask"],
            output_names=["suvr"],
            function=utils.suvr_normalization,
        ),
        name="pons_normalization",
    )
    pons_normalization.inputs.pet_tracer = acq_label

    # read_psf_information expects a list of subjects/sessions and returns a list of PSF
    psf_info = read_psf_information(pvc_psf_tsv, [subject_id], [session_id], acq_label)[
        0
    ]

    pvc = pe.Node(PETPVC(pvc="IY"), name="petpvc")
    pvc.inputs.fwhm_x = psf_info[0]
    pvc.inputs.fwhm_y = psf_info[1]
    pvc.inputs.fwhm_z = psf_info[2]

    reformat_surface_name = pe.Node(
        niu.Function(
            input_names=["hemi", "left_surface", "right_surface"],
            output_names=["out"],
            function=utils.reformat_surfname,
        ),
        name="reformat_surface_name",
    )
    reformat_surface_name.inputs.left_surface = white_surface_left
    reformat_surface_name.inputs.right_surface = white_surface_right
    reformat_surface_name.iterables = ("hemi", ["lh", "rh"])

    mris_exp = pe.Node(
        niu.Function(
            input_names=["in_surface"],
            output_names=["out_surface"],
            function=utils.mris_expand,
        ),
        name="mris_expand_white",
    )

    surf_conversion = pe.MapNode(
        niu.Function(
            input_names=[
                "in_surface",
                "reg_file",
                "gtmsegfile",
                "subject_id",
                "session_id",
                "caps_dir",
                "is_longitudinal",
            ],
            output_names=["tval"],
            function=utils.surf2surf,
        ),
        name="surf_conversion",
        iterfield=["in_surface"],
    )
    surf_conversion.inputs.subject_id = subject_id
    surf_conversion.inputs.session_id = session_id
    surf_conversion.inputs.caps_dir = caps_dir
    surf_conversion.inputs.is_longitudinal = is_longitudinal

    vol_on_surf = pe.MapNode(
        niu.Function(
            input_names=[
                "volume",
                "surface",
                "subject_id",
                "session_id",
                "caps_dir",
                "gtmsegfile",
                "is_longitudinal",
            ],
            output_names=["output"],
            function=utils.vol2surf,
        ),
        name="vol_on_surf",
        iterfield=["surface"],
    )
    vol_on_surf.inputs.subject_id = subject_id
    vol_on_surf.inputs.session_id = session_id
    vol_on_surf.inputs.caps_dir = caps_dir
    vol_on_surf.inputs.is_longitudinal = is_longitudinal

    normal_average = pe.Node(
        niu.Function(
            input_names=["in_surfaces"],
            output_names=["out_surface"],
            function=utils.weighted_mean,
        ),
        name="normal_average",
    )

    project_on_fsaverage = pe.Node(
        niu.Function(
            input_names=[
                "projection",
                "subject_id",
                "caps_dir",
                "session_id",
                "fwhm",
                "is_longitudinal",
            ],
            output_names=["out_fsaverage"],
            function=utils.fsaverage_projection,
        ),
        name="project_on_fsaverage",
    )
    project_on_fsaverage.iterables = ("fwhm", [0, 5, 10, 15, 20, 25])
    project_on_fsaverage.inputs.subject_id = subject_id
    project_on_fsaverage.inputs.session_id = session_id
    project_on_fsaverage.inputs.caps_dir = caps_dir
    project_on_fsaverage.inputs.is_longitudinal = is_longitudinal

    extract_mid_surface = pe.Node(
        niu.Function(
            input_names=["in_surfaces"],
            output_names=["mid_surface"],
            function=utils.get_mid_surface,
        ),
        name="extract_mid_surface",
    )

    surface_atlas = {
        "destrieux": {"lh": destrieux_left, "rh": destrieux_right},
        "desikan": {"lh": desikan_left, "rh": desikan_right},
    }

    gather_pet_projection = pe.JoinNode(
        niu.IdentityInterface(fields=["pet_projection_lh_rh"]),
        name="gather_pet_projection_hemisphere",
        joinsource="reformat_surface_name",
        joinfield=["pet_projection_lh_rh"],
    )

    atlas_tsv = pe.Node(
        niu.Function(
            input_names=["pet", "atlas_files"],
            output_names=["destrieux_tsv", "desikan_tsv"],
            function=utils.produce_tsv,
        ),
        name="atlas_tsv",
    )
    atlas_tsv.inputs.atlas_files = surface_atlas

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
    cross_sectional_regexp_substitutions = [
        # Mid surface
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
            r"\1/\2_\3_hemi-\4_midcorticalsurface",
        ),
        # Projection in native space
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_native\/.*_hemi_([a-z]+).*",
            r"\1/\2_\3_trc-"
            + acq_label
            + r"_pet_space-native_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\4_projection.mgh",
        ),
        # Projection in fsaverage
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
            r"\1/\2_\3_trc-"
            + acq_label
            + r"_pet_space-fsaverage_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\4_fwhm-\5_projection.mgh",
        ),
        # TSV file for Destrieux atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/destrieux_tsv\/destrieux.tsv",
            r"\1/atlas_statistics/\2_\3_trc-"
            + acq_label
            + "_pet_space-destrieux_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
        # TSV file for Desikan atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/desikan_tsv\/desikan.tsv",
            r"\1/atlas_statistics/\2_\3_trc-"
            + acq_label
            + "_pet_space-desikan_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
    ]
    longitudinal_regexp_substitutions = [
        # Mid surface
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
            r"\1/\2_\3_\4_hemi-\5_midcorticalsurface",
        ),
        # Projection in native space
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/projection_native\/.*_hemi_([a-z]+).*",
            r"\1/\2_\3_\4_trc-"
            + acq_label
            + r"_pet_space-native_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\5_projection.mgh",
        ),
        # Projection in fsaverage
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
            r"\1/\2_\3_\4_trc-"
            + acq_label
            + r"_pet_space-fsaverage_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\5_fwhm-\6_projection.mgh",
        ),
        # TSV file for Destrieux atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/destrieux_tsv\/destrieux.tsv",
            r"\1/atlas_statistics/\2_\3_\4_trc-"
            + acq_label
            + "_pet_space-destrieux_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
        # TSV file for Desikan atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/desikan_tsv\/desikan.tsv",
            r"\1/atlas_statistics/\2_\3_\4_trc-"
            + acq_label
            + "_pet_space-desikan_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
    ]
    if is_longitudinal:
        datasink.inputs.regexp_substitutions = longitudinal_regexp_substitutions
    else:
        datasink.inputs.regexp_substitutions = cross_sectional_regexp_substitutions

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
        (vol2vol, pons_normalization, [("transformed_file", "pet_path")]),
        (vol2vol_mask, pons_normalization, [("transformed_file", "mask")]),
        (convert_gtmseg, labelconversion, [("out_file", "gtmsegfile")]),
        (labelconversion, merge_volume, [("list_of_regions", "in_files")]),
        (merge_volume, pvc, [("merged_file", "mask_file")]),
        (pons_normalization, pvc, [("suvr", "in_file")]),
        (reformat_surface_name, mris_exp, [("out", "in_surface")]),
        (mris_exp, extract_mid_surface, [("out_surface", "in_surfaces")]),
        (mris_exp, surf_conversion, [("out_surface", "in_surface")]),
        (tkregister, surf_conversion, [("reg_file", "reg_file")]),
        (gtmsegmentation, surf_conversion, [("gtmseg_file", "gtmsegfile")]),
        (pvc, vol_on_surf, [("out_file", "volume")]),
        (surf_conversion, vol_on_surf, [("tval", "surface")]),
        (gtmsegmentation, vol_on_surf, [("gtmseg_file", "gtmsegfile")]),
        (vol_on_surf, normal_average, [("output", "in_surfaces")]),
        (normal_average, project_on_fsaverage, [("out_surface", "projection")]),
        (
            normal_average,
            gather_pet_projection,
            [("out_surface", "pet_projection_lh_rh")],
        ),
        (gather_pet_projection, atlas_tsv, [("pet_projection_lh_rh", "pet")]),
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
