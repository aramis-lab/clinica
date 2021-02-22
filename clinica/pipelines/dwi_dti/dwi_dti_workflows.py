# coding: utf8


def register_dti_maps_on_atlas(
    working_directory=None, name="register_dti_maps_on_atlas"
):
    """Register FA-map on a subject towards a FA atlas and apply the estimated deformation to MD, AD & RD.

    This pipelines performs the analysis of tracts using a white matter atlas
    and computes mean value of the scalar on each tracts of this atlas. The
    pipelines registers the FA-map of a subject onto the FA-map of the atlas
    thanks to antsRegistrationSyNQuick. Then, the estimated deformation is
    applied to the MD-map, AD-map and RD-map. Finally, the labelled atlas
    is used to compute the statistics of each scalar on each tract of the white
    matter atlas.

    Args:
        working_directory (Optional[str]): Directory where the temporary
            results are stored. If not specified, it is
            automatically generated (generally in /tmp/).
        name (Optional[str]): Name of the pipelines.

    Inputnode:
        in_fa (str): FA-map of the subject in native space.
        in_md (str): MD-map of the subject in native space.
        in_ad (str): AD-map of the subject in native space.
        in_rd (str): RD-map of the subject in native space.

    Outputnode:
        out_affine_transform (str): Affine transformation matrix obtained by
            antsRegistrationSyNQuick after registration towards <atlas_name>.
        out_b_spline_transform (str): BSpline transformation obtained by
            antsRegistrationSyNQuick after registration towards <atlas_name>.
        out_norm_fa (str): FA-map registered on <atlas_name>.
        out_norm_md (str): MD-map registered on <atlas_name>.
        out_norm_ad (str): AD-map registered on <atlas_name>.
        out_norm_rd (str): RD-map registered on <atlas_name>.
    """
    import os
    import tempfile

    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from clinica.utils.atlas import AtlasAbstract, JHUDTI811mm
    from clinica.utils.check_dependency import check_environment_variable
    from clinica.utils.mri_registration import (
        apply_ants_registration_syn_quick_transformation,
    )
    from nipype.interfaces.ants import RegistrationSynQuick

    atlas = JHUDTI811mm()

    if not isinstance(atlas, AtlasAbstract):
        raise Exception("Atlas element must be an AtlasAbstract type")

    if working_directory is None:
        working_directory = tempfile.mkdtemp()

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["in_fa", "in_md", "in_ad", "in_rd", "in_atlas_scalar_image"]
        ),
        name="inputnode",
    )
    fsl_dir = check_environment_variable("FSLDIR", "FSL")
    fa_map = os.path.join(fsl_dir, "data", "atlases", "JHU", "JHU-ICBM-FA-1mm.nii.gz")
    inputnode.inputs.in_atlas_scalar_image = fa_map

    register_fa = pe.Node(interface=RegistrationSynQuick(), name="register_fa")

    apply_ants_registration_for_md = pe.Node(
        interface=niu.Function(
            input_names=[
                "in_image",
                "in_reference_image",
                "in_affine_transformation",
                "in_bspline_transformation",
                "name_output_image",
            ],
            output_names=["out_deformed_image"],
            function=apply_ants_registration_syn_quick_transformation,
        ),
        name="apply_ants_registration_for_md",
    )
    apply_ants_registration_for_md.inputs.name_output_image = (
        "space-"
        + atlas.get_name_atlas()
        + "_res-"
        + atlas.get_spatial_resolution()
        + "_MD.nii.gz"
    )
    apply_ants_registration_for_ad = pe.Node(
        interface=niu.Function(
            input_names=[
                "in_image",
                "in_reference_image",
                "in_affine_transformation",
                "in_bspline_transformation",
                "name_output_image",
            ],
            output_names=["out_deformed_image"],
            function=apply_ants_registration_syn_quick_transformation,
        ),
        name="apply_ants_registration_for_ad",
    )
    apply_ants_registration_for_ad.inputs.name_output_image = (
        "space-"
        + atlas.get_name_atlas()
        + "_res-"
        + atlas.get_spatial_resolution()
        + "_AD.nii.gz"
    )
    apply_ants_registration_for_rd = pe.Node(
        interface=niu.Function(
            input_names=[
                "in_image",
                "in_reference_image",
                "in_affine_transformation",
                "in_bspline_transformation",
                "name_output_image",
            ],
            output_names=["out_deformed_image"],
            function=apply_ants_registration_syn_quick_transformation,
        ),
        name="apply_ants_registration_for_rd",
    )
    apply_ants_registration_for_rd.inputs.name_output_image = (
        "space-"
        + atlas.get_name_atlas()
        + "_res-"
        + atlas.get_spatial_resolution()
        + "_RD.nii.gz"
    )

    thres_map = pe.Node(
        fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="RemoveNegative"
    )
    thres_fa = thres_map.clone("RemoveNegative_FA")
    thres_md = thres_map.clone("RemoveNegative_MD")
    thres_ad = thres_map.clone("RemoveNegative_AD")
    thres_rd = thres_map.clone("RemoveNegative_RD")

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "out_norm_fa",
                "out_norm_md",
                "out_norm_ad",
                "out_norm_rd",
                "out_affine_matrix",
                "out_b_spline_transform",
            ]
        ),
        name="outputnode",
    )

    wf = pe.Workflow(name=name, base_dir=working_directory)
    wf.connect(
        [
            # Registration of FA-map onto the atlas:
            (
                inputnode,
                register_fa,
                [
                    ("in_fa", "moving_image"),
                    ("in_atlas_scalar_image", "fixed_image"),
                ],
            ),
            # Apply deformation field on MD, AD & RD:
            (
                inputnode,
                apply_ants_registration_for_md,
                [("in_md", "in_image")],
            ),
            (
                inputnode,
                apply_ants_registration_for_md,
                [("in_atlas_scalar_image", "in_reference_image")],
            ),
            (
                register_fa,
                apply_ants_registration_for_md,
                [("out_matrix", "in_affine_transformation")],
            ),
            (
                register_fa,
                apply_ants_registration_for_md,
                [("forward_warp_field", "in_bspline_transformation")],
            ),
            (
                inputnode,
                apply_ants_registration_for_ad,
                [("in_ad", "in_image")],
            ),
            (
                inputnode,
                apply_ants_registration_for_ad,
                [("in_atlas_scalar_image", "in_reference_image")],
            ),
            (
                register_fa,
                apply_ants_registration_for_ad,
                [("out_matrix", "in_affine_transformation")],
            ),
            (
                register_fa,
                apply_ants_registration_for_ad,
                [("forward_warp_field", "in_bspline_transformation")],
            ),
            (
                inputnode,
                apply_ants_registration_for_rd,
                [("in_rd", "in_image")],
            ),
            (
                inputnode,
                apply_ants_registration_for_rd,
                [("in_atlas_scalar_image", "in_reference_image")],
            ),
            (
                register_fa,
                apply_ants_registration_for_rd,
                [("out_matrix", "in_affine_transformation")],
            ),
            (
                register_fa,
                apply_ants_registration_for_rd,
                [("forward_warp_field", "in_bspline_transformation")],
            ),
            # Remove negative values from the DTI maps:
            (register_fa, thres_fa, [("warped_image", "in_file")]),
            (
                apply_ants_registration_for_md,
                thres_md,
                [("out_deformed_image", "in_file")],
            ),
            (
                apply_ants_registration_for_rd,
                thres_rd,
                [("out_deformed_image", "in_file")],
            ),
            (
                apply_ants_registration_for_ad,
                thres_ad,
                [("out_deformed_image", "in_file")],
            ),
            # Outputnode:
            (thres_fa, outputnode, [("out_file", "out_norm_fa")]),
            (
                register_fa,
                outputnode,
                [
                    ("out_matrix", "out_affine_matrix"),
                    ("forward_warp_field", "out_b_spline_transform"),
                    ("inverse_warp_field", "out_inverse_warp"),
                ],
            ),
            (thres_md, outputnode, [("out_file", "out_norm_md")]),
            (thres_ad, outputnode, [("out_file", "out_norm_ad")]),
            (thres_rd, outputnode, [("out_file", "out_norm_rd")]),
        ]
    )

    return wf
