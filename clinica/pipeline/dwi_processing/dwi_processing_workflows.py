

def register_dti_maps_on_atlas(working_directory=None, name="register_dti_maps_on_atlas"):
    """
    Register FA-map on a subject towards a FA atlas and apply the deformation to MD, AD & RD.

    This pipeline performs the analysis of tracts using a white matter atlas and compute mean value of the scalar
    on each tracts of this atlas. The pipeline registers the FA-map of a subject onto the FA-map of the atlas thanks to
    antsRegistrationSyNQuick. Then, the estimated deformation is applied to the MD-map, AD-map and RD-map. Finally,
    the labelled atlas is used to compute the statistics of each scalar on each tract of the white matter atlas.

    Args:
        working_directory (Optional[str]): Directory where the temporary results are stored. If not specified, it is
            automatically generated (generally in /tmp/).
        atlas(str): An instance of AbstractClass with a FA map in the get_name() method.
        name (Optional[str]): Name of the pipeline.

    Inputnode:
        in_fa (str): FA-map of the subject in native space.
        in_md (str): MD-map of the subject in native space.
        in_ad (str): AD-map of the subject in native space.
        in_rd (str): RD-map of the subject in native space.

    Outputnode:
        out_affine_transform (str): Affine transformation matrix obtained by antsRegistrationSyNQuick
            after registration towards <atlas_name>.
        out_b_spline_transform (str): BSpline transformation obtained by antsRegistrationSyNQuick
            after registration towards <atlas_name>.
        out_registered_fa (str): FA-map registered on <atlas_name>.
        out_registered_md (str): MD-map registered on <atlas_name>.
        out_registered_ad (str): AD-map registered on <atlas_name>.
        out_registered_rd (str): RD-map registered on <atlas_name>.
    """
    import tempfile
    import os
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.mri_registration import ants_registration_syn_quick
    from clinica.utils.mri_registration import apply_ants_registration_syn_quick_transformation

    from clinica.utils.atlas import JHUTracts0_1mm
    atlas = JHUTracts0_1mm()

    if not isinstance(atlas, AtlasAbstract):
        raise Exception("Atlas element must be an AtlasAbstract type")

    if working_directory is None:
        working_directory = tempfile.mkdtemp()

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_fa', 'in_md', 'in_ad', 'in_rd', 'in_atlas_scalar_image']),
        name='inputnode')
    inputnode.inputs.in_atlas_scalar_image = atlas.get_atlas_map()

    ants_registration = pe.Node(interface=niu.Function(
        input_names=['fixe_image', 'moving_image', 'prefix_output'],
        output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
        function=ants_registration_syn_quick), name='ants_registration')

    apply_ants_registration_for_md = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation',
                     'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_md')
    apply_ants_registration_for_md.inputs.name_output_image = 'md_map_registered_to_atlas.nii.gz'
    apply_ants_registration_for_ad = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation',
                     'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_ad')
    apply_ants_registration_for_ad.inputs.name_output_image = 'ad_map_registered_to_atlas.nii.gz'
    apply_ants_registration_for_rd = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation',
                     'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_rd')
    apply_ants_registration_for_rd.inputs.name_output_image = 'rd_map_registered_to_atlas.nii.gz'

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_registered_fa', 'out_registered_md', 'out_registered_ad', 'out_registered_rd',
                'out_affine_transform', 'out_b_spline_transform']),
        name='outputnode')

    wf = pe.Workflow(name=name, base_dir=working_directory)
    wf.connect([
        # Registration of FA-map onto the atlas:
        (inputnode,                      ants_registration,     [('in_fa', 'moving_image'),
                                                                 ('in_atlas_scalar_image', 'fixe_image')]),
        # Apply deformation field on MD, AD & RD and compute statistics:
        (inputnode,            apply_ants_registration_for_md, [('in_md', 'in_image')]),
        (inputnode,            apply_ants_registration_for_md, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,    apply_ants_registration_for_md, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,    apply_ants_registration_for_md, [('warp', 'in_bspline_transformation')]),

        (inputnode,            apply_ants_registration_for_ad, [('in_ad', 'in_image')]),
        (inputnode,            apply_ants_registration_for_ad, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,    apply_ants_registration_for_ad, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,    apply_ants_registration_for_ad, [('warp', 'in_bspline_transformation')]),

        (inputnode,            apply_ants_registration_for_rd, [('in_rd', 'in_image')]),
        (inputnode,            apply_ants_registration_for_rd, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,    apply_ants_registration_for_rd, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,    apply_ants_registration_for_rd, [('warp', 'in_bspline_transformation')]),
        # Outputnode:
        (ants_registration,                        outputnode,  [('image_warped', 'out_registered_fa'),
                                                                 ('affine_matrix', 'out_affine_matrix'),
                                                                 ('warp',         'out_b_spline_transform'),
                                                                 ('inverse_warp',      'out_inverse_warp')]),
        (apply_ants_registration_for_md,           outputnode,  [('out_deformed_image', 'out_registered_md')]),
        (apply_ants_registration_for_ad,           outputnode,  [('out_deformed_image', 'out_registered_ad')]),
        (apply_ants_registration_for_rd,           outputnode,  [('out_deformed_image', 'out_registered_rd')])
    ])

    return wf
