#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains tools for the normalization of DTI-scalar measures onto an atlas with labelled tracts."""




def dti_based_analysis_pipeline(
        participant_id, session_id, caps_directory,
        working_directory=None, atlas_name=None, name="dti_based_analysis_pipeline"):
    """
    Perform tracts analysis according to a white matter atlas using a tensor-derived scalar image.

    This pipeline performs the analysis of tracts using a white matter atlas and compute mean value of the scalar
    on each tracts of this atlas. The pipeline registers the FA-map of a subject onto the FA-map of the atlas thanks to
    antsRegistrationSyNQuick. Then, the estimated deformation is applied to the MD-map, AD-map and RD-map. Finally,
    the labelled atlas is used to compute the statistics of each scalar on each tract of the white matter atlas.

    Args:
        participant_id (str): Subject (participant) ID in a BIDS format ('sub-<participant_label>').
        session_id (str): Session ID in a BIDS format ('ses-<session_label>').
        caps_directory (str): Directory where the results are stored in a CAPS hierarchy.
        working_directory (Optional[str]): Directory where the temporary results are stored. If not specified, it is
            automatically generated (generally in /tmp/).
        atlas_name(str): Name of the white matter atlas.
        name (Optional[str]): Name of the pipeline.

    Inputnode:
        in_fa (str): FA-map of the subject.
        in_md (str): MD-map of the subject.
        in_ad (str): AD-map of the subject.
        in_rd (str): RD-map of the subject.

    Outputnode:
        out_affine_transform (str): Affine transformation matrix obtained by antsRegistrationSyNQuick
            after registration towards <atlas_name>.
        out_b_spline_transform (str): BSpline transformation obtained by antsRegistrationSyNQuick
            after registration towards <atlas_name>.
        out_registered_fa (str): FA-map registered on <atlas_name>.
        out_registered_md (str): MD-map registered on <atlas_name>.
        out_registered_ad (str): AD-map registered on <atlas_name>.
        out_registered_rd (str): RD-map registered on <atlas_name>.
        out_fa_statistics (str): Statistics of the FA-map on <atlas_name>.
        out_md_statistics (str): Statistics of the MD-map on <atlas_name>.
        out_ad_statistics (str): Statistics of the AD-map on <atlas_name>.
        out_rd_statistics (str): Statistics of the RD-map on <atlas_name>.
    """
    import tempfile
    import os
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from clinica.utils.mri_registration import ants_registration_syn_quick
    from clinica.utils.mri_registration import apply_ants_registration_syn_quick_transformation

    try:
        fsl_dir = os.environ.get('FSLDIR', '')
        if not fsl_dir:
            raise RuntimeError('FSLDIR variable is not set')
    except Exception as e:
        print(str(e))
        exit(1)

    if working_directory is None:
        working_directory = tempfile.mkdtemp()


    if atlas_name is None or atlas_name.lower() == 'jhu-icbm-tracts-maxprob-thr25' :
        atlas_name='JHU-ICBM-tracts-maxprob-thr25'
        atlas_scalar_image=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')
        atlas_labels=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', atlas_name + '-1mm.nii.gz')
    elif atlas_name.lower() == 'jhu-icbm-labels':
        atlas_name='JHU-ICBM-labels'
        atlas_scalar_image=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')
        atlas_labels=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', atlas_name + '-1mm.nii.gz')
    else:
        raise IOError('atlas_name should be jhu-icbm-tracts-maxprob-thr25 or jhu-icbm-labels (it is ' + atlas_name.lower()  + ')')


    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_fa', 'in_md', 'in_ad', 'in_rd', 'in_atlas_scalar_image', 'in_atlas_labels']),
        name='inputnode')
    inputnode.inputs.in_atlas_scalar_image=atlas_scalar_image
    inputnode.inputs.in_atlas_labels=atlas_labels

    ants_registration = pe.Node(interface=niu.Function(
        input_names=['fixe_image', 'moving_image', 'prefix_output'],
        output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
        function=ants_registration_syn_quick), name='ants_registration')

    apply_ants_registration_for_md = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation', 'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_md')
    apply_ants_registration_for_md.inputs.name_output_image = 'md_map_registered_to_atlas.nii.gz'
    apply_ants_registration_for_ad = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation', 'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_ad')
    apply_ants_registration_for_ad.inputs.name_output_image = 'ad_map_registered_to_atlas.nii.gz'
    apply_ants_registration_for_rd = pe.Node(interface=niu.Function(
        input_names=['in_image', 'in_reference_image', 'in_affine_transformation', 'in_bspline_transformation', 'name_output_image'],
        output_names=['out_deformed_image'],
        function=apply_ants_registration_syn_quick_transformation), name='apply_ants_registration_for_rd')
    apply_ants_registration_for_rd.inputs.name_output_image = 'rd_map_registered_to_atlas.nii.gz'

    scalar_analysis_fa = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_fa')
    scalar_analysis_fa.inputs.name_output_file = 'stats_fa.tsv'
    scalar_analysis_md = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_md')
    scalar_analysis_md.inputs.name_output_file = 'stats_md.tsv'
    scalar_analysis_ad = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_ad')
    scalar_analysis_ad.inputs.name_output_file = 'stats_ad.tsv'
    scalar_analysis_rd = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_rd')
    scalar_analysis_rd.inputs.name_output_file = 'stats_rd.tsv'

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_fa_statistics', 'out_md_statistics', 'out_ad_statistics', 'out_rd_statistics',
                'out_registered_fa', 'out_registered_md', 'out_registered_ad', 'out_registered_rd',
                'out_affine_transform', 'out_b_spline_transform']),
        name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    caps_identifier = participant_id + '_' + session_id
    datasink.inputs.base_directory = os.path.join(caps_directory, 'subjects', participant_id, session_id, 'dwi')
    datasink.inputs.substitutions = [('SyN_Quick0GenericAffine.mat', caps_identifier + '_transform-affine_' + atlas_name + '.mat'),
                                     ('SyN_Quick1Warp.nii.gz', caps_identifier + '_transform-bspline_' + atlas_name + '.nii.gz'),
                                     ('SyN_QuickWarped.nii.gz', caps_identifier + '_map-fa_registeredOn' + atlas_name + '.nii.gz'),
                                     ('md_map_registered_to_atlas.nii.gz', caps_identifier + '_map-md_registeredOn' + atlas_name + '.nii.gz'),
                                     ('ad_map_registered_to_atlas.nii.gz', caps_identifier + '_map-ad_registeredOn' + atlas_name + '.nii.gz'),
                                     ('rd_map_registered_to_atlas.nii.gz', caps_identifier + '_map-rd_registeredOn' + atlas_name + '.nii.gz'),
                                     ('stats_fa.csv', caps_identifier + '_map-fa_statisticsOn' + atlas_name + '.tsv'),
                                     ('stats_md.csv', caps_identifier + '_map-md_statisticsOn' + atlas_name + '.tsv'),
                                     ('stats_ad.csv', caps_identifier + '_map-ad_statisticsOn' + atlas_name + '.tsv'),
                                     ('stats_rd.csv', caps_identifier + '_map-rd_statisticsOn' + atlas_name + '.tsv')
                                     ]

    wf = pe.Workflow(name=name, base_dir=working_directory)
    wf.connect([
        # Registration of FA-map onto the atlas:
        (inputnode,                      ants_registration,              [('in_fa', 'moving_image'),
                                                                         ('in_atlas_scalar_image', 'fixe_image')]),
        # Statistics of FA:
        (inputnode,                                  scalar_analysis_fa, [('in_atlas_labels', 'atlas_labels_image')]),
        (ants_registration,                          scalar_analysis_fa, [('image_warped', 'input_image')]),
        # Apply deformation field on MD, AD & RD and compute statistics:
        (inputnode,                      apply_ants_registration_for_md, [('in_md', 'in_image')]),
        (inputnode,                      apply_ants_registration_for_md, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,              apply_ants_registration_for_md, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,              apply_ants_registration_for_md, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_md,             [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_md, scalar_analysis_md,             [('out_deformed_image', 'input_image')]),

        (inputnode,                      apply_ants_registration_for_ad, [('in_ad', 'in_image')]),
        (inputnode,                      apply_ants_registration_for_ad, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,              apply_ants_registration_for_ad, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,              apply_ants_registration_for_ad, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_ad,              [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_ad, scalar_analysis_ad,              [('out_deformed_image', 'input_image')]),

        (inputnode,                      apply_ants_registration_for_rd, [('in_rd', 'in_image')]),
        (inputnode,                      apply_ants_registration_for_rd, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration,              apply_ants_registration_for_rd, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration,              apply_ants_registration_for_rd, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_rd,             [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_rd, scalar_analysis_rd,             [('out_deformed_image', 'input_image')]),
        # Outputnode:
        (ants_registration,                                 outputnode,  [('image_warped', 'out_deformed_fa'),
                                                                          ('affine_matrix', 'out_affine_matrix'),
                                                                          ('warp',         'out_b_spline_transform'),
                                                                          ('inverse_warp',      'out_inverse_warp')]),
        (apply_ants_registration_for_md,                    outputnode,  [('out_deformed_image', 'out_deformed_md')]),
        (apply_ants_registration_for_ad,                    outputnode,  [('out_deformed_image', 'out_deformed_ad')]),
        (apply_ants_registration_for_rd,                    outputnode,  [('out_deformed_image', 'out_deformed_rd')]),
        (scalar_analysis_fa,                                outputnode,  [('outfile',           'out_stats_file_fa')]),
        (scalar_analysis_md,                                outputnode,  [('outfile',           'out_stats_file_md')]),
        (scalar_analysis_ad,                                outputnode,  [('outfile',           'out_stats_file_ad')]),
        (scalar_analysis_rd,                                outputnode,  [('outfile',           'out_stats_file_rd')]),
        # Saving files with datasink:
        (ants_registration,                                    datasink, [('image_warped',              'atlas-registration.@out_deformed_fa'),
                                                                          ('affine_matrix',           'atlas-registration.@out_affine_matrix'),
                                                                          ('warp',               'atlas-registration.@out_b_spline_transform')]),
        (apply_ants_registration_for_md,                       datasink, [('out_deformed_image',        'atlas-registration.@out_deformed_md')]),
        (apply_ants_registration_for_ad,                       datasink, [('out_deformed_image',        'atlas-registration.@out_deformed_ad')]),
        (apply_ants_registration_for_rd,                       datasink, [('out_deformed_image',        'atlas-registration.@out_deformed_rd')]),
        (scalar_analysis_fa,                                   datasink, [('outfile',                   'atlas-statistics.@out_stats_file_fa')]),
        (scalar_analysis_md,                                   datasink, [('outfile',                   'atlas-statistics.@out_stats_file_md')]),
        (scalar_analysis_ad,                                   datasink, [('outfile',                   'atlas-statistics.@out_stats_file_ad')]),
        (scalar_analysis_rd,                                   datasink, [('outfile',                   'atlas-statistics.@out_stats_file_rd')])
    ])

    return wf





def dti_atlas_scalar_analysis(input_image, atlas_labels_image, name_output_file=None):
    """
    Compute statistics.

    Given a set a labeled tracts, this function compute statistics of a scalar image from DTI (e.g. FA, MD, etc.)
    in each tract.

    Args:
        input_image (str): File containing a scalar image (e.g. FA, MD, etc.).
        atlas_labels_image (str): File containing labels. These labels are used to compute statistics
        name_output_file (Optional[str]): Name of the output statistics file (default=scalar_stats.tsv).

    Returns:
        outfile (str): TSV files containing the statistics (content of the columns:
            label, mean scalar, std of the scalar', number of voxels).
    """

    import nibabel as nib
    import numpy as np
    import pandas
    import os.path as op

    if name_output_file is None:
        outfile = op.abspath('scalar_stats.tsv')
    else:
        outfile = op.abspath(name_output_file)

    dti_atlas = nib.load(atlas_labels_image)
    atlas_image_data = dti_atlas.get_data()

    labels = list(set(atlas_image_data.ravel()))
    stats_scalar = np.zeros((len(labels),4))

    in_image = nib.load(input_image)
    scalar_image_data = in_image.get_data()

    for index, label in enumerate(labels):
        stats_scalar[index, 0] = label
        atlas_label_index = np.array(np.where(atlas_image_data==label))
        nb_voxel = atlas_label_index.shape[1]
        stats_scalar[index, 3] = nb_voxel
        labeled_voxel = scalar_image_data[atlas_label_index[0,:], atlas_label_index[1,:], atlas_label_index[2,:]]
        mean_scalar = labeled_voxel.mean()
        stats_scalar[index, 1] = mean_scalar
        std_scalar = labeled_voxel.std()
        stats_scalar[index, 2] = std_scalar

    columns = np.array(['label', 'mean_scalar', 'std_scalar', 'nb_of_voxel'])

    data = pandas.DataFrame(stats_scalar, columns=columns, encoding='utf8')

    data.to_csv(outfile, sep='\t', index=False)

    return outfile
