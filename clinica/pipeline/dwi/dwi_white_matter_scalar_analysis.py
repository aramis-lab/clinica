#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the post-processing pipeline."""


def dti_based_analysis_pipeline(subject_id, session_id, caps_directory, analysis_series_id='default',
        working_directory=None, atlas_name=None, name="dti_based_analysis_pipeline"):
    """
    TODO

    Args:
        subject_id (str): Subject ID in a BIDS format ('sub-<participant_label>').
        session_id (str): Session ID in a BIDS format ('ses-<session_label>').
        analysis_series_id (str): Analysis series ID (will create the 'analysis-series-<analysis_series_id>/' folder
            for the CAPS hierarchy)
        caps_directory (str): Directory where the results are stored in a CAPS hierarchy.
        working_directory (Optional[str]): Directory where the temporary results are stored. If not specified, it is
            automatically generated (generally in /tmp/).
        name (Optional[str]): Name of the pipeline.

    Inputnode:
        TODO

    Outputnode:
        TODO
    """
    import tempfile
    import os
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    from clinica.pipeline.registration.mri_utils import ants_registration_syn_quick
    from clinica.pipeline.registration.mri_utils import apply_ants_registration_syn_quick_transformation

    try:
        fsl_dir = os.environ.get('FSLDIR', '')
        if not fsl_dir:
            raise RuntimeError('FSLDIR variable is not set')
    except Exception as e:
        print(str(e))
        exit(1)

    try:
        if fsl.Info.version().split(".") < ['5', '0', '5']:
            raise RuntimeError('FSL version must be greater than 5.0.5')
    except Exception as e:
        print(str(e))
        exit(1)

    if working_directory is None:
        working_directory = tempfile.mkdtemp()


    if atlas_name is None:
        atlas_name='JHU-ICBM-labels'
        atlas_scalar_image=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')
        atlas_labels=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', atlas_name + '-1mm.nii.gz')
#        atlas_name='JHU-ICBM-tracts-maxprob-thr25'
#        atlas_scalar_image=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')
#        atlas_labels=os.path.join(fsl_dir, 'data', 'atlases', 'JHU', atlas_name + '-1mm.nii.gz')

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
    scalar_analysis_fa.inputs.name_output_file = 'stats_fa.csv'
    scalar_analysis_md = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_md')
    scalar_analysis_md.inputs.name_output_file = 'stats_md.csv'
    scalar_analysis_ad = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_ad')
    scalar_analysis_ad.inputs.name_output_file = 'stats_ad.csv'
    scalar_analysis_rd = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'],
        output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis_rd')
    scalar_analysis_rd.inputs.name_output_file = 'stats_rd.csv'

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_stats_file', 'image_warped', 'affine_matrix', 'warp', 'inverse_warp', 'inverse_warped']),
        name='outputnode')

#    outputnode = pe.Node(niu.IdentityInterface(
#        fields=['out_fa_statistics', 'out_md_statistics', 'out_ad_statistics', 'out_rd_statistics',
#                'out_registered_fa', 'out_registered_md', 'out_registered_ad', 'out_registered_rd',
#                'out_affine_transform', 'out_b_spline_transform']),
#        name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    caps_identifier = 'sub-' + subject_id + '_ses­' + session_id
    datasink.inputs.base_directory = os.path.join(caps_directory, 'analysis-series-' + analysis_series_id, 'subjects',
                                          subject_id, session_id, 'dwi')
#                                          subject_id, session_id, 'dwi', atlas_name + '-registration')
    datasink.inputs.substitutions = [('SyN_Quick0GenericAffine.mat', caps_identifier + '_affine-transform.mat'),
                                     ('SyN_Quick1Warp.nii.gz', caps_identifier + '_affine-transform.mat'),
                                     ('SyN_QuickWarped.nii.gz', caps_identifier + '_fa-map-to-' + atlas_name + '.nii.gz'),
                                     ('md_map_registered_to_atlas.nii.gz', caps_identifier + '_md-map-to-' + atlas_name + '.nii.gz'),
                                     ('ad_map_registered_to_atlas.nii.gz', caps_identifier + '_ad-map-to-' + atlas_name + '.nii.gz'),
                                     ('rd_map_registered_to_atlas.nii.gz', caps_identifier + '_rd-map-to-' + atlas_name + '.nii.gz'),
                                     ('stats_fa.csv', caps_identifier + '_fa-map-statistics-on-' + atlas_name + '.csv'),
                                     ('stats_md.csv', caps_identifier + '_md-map-statistics-on-' + atlas_name + '.csv'),
                                     ('stats_ad.csv', caps_identifier + '_ad-map-statistics-on-' + atlas_name + '.csv'),
                                     ('stats_rd.csv', caps_identifier + '_rd-map-statistics-on-' + atlas_name + '.csv')
                                     ]

    wf = pe.Workflow(name=name)
    wf.base_dir = working_directory
    wf.connect([
        # Registration of FA-map onto the atlas:
        (inputnode, ants_registration, [('in_fa', 'moving_image'),
                                        ('in_atlas_scalar_image', 'fixe_image')]),
        # Statistics of FA:
        (inputnode,         scalar_analysis_fa, [('in_atlas_labels', 'atlas_labels_image')]),
        (ants_registration, scalar_analysis_fa, [('image_warped', 'input_image')]),
        # Apply deformation field on MD, AD & RD and compute statistics:
        (inputnode,         apply_ants_registration_for_md, [('in_md', 'in_image')]),
        (inputnode,         apply_ants_registration_for_md, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration, apply_ants_registration_for_md, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration, apply_ants_registration_for_md, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_md, [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_md, scalar_analysis_md, [('out_deformed_image', 'input_image')]),
        (inputnode,         apply_ants_registration_for_ad, [('in_ad', 'in_image')]),
        (inputnode,         apply_ants_registration_for_ad, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration, apply_ants_registration_for_ad, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration, apply_ants_registration_for_ad, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_ad, [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_ad, scalar_analysis_ad, [('out_deformed_image', 'input_image')]),
        (inputnode,         apply_ants_registration_for_rd, [('in_rd', 'in_image')]),
        (inputnode,         apply_ants_registration_for_rd, [('in_atlas_scalar_image', 'in_reference_image')]),
        (ants_registration, apply_ants_registration_for_rd, [('affine_matrix', 'in_affine_transformation')]),
        (ants_registration, apply_ants_registration_for_rd, [('warp', 'in_bspline_transformation')]),
        (inputnode,                      scalar_analysis_rd, [('in_atlas_labels', 'atlas_labels_image')]),
        (apply_ants_registration_for_rd, scalar_analysis_rd, [('out_deformed_image', 'input_image')]),
        # Outputnode:
        (ants_registration,  outputnode, [('image_warped', 'image_warped'),
                                          ('affine_matrix', 'affine_matrix'),
                                          ('warp', 'warp'),
                                          ('inverse_warp', 'inverse_warp'),
                                          ('inverse_warped', 'inverse_warped')]),
        (scalar_analysis_fa, outputnode, [('outfile', 'out_stats_file')]),
        # Saving files with datasink:
        (ants_registration,  datasink, [('image_warped', 'atlas-registration.@image_warped'),
                                        ('affine_matrix', 'atlas-registration.@affine_matrix'),
                                        ('warp', 'atlas-registration.@warp'),
                                        ('inverse_warp', 'inverse_warp'),
                                        ('inverse_warped', 'inverse_warped')]),
        (apply_ants_registration_for_md, datasink, [('out_deformed_image', 'atlas-registration.@out_deformed_md')]),
        (apply_ants_registration_for_ad, datasink, [('out_deformed_image', 'atlas-registration.@out_deformed_ad')]),
        (apply_ants_registration_for_rd, datasink, [('out_deformed_image', 'atlas-registration.@out_deformed_rd')]),
        (scalar_analysis_fa, datasink, [('outfile', 'atlas-registration.@out_stats_file_fa')]),
        (scalar_analysis_md, datasink, [('outfile', 'atlas-registration.@out_stats_file_md')]),
        (scalar_analysis_ad, datasink, [('outfile', 'atlas-registration.@out_stats_file_ad')])
    ])

    return wf




def dti_atlas_scalar_analysis_pipeline(
                in_scalar_image, atlas_labels, atlas_scalar_image, datasink_directory, working_directory=None, name="dti_atlas_scalar_analysis_pipeline"):
    """
    Perform tracts analysis according to a white matter atlas using a tensor-derived scalar image.

    This function performs the analysis of tracts using a white matter atlas and compute mean value of the scalar on each tracts of this atlas. The function first coregister the subject scalar image on the equivalent scalar image of the atlas and then use the labels to computes the statistics of the scalar on each tracks of the white matter atlas.

    Args:
        in_scalar_image (str): 3D image of the scalar obtained from the tensor
        atlas_labels (str): 3D Image of the white matter labels from the atlas
        atlas_scalar_image (str): 3D image of the same scalar as in "in_scalar_image" but from the atlas
        datasink_directory (str): Directory where the results are stored.
        working_directory (Optional[str]): Directory where the temporary results are stored. If not specified, it is automatically generated (generally in /tmp/).

    Outputnode:
        out_stats_file (str): File containing for each tract, the mean value of the scalar, the standard deviation and the nb of voxels.
    """
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    from clinica.pipeline.registration.mri_registration import antsRegistrationSyNQuick
    import tempfile

    if working_directory is None:
        working_directory = tempfile.mkdtemp()

    inputs=[in_scalar_image, atlas_labels, atlas_scalar_image, working_directory, datasink_directory]

    for input_file in inputs:
        if not op.exists(input_file):
            raise IOError('file {} does not exist'.format(input_file))

    datasource = pe.Node(interface=nio.DataGrabber(infields=[], outfields=['in_scalar_image', 'atlas_labels', 'atlas_scalar_image']), name='datasource')
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(in_scalar_image= in_scalar_image,
                                            atlas_labels=atlas_labels,
                                            atlas_scalar_image=atlas_scalar_image)
    datasource.inputs.template_args = dict(in_scalar_image=[[]],
                                           atlas_labels=[[]],
                                           atlas_scalar_image=[[]])
    datasource.inputs.sort_filelist = True

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_scalar_image', 'atlas_labels', 'atlas_scalar_image']),
        name='inputnode')

    antsRegistrationSyNQuick = pe.Node(interface=niu.Function(
        input_names=['fixe_image', 'moving_image'], output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
        function=antsRegistrationSyNQuick), name='antsRegistrationSyNQuick')

    scalar_analysis = pe.Node(interface=niu.Function(
        input_names=['input_image', 'atlas_labels_image', 'name_output_file'], output_names=['outfile'],
        function=dti_atlas_scalar_analysis), name='scalar_analysis')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_stats_file', 'image_warped', 'affine_matrix', 'warp', 'inverse_warp', 'inverse_warped']),
        name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'dti_scalar_analysis/')


    wf = pe.Workflow(name='dti_scalar_analysis')
    wf.base_dir = working_directory

    wf.connect([
        (datasource, inputnode, [('in_scalar_image', 'in_scalar_image'),
                                 ('atlas_labels', 'atlas_labels'),
                                 ('atlas_scalar_image', 'atlas_scalar_image')]),
        (inputnode, antsRegistrationSyNQuick, [('in_scalar_image', 'moving_image'),
                                               ('atlas_scalar_image', 'fixe_image')]),
        (inputnode, scalar_analysis, [('atlas_labels', 'atlas_labels_image')]),
        (antsRegistrationSyNQuick, scalar_analysis, [('image_warped', 'input_image')]),
        (antsRegistrationSyNQuick, outputnode, [('image_warped', 'image_warped'),
                                                ('affine_matrix', 'affine_matrix'),
                                                ('warp', 'warp'),
                                                ('inverse_warp', 'inverse_warp'),
                                                ('inverse_warped', 'inverse_warped')]),
        (scalar_analysis,          outputnode, [('outfile', 'out_stats_file')]),
        (antsRegistrationSyNQuick, datasink, [('image_warped', 'image_warped'),
                                              ('affine_matrix', 'affine_matrix'),
                                              ('warp', 'warp'),
                                              ('inverse_warp', 'inverse_warp'),
                                              ('inverse_warped', 'inverse_warped')]),
        (scalar_analysis,          datasink, [('outfile', 'out_stats_file')])
    ])

    return wf



def dti_atlas_scalar_analysis(input_image, atlas_labels_image, name_output_file=None):
    """
    Compute statistics.

    Given a set a labeled tracts, this function compute statistics of a scalar image from DTI (e.g. FA, MD, etc.) in each tract.

    Args:
        input_image (str): File containing a scalar image (e.g. FA, MD, etc.).
        atlas_labels_image (str): File containing labels. These labels are used to compute statistics
        name_output_file (Optional[str]): Name of the output statistics file (default=scalar_stats.csv).

    Returns:
        outfile (str): CSV files containing the statistics (content of the columns: #Label, Mean, Standard deviation, #Voxels).
    """

    import nibabel as nib
    import numpy as np
    import pandas as pd
    import os.path as op

    if name_output_file is None:
        outfile = op.abspath('scalar_stats.csv')
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
        labeled_voxel = labeled_voxel = scalar_image_data[atlas_label_index[0,:], atlas_label_index[1,:], atlas_label_index[2,:]]
        mean_scalar = labeled_voxel.mean()
        stats_scalar[index, 1] = mean_scalar
        std_scalar = labeled_voxel.std()
        stats_scalar[index, 2] = std_scalar

    columns = np.array(['Label', 'Mean scalar', 'Std scalar', 'Nb of voxel'])
    data = pd.DataFrame(stats_scalar, columns=columns)

    data.to_csv(outfile, sep=',', index=False)

    return outfile
