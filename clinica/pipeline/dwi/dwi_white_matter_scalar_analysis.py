#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains the DTI scalar analysis pipeline."""

def create_dti_atlas_scalar_analysis(
                in_scalar_image, atlas_labels, atlas_scalar_image, working_directory, datasink_directory, name="create_dti_atlas_scalar_analysis")):
    """
    Perform tracts analysis according to a white matter atlas using a tensor-derived scalar image.

    This function perform the analysis of tracts using a white matter atlas and compute mean value of
    the scalar on each tracts of this atlas. The function first coregister the subject scalar image
    on the equivalent scalar image of the atlas and then use the labels to computs the statistics of the
    scalar on each tracks of the whie matter atlas.

    Args:
        in_scalar_image (str): 3D image of the scalar obtained from the tensor
        atlas_labels (str): 3D Image of the white matter labels from the atlas
        atlas_scalar_image (str): 3D image of the same scalar as in "in_scalar_image" but from the atlas
        datasink_directory (str): Directory where the results are stored.
        working_directory (Optional[str]): Directory where the temporary
            results are stored. If not specified, it is automatically
            generated (generally in /tmp/).


    Outputnode:
        out_stats_file (str): File containing for each tract, the mean value of the scalar, the standard deviation
            and the nb of voxels.
    """
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    from clinica.pipeline.registration.mri_registration import antsRegistrationSyNQuick
    from clinica.pipeline.postprocessing.dwi_utils import dti_atlas_scalar_analysis
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
    datasink.inputs.base_directory = op.join(datasink_directory,'dti_scalar_analysis/')


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

    return wf
