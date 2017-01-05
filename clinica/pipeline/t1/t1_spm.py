

import csv
from clinica.pipeline.t1.t1_spm_workflows import segmentation_pipeline, dartel_pipeline, t1_spm_full_pipeline
from clinica.pipeline.t1.t1_spm_utils import group_nested_images_by_subject, group_images_list_by_subject
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import pandas as pd
import os.path as op
from os import makedirs
import errno
from shutil import copyfile


def datagrabber_t1_spm_segment_pipeline(input_directory,
                                        output_directory,
                                        subjects_visits_tsv,
                                        analysis_series_id='default',
                                        working_directory=None,
                                        tissue_classes=[1, 2, 3],
                                        dartel_tissues=[1],
                                        save_warped_unmodulated=False,
                                        save_warped_modulated=False,
                                        in_write_deformation_fields=None):
    """

    T1 segmentation workflow with DataGrabber as input.


    :param input_directory: Directory where the input NIFTI images are stored
    :param output_directory: Directory to save the resulting images
    :param working_directory: Temporary directory to run the workflow
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'

    NewSegment parameters
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :param in_write_deformation_fields: Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]

    :return: SPM Segmentation workflow
    """

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)


    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = input_directory
    selectfiles.inputs.template = '%s/%s/anat/%s_%s_T1w.nii*'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.session = sessions
    selectfiles.inputs.subject_repeat = subjects
    selectfiles.inputs.session_repeat = sessions
    selectfiles.inputs.sort_filelist = False

    # Creating T1 preprocessing workflow
    seg_prep_wf = segmentation_pipeline(working_directory=working_directory,
                                        tissue_classes=tissue_classes,
                                        dartel_tissues=dartel_tissues,
                                        save_warped_unmodulated=save_warped_unmodulated,
                                        save_warped_modulated=save_warped_modulated,
                                        in_write_deformation_fields=in_write_deformation_fields)

    datasink_infields = ['native_space', 'dartel_input']
    datasink_connections = [(('new_segment.native_class_images', group_nested_images_by_subject), 'native_space'),
                            (('new_segment.dartel_input_images', group_nested_images_by_subject), 'dartel_input')]

    if save_warped_unmodulated:
        datasink_connections.append((('new_segment.normalized_class_images', group_nested_images_by_subject), 'normalized'))
        datasink_infields.append('normalized')

    if save_warped_modulated:
        datasink_connections.append((('new_segment.modulated_class_images', group_nested_images_by_subject), 'modulated_normalized'))
        datasink_infields.append('modulated_normalized')

    if in_write_deformation_fields is not None:
        if in_write_deformation_fields[0]:
            datasink_connections.append(('new_segment.inverse_deformation_field', 'inverse_deformation_field'))
            datasink_infields.append('inverse_deformation_field')
        if in_write_deformation_fields[1]:
            datasink_connections.append(('new_segment.forward_deformation_field', 'forward_deformation_field'))
            datasink_infields.append('forward_deformation_field')

    datasink_iterfields = ['container'] + datasink_infields

    print 'Datasink infields'
    print datasink_infields
    print
    print 'Datasink iterfields'
    print datasink_iterfields

    datasink = pe.MapNode(nio.DataSink(infields=datasink_infields), name='datasink', iterfield=datasink_iterfields)
    datasink.inputs.parameterization = False
    datasink.inputs.base_directory = output_directory
    datasink.inputs.container = ['analysis-series-' + analysis_series_id + '/subjects/sub-' + subjects[i] + '/ses-' + sessions[i] + '/t1/spm/segmentation' for i in range(len(subjects))]

    seg_wf = pe.Workflow(name='seg_wf')
    seg_wf.base_dir = working_directory
    seg_wf.connect([
        (selectfiles, seg_prep_wf, [('out_files', 'unzip.in_file')]),
        (seg_prep_wf, datasink, datasink_connections)
    ])

    return seg_wf


def datagrabber_t1_spm_full_pipeline(input_directory,
                                     output_directory,
                                     subjects_visits_tsv,
                                     group_id,
                                     analysis_series_id='default',
                                     working_directory=None,
                                     tissue_classes=[1, 2, 3],
                                     dartel_tissues=[1],
                                     save_warped_unmodulated=False,
                                     save_warped_modulated=False,
                                     in_write_deformation_fields=None,
                                     in_fwhm=None,
                                     in_modulate=True,
                                     in_voxel_size=None):
    """

    T1 preprocessing workflow with DataGrabber as input.

    :param input_directory: Directory where the input NIFTI images are stored
    :param output_directory: Directory to save the resulting images
    :param working_directory: Temporary directory to run the workflow
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'

    NewSegment parameters
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :param in_write_deformation_fields: Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]

    DARTELNorm2MNI parameters
    :param in_fwhm: A list of 3 floats specifying the FWHM for each dimension
    :param in_modulate: A boolean. Modulate output images - no modulation preserves concentrations
    :param in_voxel_size: A list of 3 floats specifying voxel sizes for each dimension of output image

    :return: SPM Segmentation and registration workflow
    """

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)


    dest_group_file = op.join(output_directory, 'analysis-series-' + analysis_series_id + '/group-' + group_id + '/group-' + group_id + '_subjects_visits_list.tsv')

    if op.isfile(dest_group_file):

        existing_subjects = []
        existing_sessions = []
        with open(dest_group_file, 'rb') as tsvin:
            tsv_reader = csv.reader(tsvin, delimiter='\t')

            for row in tsv_reader:
                existing_subjects.append(row[0])
                existing_sessions.append(row[1])

        if subjects != existing_subjects or sessions != existing_sessions:
            raise Exception('A different list of subjects and visits already exists for the same group definition: ' + group_id)

    else:
        try:
            makedirs(op.dirname(dest_group_file))

        except OSError as exc:
            if exc.errno == errno.EEXIST and op.isdir(op.dirname(dest_group_file)):
                pass
            else:
                raise

        copyfile(subjects_visits_tsv, dest_group_file)

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                          outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = input_directory
    selectfiles.inputs.template = '%s/%s/anat/%s_%s_T1w.nii*'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.session = sessions
    selectfiles.inputs.subject_repeat = subjects
    selectfiles.inputs.session_repeat = sessions
    selectfiles.inputs.sort_filelist = False

    t1_spm_prep_wf = t1_spm_full_pipeline(working_directory=working_directory,
                                          name='t1_spm_full_wf',
                                          tissue_classes=tissue_classes,
                                          dartel_tissues=dartel_tissues,
                                          save_warped_unmodulated=save_warped_unmodulated,
                                          save_warped_modulated=save_warped_modulated,
                                          in_write_deformation_fields=in_write_deformation_fields,
                                          in_fwhm=in_fwhm,
                                          in_modulate=in_modulate,
                                          in_voxel_size=in_voxel_size)

    seg_datasink_infields = ['native_space', 'dartel_input']
    seg_datasink_connections = [(('segmentation_wf.new_segment.native_class_images', group_nested_images_by_subject), 'native_space'), (('segmentation_wf.new_segment.dartel_input_images', group_nested_images_by_subject), 'dartel_input')]

    if save_warped_unmodulated:
        seg_datasink_connections.append((('segmentation_wf.new_segment.normalized_class_images', group_nested_images_by_subject), 'normalized'))
        seg_datasink_infields.append('normalized')

    if save_warped_modulated:
        seg_datasink_connections.append((('segmentation_wf.new_segment.modulated_class_images', group_nested_images_by_subject), 'modulated_normalized'))
        seg_datasink_infields.append('modulated_normalized')

    if in_write_deformation_fields is not None:
        if in_write_deformation_fields[0]:
            seg_datasink_connections.append(('segmentation_wf.new_segment.inverse_deformation_field', 'inverse_deformation_field'))
            seg_datasink_infields.append('inverse_deformation_field')
        if in_write_deformation_fields[1]:
            seg_datasink_connections.append(('segmentation_wf.new_segment.forward_deformation_field', 'forward_deformation_field'))
            seg_datasink_infields.append('forward_deformation_field')

    seg_datasink_iterfields = ['container'] + seg_datasink_infields

    print 'Datasink infields'
    print seg_datasink_infields
    print
    print 'Datasink iterfields'
    print seg_datasink_iterfields

    seg_datasink = pe.MapNode(nio.DataSink(infields=seg_datasink_infields), name='datasink', iterfield=seg_datasink_iterfields)
    seg_datasink.inputs.parameterization = False
    seg_datasink.inputs.base_directory = output_directory
    seg_datasink.inputs.container = ['analysis-series-' + analysis_series_id + '/subjects/sub-' + subjects[i] + '/ses-' + sessions[i] + '/t1/spm/segmentation' for i in range(len(subjects))]

    template_datasink = pe.Node(nio.DataSink(), name='template_datasink')
    template_datasink.inputs.parameterization = False
    template_datasink.inputs.base_directory = op.join(output_directory, 'analysis-series-' + analysis_series_id + '/group-' + group_id + '/t1/spm')
    # template_datasink.container = 'analysis-series-' + analysis_series_id + '/group-' + group_id + '/t1/spm'

    dartel_datasink = pe.MapNode(nio.DataSink(infields=['flow_fields', 'registered']), name='dartel_datasink', iterfield=['container', 'flow_fields', 'registered'])
    dartel_datasink.inputs.parameterization = False
    dartel_datasink.inputs.base_directory = output_directory
    dartel_datasink.inputs.container = ['analysis-series-' + analysis_series_id + '/subjects/sub-' + subjects[i] + '/ses-' + sessions[i] + '/t1/spm/dartel/group-' + group_id for i in range(len(subjects))]

    preproc_wf = pe.Workflow(name='preproc_wf')
    if working_directory is not None:
        preproc_wf.base_dir = working_directory

    preproc_wf.connect([
        (selectfiles, t1_spm_prep_wf, [('out_files', 'segmentation_wf.unzip.in_file')]),
        # (selectfiles, t1_spm_prep_wf, [('out_files', 'segmentation_wf.new_segment.channel_files')]),
        (t1_spm_prep_wf, seg_datasink, seg_datasink_connections),
        (t1_spm_prep_wf, template_datasink, [('outputnode.out_final_template_file', 'final_template'),
                                             ('outputnode.out_template_files', 'all_templates')]),
        (t1_spm_prep_wf, dartel_datasink, [('dartel_wf.dartelTemplate.dartel_flow_fields', 'flow_fields'),
                                           (('dartel_wf.dartel2MNI.normalized_files', group_images_list_by_subject, tissue_classes), 'registered')])
    ])

    print seg_datasink.outputs

    return preproc_wf
