

import csv
from clinica.pipeline.t1.t1_spm_workflows import segmentation_pipeline, dartel_pipeline, t1_spm_full_pipeline
from clinica.pipeline.t1.t1_spm_utils import group_images_by_subject
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import os.path as op


def datagrabber_t1_spm_segment_pipeline_bids(input_directory,
                                        output_directory,
                                        subjects_visits_tsv,
                                        analysis_series_id,
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

    subjects = []
    sessions = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            subjects.append(row[0])
            sessions.append(row[1])

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = input_directory
    selectfiles.inputs.template = 'sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii'
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

    # def segment_to_caps(output_directory, analysis_series_id, group_id, subject_id, visit_id, native_space, dartel_input, normalized=[], modulated_normalized=[], inverse_deformation_field=[], forward_deformation_field=[]):
    #     """
    #     Function receives a list of lists of native class images and a list of flow fields.
    #     It returns two lists of the same length to give as input to DARTEL2MNI
    #     This will allow to process each pair of files in parallel.
    #     :param native_space_images: list of lists of native space images
    #     :param flowfield_files: list of flow fields files
    #     :return: expanded list of native images,list of flow fields files of the same length
    #     """
    #
    #     native_files = [p in native_space]
    #
    #     return
    #
    # write_segment_to_caps = pe.Node(niu.Function(input_names=['output_directory', 'analysis_series_id', 'group_id',
    #                                                           'subject_id', 'visit_id', 'native_space', 'dartel_input',
    #                                                           'normalized', 'modulated_normalized',
    #                                                           'inverse_deformation_field', 'forward_deformation_field'],
    #                                              function=segment_to_caps),
    #                                 name='write_segment_to_caps',
    #                                 iterfield=['subject_id', 'visit_id'])
    #
    # write_segment_to_caps.inputs.output_directory = output_directory
    # write_segment_to_caps.inputs.analysis_series_id = analysis_series_id
    # write_segment_to_caps.inputs.group_id = group_id
    # write_segment_to_caps.inputs.subject_id = subjects
    # write_segment_to_caps.inputs.visit_id = sessions

    # connections_segment_to_caps = [(('outputnode.out_native_class_images', get_class_images, tissue_classes), 'native_space'),
    #                                 (('outputnode.out_dartel_input_images', get_class_images, dartel_tissues), 'dartel_input')]
    # if save_warped_unmodulated:
    #     connections_segment_to_caps.append((('outputnode.out_normalized_class_images', get_class_images, tissue_classes), 'normalized'))
    #     datasink_connections.append(('normalized', 'normalized'))
    #
    # if save_warped_modulated:
    #     connections_segment_to_caps.append((('outputnode.out_modulated_class_images', get_class_images, tissue_classes), 'modulated_normalized'))
    #     datasink_connections.append(('modulated_normalized', 'modulated_normalized'))
    #
    # if in_write_deformation_fields is not None:
    #     if in_write_deformation_fields[0]:
    #         connections_segment_to_caps.append(('outputnode.out_inverse_deformation_field', 'inverse_deformation_field'))
    #         datasink_connections.append(('inverse_deformation_field', 'inverse_deformation_field'))
    #     if in_write_deformation_fields[1]:
    #         connections_segment_to_caps.append(('outputnode.out_forward_deformation_field', 'forward_deformation_field'))
    #         datasink_connections.append(('forward_deformation_field', 'forward_deformation_field'))

    datasink_infields = ['native_space', 'dartel_input']
    datasink_connections = [(('new_segment.native_class_images', group_images_by_subject), 'native_space'), (('new_segment.dartel_input_images', group_images_by_subject), 'dartel_input')]

    if save_warped_unmodulated:
        datasink_connections.append((('new_segment.normalized_class_images', group_images_by_subject), 'normalized'))
        datasink_infields.append('normalized')

    if save_warped_modulated:
        datasink_connections.append((('new_segment.modulated_class_images', group_images_by_subject), 'modulated_normalized'))
        datasink_infields.append('modulated_normalized')

    if in_write_deformation_fields is not None:
        if in_write_deformation_fields[0]:
            datasink_connections.append(('new_segment.inverse_deformation_field', 'inverse_deformation_field'))
        datasink_infields.append('inverse_deformation_field')
        if in_write_deformation_fields[1]:
            datasink_connections.append(('new_segment.forward_deformation_field', 'forward_deformation_field'))
        datasink_infields.append('forward_deformation_field')

    datasink_iterfields = ['substitutions'] + datasink_infields

    datasink = pe.MapNode(nio.DataSink(infields=datasink_infields), name='datasink', iterfield=datasink_iterfields)
    datasink.inputs.parameterization = False
    datasink.inputs.base_directory = output_directory
    datasink.inputs.container = 'analysis-series-' + analysis_series_id + '/subjects/sub-SUBJECT_ID/ses-SESSION_ID/t1/spm/segmentation'
    datasink.inputs.substitutions = [[('SUBJECT_ID', subjects[i]), ('SESSION_ID', sessions[i])] for i in range(len(subjects))]

    seg_wf = pe.Workflow(name='seg_wf')
    seg_wf.base_dir = working_directory
    seg_wf.connect([
        (selectfiles, seg_prep_wf, [('out_files', 'new_segment.channel_files')]),
        (seg_prep_wf, datasink, datasink_connections)
        # (seg_prep_wf, write_segment_to_caps, connections_segment_to_caps),
        # (write_segment_to_caps, datasink, datasink_connections)
    ])

    return seg_wf


def datagrabber_t1_spm_full_pipeline_bids(input_directory,
                                     output_directory,
                                     subjects_visits_tsv,
                                     analysis_series_id,
                                     group_id,
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

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(input_directory):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove '.nii'
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = input_directory
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
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

    preproc_wf = pe.Workflow(name='preproc_wf')
    if working_directory is not None:
        preproc_wf.base_dir = working_directory

    preproc_wf.connect([
        (selectfiles, t1_spm_prep_wf, [('out_files', 'segmentation_wf.new_segment.channel_files')])
    ])

    return preproc_wf

