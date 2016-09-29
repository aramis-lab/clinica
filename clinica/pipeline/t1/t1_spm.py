

from os import walk
from clinica.pipeline.t1.t1_spm_workflows import segmentation_pipeline, dartel_pipeline, t1_spm_full_pipeline
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio


def datagrabber_t1_spm_segment_pipeline(input_directory,
                                        output_directory,
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

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(input_directory):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove .nii
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = input_directory
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    # Creating T1 preprocessing workflow
    seg_prep_wf = segmentation_pipeline(output_directory,
                                        working_directory=working_directory,
                                        tissue_classes=tissue_classes,
                                        dartel_tissues=dartel_tissues,
                                        save_warped_unmodulated=save_warped_unmodulated,
                                        save_warped_modulated=save_warped_modulated,
                                        in_write_deformation_fields=in_write_deformation_fields)

    seg_wf = pe.Workflow(name='seg_wf')
    seg_wf.base_dir = working_directory
    seg_wf.connect([
        (selectfiles, seg_prep_wf, [('out_files', 'new_segment.channel_files')])
    ])

    return seg_wf


def datagrabber_t1_spm_full_pipeline(input_directory,
                                     output_directory,
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

    t1_spm_prep_wf = t1_spm_full_pipeline(output_directory,
                                          working_directory=working_directory,
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

