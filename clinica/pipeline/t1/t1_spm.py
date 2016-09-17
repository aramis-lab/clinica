

from os import walk
from clinica.pipeline.t1.t1_spm_workflows import segmentation_pipeline, dartel_pipeline, t1_spm_full_pipeline
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio


def datagrabber_t1_spm_segment_pipeline(data_dir, experiment_dir, datasink_directory, class_images=[1, 2, 3],
                                        in_affine_regularization=None, in_channel_info=None, in_sampling_distance=None,
                                        in_tissues_to_save=None, in_warping_regularization=None,
                                        in_write_deformation_fields=None):
    """

    T1 segmentation workflow with DataGrabber as input.


    :param data_dir: Directory containing the NIFTI files
    :param experiment_dir: Directory to run the workflow
    :param datasink_directory: Directory to save the resulting images of segmentation and registration processes
    :param class_images: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF

    NewSegment parameters
    :param in_affine_regularization: ('mni' or 'eastern' or 'subj' or 'none')
    :param in_channel_info: a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
        A tuple with the following fields:
         - bias reguralisation (0-10)
         - FWHM of Gaussian smoothness of bias
         - which maps to save (Corrected, Field) - a tuple of two boolean values
    :param in_sampling_distance: (a float) Sampling distance on data for parameter estimation
    :param in_tissues_to_save: A list of tuples (one per tissue) of the form:
        ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated)) with the boolen value for the type of images to save
        Ex.: [((True, True), (False, False)),
              ((True, False), (False, False))]
    :param in_warping_regularization: (a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)
    :param in_write_deformation_fields: (a list of from 2 to 2 items which are a boolean)
        Which deformation fields to write:[Inverse, Forward]

    :return: SPM Segmentation workflow
    """

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove .nii
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = data_dir
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    # Creating T1 preprocessing workflow
    seg_prep_wf = segmentation_pipeline(experiment_dir, datasink_directory, class_images=class_images,
                                        in_affine_regularization=in_affine_regularization, in_channel_info=in_channel_info,
                                        in_sampling_distance=in_sampling_distance, in_tissues_to_save=in_tissues_to_save,
                                        in_warping_regularization=in_warping_regularization,
                                        in_write_deformation_fields=in_write_deformation_fields)

    seg_wf = pe.Workflow(name='seg_wf')
    seg_wf.base_dir = experiment_dir
    seg_wf.connect([
        (selectfiles, seg_prep_wf, [('out_files', 'new_segment.channel_files')])
    ])

    return seg_wf


def datagrabber_t1_spm_full_pipeline(data_dir, experiment_dir, datasink_directory, class_images=[1, 2, 3],
                         dartel_class_images=[1], in_affine_regularization=None, in_channel_info=None,
                         in_sampling_distance=None, in_tissues_to_save=None, in_warping_regularization=None,
                         in_write_deformation_fields=None, in_iteration_parameters=None, in_optimization_parameters=None,
                         in_regularization_form=None, in_template_prefix=None, in_bounding_box=None,
                         in_fwhm=0, in_modulate=True, in_voxel_size=None):
    """

    T1 preprocessing workflow with DataGrabber as input.

    :param data_dir: Directory containing the NIFTI files
    :param experiment_dir: Directory to run the workflow
    :param datasink_directory: Directory to save the resulting images of segmentation and registration processes
    :param class_images: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_class_images: Classes of images to use for DARTEL template calculation. Ex: [1] is only GM'

    NewSegment parameters
    :param in_affine_regularization: ('mni' or 'eastern' or 'subj' or 'none')
    :param in_channel_info: a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
        A tuple with the following fields:
         - bias reguralisation (0-10)
         - FWHM of Gaussian smoothness of bias
         - which maps to save (Corrected, Field) - a tuple of two boolean values
    :param in_sampling_distance: (a float) Sampling distance on data for parameter estimation
    :param in_tissues_to_save: A list of tuples (one per tissue) of the form:
        ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated)) with the boolen value for the type of images to save
        Ex.: [((True, True), (False, False)),
              ((True, False), (False, False))]
    :param in_warping_regularization: (a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)
    :param in_write_deformation_fields: (a list of from 2 to 2 items which are a boolean)
        Which deformation fields to write:[Inverse, Forward]

    DARTEL parameters
    :param in_iteration_parameters: (a list of from 3 to 12 items which are a tuple
         of the form: (1 <= an integer <= 10, a tuple of the form: (a float,
         a float, a float), 1 or 2 or 4 or 8 or 16 or 32 or 64 or 128 or 256
         or 512, 0 or 0.5 or 1 or 2 or 4 or 8 or 16 or 32))
        List of tuples for each iteration
         - Inner iterations
         - Regularization parameters
         - Time points for deformation model
         - smoothing parameter
    :param in_optimization_parameters: (a tuple of the form: (a float, 1 <= an
         integer <= 8, 1 <= an integer <= 8))
         Optimization settings a tuple
         - LM regularization
         - cycles of multigrid solver
         - relaxation iterations
    :param in_regularization_form: ('Linear' or 'Membrane' or 'Bending')
        Form of regularization energy term
    :param in_template_prefix: (a string, nipype default value: Template)
        Prefix for template

    DARTELNorm2MNI parameters
    :param in_bounding_box: (a tuple of the form: (a float, a float, a float, a
         float, a float, a float))
        Voxel sizes for output file
    :param in_fwhm: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :param in_modulate: (a boolean)
        Modulate out images - no modulation preserves concentrations
    :param in_voxel_size: (a tuple of the form: (a float, a float, a float))
        Voxel sizes for output file
    :return: SPM Segmentation and registration workflow
    """

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove '.nii'
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = data_dir
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    t1_spm_prep_wf = t1_spm_full_pipeline(experiment_dir, datasink_directory, name='t1_spm_full_wf',
                                          class_images=class_images, dartel_class_images=dartel_class_images,
                                          in_affine_regularization=in_affine_regularization, in_channel_info=in_channel_info,
                                          in_sampling_distance=in_sampling_distance, in_tissues_to_save=in_tissues_to_save,
                                          in_warping_regularization=in_warping_regularization, in_write_deformation_fields=in_write_deformation_fields,
                                          in_iteration_parameters=in_iteration_parameters, in_optimization_parameters=in_optimization_parameters,
                                          in_regularization_form=in_regularization_form, in_template_prefix=in_template_prefix,
                                          in_bounding_box=in_bounding_box, in_fwhm=in_fwhm, in_modulate=in_modulate,
                                          in_voxel_size=in_voxel_size)

    preproc_wf = pe.Workflow(name='preproc_wf')
    preproc_wf.base_dir = experiment_dir

    preproc_wf.connect([
        (selectfiles, t1_spm_prep_wf, [('out_files', 'segmentation_wf.new_segment.channel_files')])
    ])

    return preproc_wf

