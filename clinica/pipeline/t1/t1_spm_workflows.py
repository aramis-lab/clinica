import os
import os.path as op
import nipype.interfaces.spm as spm
import nipype.interfaces.matlab as mlab
import nipype.pipeline.engine as npe
import nipype.interfaces.utility as niu
import nipype.interfaces.io as nio
from clinica.pipeline.t1.t1_spm_utils import get_tissue_tuples, get_class_images, DARTELExistingTemplate
import clinica.pipeline.engine as cpe
from clinica.pipeline.t1.t1_spm_utils import select_image
from clinica.utils.io import unzip_nii, zip_nii


def segmentation_pipeline(working_directory=None,
                          name='segmentation_wf',
                          tissue_classes=[1, 2, 3],
                          dartel_tissues=[1, 2, 3],
                          save_warped_unmodulated=False,
                          save_warped_modulated=False,
                          affine_regularization=None,
                          channel_info=None,
                          sampling_distance=None,
                          warping_regularization=None,
                          write_deformation_fields=None):
    """
    Creates a pipeline that performs SPM unified segmentation.
    It takes a series of MRI T1 images and extracts gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images.

    new_segment node inputs:

    [Mandatory]
    new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented.


    outputnode node outputs:

    bias_corrected_images: (a list of items which are an existing file name)
        bias corrected images
    bias_field_images: (a list of items which are an existing file name)
            bias field images
    dartel_input_images: (a list of items which are a list of items which are an existing file name)
            dartel imported class images
    forward_deformation_field: (a list of items which are an existing file name)
    inverse_deformation_field: (a list of items which are an existing file name)
    modulated_class_images: (a list of items which are a list of items which are an existing file name)
            modulated+normalized class images
    native_tissue_classes: (a list of items which are a list of items which are an existing file name)
            native space probability maps
    normalized_class_images: (a list of items which are a list of items which are an existing file name)
            normalized class images
    transformation_mat: (a list of items which are an existing file name)
            Normalization transformation

    For more optional new_segment parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment

    :param output_directory: Directory to save the resulting images of the segmentation process.
    :param working_directory: Directory to run the workflow.
    :param name: Workflow name
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :param affine_regularization: ('mni' or 'eastern' or 'subj' or 'none')
    :param channel_info: a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
            A tuple with the following fields:
             - bias reguralisation (0-10)
             - FWHM of Gaussian smoothness of bias
             - which maps to save (Corrected, Field) - a tuple of two boolean values
    :param sampling_distance: (a float) Sampling distance on data for parameter estimation
    :param warping_regularization: (a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)
    :param write_deformation_fields: Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]

    :return: Segmentation workflow
    """

    spm_home = os.getenv("SPM_HOME")
    mlab_home = os.getenv("MATLABCMD")
    mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
    mlab.MatlabCommand.set_default_paths(spm_home)
    
    version = spm.Info.version()

    if version:
        spm_path = version['path']
        if version['name'] == 'SPM8':
            print 'You are using SPM version 8. The recommended version to use with Clinica is SPM 12. Please upgrade your SPM toolbox.'
            tissue_map = op.join(spm_path,'toolbox/Seg/TPM.nii')
        elif version['name'] == 'SPM12':
            tissue_map = op.join(spm_path,'tpm/TPM.nii')
        else:
            raise RuntimeError('SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
    else:
        raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')


    unzip = npe.MapNode(niu.Function(input_names=['in_file'],
                                            output_names=['out_file'],
                                            function=unzip_nii),
                        name='unzip', iterfield=['in_file'])

    new_segment = npe.MapNode(spm.NewSegment(), name='new_segment', iterfield=['channel_files'])

    if affine_regularization is not None:
        new_segment.inputs.affine_regularization = affine_regularization
    if channel_info is not None:
        new_segment.inputs.channel_info = channel_info
    if sampling_distance is not None:
        new_segment.inputs.sampling_distance = sampling_distance
    if warping_regularization is not None:
        new_segment.inputs.warping_regularization = warping_regularization
    if write_deformation_fields is not None:
        new_segment.inputs.write_deformation_fields = write_deformation_fields

    new_segment.inputs.tissues = get_tissue_tuples(tissue_map, tissue_classes, dartel_tissues, save_warped_unmodulated, save_warped_modulated)

    outputnode = npe.Node(niu.IdentityInterface(fields=['bias_corrected_images', 'bias_field_images',
                                                       'dartel_input_images', 'forward_deformation_field',
                                                       'inverse_deformation_field', 'modulated_class_images',
                                                       'native_class_images', 'normalized_class_images',
                                                       'transformation_mat']),
                         name='outputnode')

    wf = npe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory


    wf.connect([
        (unzip, new_segment, [('out_file', 'channel_files')]),
        (new_segment, outputnode, [('bias_corrected_images','bias_corrected_images'),
                                   ('bias_field_images', 'bias_field_images'),
                                   (('dartel_input_images', get_class_images, tissue_classes), 'dartel_input_images'),
                                   ('forward_deformation_field', 'forward_deformation_field'),
                                   ('inverse_deformation_field', 'inverse_deformation_field'),
                                   (('modulated_class_images', get_class_images, tissue_classes), 'modulated_class_images'),
                                   (('native_class_images', get_class_images, tissue_classes), 'native_class_images'),
                                   (('normalized_class_images', get_class_images, tissue_classes), 'normalized_class_images'),
                                   ('transformation_mat', 'transformation_mat')])
    ])

    return wf


def dartel_pipeline(working_directory=None,
                    name='dartel_wf',
                    iteration_parameters=None,
                    optimization_parameters=None,
                    regularization_form=None,
                    template_prefix=None,
                    bounding_box=None,
                    fwhm=None,
                    modulate=True,
                    voxel_size=None):
    """
    Creates a pipeline that performs SPM DARTEL registration for T1 images.

    It takes a series of DARTEL input images and gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images in native space
    and returns the registered images into a DARTEL template obtained for the group.

    dartelTemplate node inputs:

    [Mandatory]
    dartelTemplate.image_files: (a list of items which are a list of items which are an existing file name)
        A list of files to be registered

    dartel2mni_input node inputs:

    [Mandatory]
    dartel2mni_input.native_class_images: (a list of items which are an existing file name)
        Files to apply the transform to


    outputnode node outputs:

    dartel_flow_fields: (a list of items which are an existing file name)
        DARTEL flow fields
    final_template_file: (an existing file name)
            final DARTEL template
    template_files: (a list of items which are an existing file name)
            Templates from different stages of iteration
    normalization_parameter_file: (an existing file name)
        Transform parameters to MNI space
    normalized_files: (a list of items which are an existing file name)
        Normalized files in MNI space



    For more dartelTemplate parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartel

    For more dartel2mni parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartelnorm2mni

    :param output_directory: Directory to save the resulting images of the registration process.
    :param working_directory: Directory to run the workflow.
    :param name: Workflow name

    DARTEL parameters
    :param iteration_parameters: (a list of from 3 to 12 items which are a tuple
         of the form: (1 <= an integer <= 10, a tuple of the form: (a float,
         a float, a float), 1 or 2 or 4 or 8 or 16 or 32 or 64 or 128 or 256
         or 512, 0 or 0.5 or 1 or 2 or 4 or 8 or 16 or 32))
        List of tuples for each iteration
         - Inner iterations
         - Regularization parameters
         - Time points for deformation model
         - smoothing parameter
    :param optimization_parameters: (a tuple of the form: (a float, 1 <= an
         integer <= 8, 1 <= an integer <= 8))
         Optimization settings a tuple
         - LM regularization
         - cycles of multigrid solver
         - relaxation iterations
    :param regularization_form: ('Linear' or 'Membrane' or 'Bending')
        Form of regularization energy term
    :param template_prefix: (a string, nipype default value: Template)
        Prefix for template

    DARTELNorm2MNI parameters
    :param bounding_box: (a tuple of the form: (a float, a float, a float, a
         float, a float, a float))
        Voxel sizes for output file
    :param fwhm: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :param modulate: (a boolean)
        Modulate out images - no modulation preserves concentrations
    :param voxel_size: (a tuple of the form: (a float, a float, a float))
        Voxel sizes for output file
    :return: Registration workflow
    """
    spm_home = os.getenv("SPM_HOME")
    mlab_home = os.getenv("MATLABCMD")
    mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
    mlab.MatlabCommand.set_default_paths(spm_home)

    version = spm.Info.version()

    if version:
        if version['name'] == 'SPM8':
            print 'You are using SPM version 8. The recommended version to use with Clinica is SPM 12. Please upgrade your SPM toolbox.'
        elif version['name'] != 'SPM12':
            raise RuntimeError('SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
    else:
        raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')

    # DARTEL Template creation node
    dartelTemplate = npe.Node(spm.DARTEL(), name='dartelTemplate')

    if iteration_parameters is not None:
        dartelTemplate.inputs.iteration_parameters = iteration_parameters
    if optimization_parameters is not None:
        dartelTemplate.inputs.optimization_parameters = optimization_parameters
    if regularization_form is not None:
        dartelTemplate.inputs.regularization_form = regularization_form
    if template_prefix is not None:
        dartelTemplate.inputs.template_prefix = template_prefix


    def prepare_dartel2mni_input(native_space_images, flowfield_files):
        """
        Function receives a list of lists of native class images and a list of flow fields.
        It returns two lists of the same length to give as input to DARTEL2MNI
        This will allow to process each pair of files in parallel.
        :param native_space_images: list of lists of native space images
        :param flowfield_files: list of flow fields files
        :return: expanded list of native images,list of flow fields files of the same length
        """

        native_files = [image for nat_class in native_space_images for image in nat_class]

        if len(native_files)%len(flowfield_files) != 0:
            raise ValueError('Length of the list of native space images is not a multiple of the length of the list of flow fields images')

        ffield_files = flowfield_files * int(len(native_files)/len(flowfield_files))

        return native_files, ffield_files

    dartel2mni_input = npe.Node(niu.Function(input_names=['native_space_images', 'flowfield_files'],
                         output_names=['native_files', 'ffield_files'],
                         function=prepare_dartel2mni_input),
                name='dartel2mni_input')

    # DARTEL2MNI
    dartel2mni = npe.MapNode(spm.DARTELNorm2MNI(), name='dartel2MNI', iterfield=['apply_to_files', 'flowfield_files'])

    if bounding_box is not None:
        dartel2mni.inputs.bounding_box = bounding_box
    if voxel_size is not None:
        dartel2mni.inputs.voxel_size = voxel_size

    #Modulation
    dartel2mni.inputs.modulate = modulate

    #Smoothing
    if fwhm is not None:
        if fwhm == [0, 0, 0]:
            fwhm = 0
        dartel2mni.inputs.fwhm = fwhm

    outputnode = npe.Node(niu.IdentityInterface(fields=['dartel_flow_fields', 'final_template_file', 'template_files',
                                                       'normalization_parameter_file', 'normalized_files']), name='outputnode')

    wf = npe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory

    wf.connect([(dartelTemplate, dartel2mni_input,[('dartel_flow_fields', 'flowfield_files')]),
                (dartelTemplate, dartel2mni, [('final_template_file', 'template_file')]),
                (dartel2mni_input, dartel2mni, [('native_files', 'apply_to_files'),
                                                ('ffield_files', 'flowfield_files')]),
                (dartelTemplate, outputnode, [('dartel_flow_fields','dartel_flow_fields'),
                                              ('final_template_file', 'final_template_file'),
                                              ('template_files', 'template_files')]),
                (dartel2mni, outputnode, [('normalization_parameter_file', 'normalization_parameter_file'),
                                          ('normalized_files', 'normalized_files')])
                ])
    return wf


def dartel_existing_template_pipeline(working_directory=None,
                                      name='dartel_existing_template_wf',
                                      iteration_parameters=None,
                                      optimization_parameters=None,
                                      regularization_form=None,
                                      template_prefix=None,
                                      bounding_box=None,
                                      fwhm=None,
                                      modulate=True,
                                      voxel_size=None):
    """
    Creates a pipeline that performs SPM DARTEL registration for T1 images.

    It takes a series of DARTEL input images and gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images in native space
    and returns the registered images into a DARTEL template obtained for the group.

    dartelTemplate node inputs:

    [Mandatory]
    dartelTemplate.image_files: (a list of items which are a list of items which are an existing file name)
        A list of files to be registered

    dartel2mni_input node inputs:

    [Mandatory]
    dartel2mni_input.native_class_images: (a list of items which are an existing file name)
        Files to apply the transform to


    outputnode node outputs:

    dartel_flow_fields: (a list of items which are an existing file name)
        DARTEL flow fields
    final_template_file: (an existing file name)
            final DARTEL template
    template_files: (a list of items which are an existing file name)
            Templates from different stages of iteration
    normalization_parameter_file: (an existing file name)
        Transform parameters to MNI space
    normalized_files: (a list of items which are an existing file name)
        Normalized files in MNI space



    For more dartelTemplate parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartel

    For more dartel2mni parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartelnorm2mni

    :param output_directory: Directory to save the resulting images of the registration process.
    :param working_directory: Directory to run the workflow.
    :param name: Workflow name

    DARTEL parameters
    :param iteration_parameters: (a list of from 3 to 12 items which are a tuple
         of the form: (1 <= an integer <= 10, a tuple of the form: (a float,
         a float, a float), 1 or 2 or 4 or 8 or 16 or 32 or 64 or 128 or 256
         or 512, 0 or 0.5 or 1 or 2 or 4 or 8 or 16 or 32))
        List of tuples for each iteration
         - Inner iterations
         - Regularization parameters
         - Time points for deformation model
         - smoothing parameter
    :param optimization_parameters: (a tuple of the form: (a float, 1 <= an
         integer <= 8, 1 <= an integer <= 8))
         Optimization settings a tuple
         - LM regularization
         - cycles of multigrid solver
         - relaxation iterations
    :param regularization_form: ('Linear' or 'Membrane' or 'Bending')
        Form of regularization energy term
    :param template_prefix: (a string, nipype default value: Template)
        Prefix for template

    DARTELNorm2MNI parameters
    :param bounding_box: (a tuple of the form: (a float, a float, a float, a
         float, a float, a float))
        Voxel sizes for output file
    :param fwhm: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :param modulate: (a boolean)
        Modulate out images - no modulation preserves concentrations
    :param voxel_size: (a tuple of the form: (a float, a float, a float))
        Voxel sizes for output file
    :return: Registration workflow
    """
    spm_home = os.getenv("SPM_HOME")
    mlab_home = os.getenv("MATLABCMD")
    mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
    mlab.MatlabCommand.set_default_paths(spm_home)

    version = spm.Info.version()

    if version:
        if version['name'] == 'SPM8':
            print 'You are using SPM version 8. The recommended version to use with Clinica is SPM 12. Please upgrade your SPM toolbox.'
        elif version['name'] != 'SPM12':
            raise RuntimeError('SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
    else:
        raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')

    # DARTEL Template creation node
    dartelTemplate = npe.Node(spm.DARTEL(), name='dartelTemplate')

    if iteration_parameters is not None:
        dartelTemplate.inputs.iteration_parameters = iteration_parameters
    if optimization_parameters is not None:
        dartelTemplate.inputs.optimization_parameters = optimization_parameters
    if regularization_form is not None:
        dartelTemplate.inputs.regularization_form = regularization_form
    if template_prefix is not None:
        dartelTemplate.inputs.template_prefix = template_prefix


    def prepare_dartel2mni_input(native_space_images, flowfield_files):
        """
        Function receives a list of lists of native class images and a list of flow fields.
        It returns two lists of the same length to give as input to DARTEL2MNI
        This will allow to process each pair of files in parallel.
        :param native_space_images: list of lists of native space images
        :param flowfield_files: list of flow fields files
        :return: expanded list of native images,list of flow fields files of the same length
        """

        native_files = [image for nat_class in native_space_images for image in nat_class]

        if len(native_files)%len(flowfield_files) != 0:
            raise ValueError('Length of the list of native space images is not a multiple of the length of the list of flow fields images')

        ffield_files = flowfield_files * (len(native_files)/len(flowfield_files))

        return native_files, ffield_files

    dartel2mni_input = npe.Node(niu.Function(input_names=['native_space_images', 'flowfield_files'],
                         output_names=['native_files', 'ffield_files'],
                         function=prepare_dartel2mni_input),
                name='dartel2mni_input')

    # DARTEL2MNI
    dartel2mni = npe.MapNode(spm.DARTELNorm2MNI(), name='dartel2MNI', iterfield=['apply_to_files', 'flowfield_files'])

    if bounding_box is not None:
        dartel2mni.inputs.bounding_box = bounding_box
    if voxel_size is not None:
        dartel2mni.inputs.voxel_size = voxel_size

    #Modulation
    dartel2mni.inputs.modulate = modulate

    #Smoothing
    if fwhm is not None:
        if fwhm == [0, 0, 0]:
            fwhm = 0
        dartel2mni.inputs.fwhm = fwhm

    outputnode = npe.Node(niu.IdentityInterface(fields=['dartel_flow_fields', 'final_template_file', 'template_files',
                                                       'normalization_parameter_file', 'normalized_files']), name='outputnode')

    wf = npe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory

    wf.connect([(dartelTemplate, dartel2mni_input,[('dartel_flow_fields', 'flowfield_files')]),
                (dartelTemplate, dartel2mni, [('final_template_file', 'template_file')]),
                (dartel2mni_input, dartel2mni, [('native_files', 'apply_to_files'),
                                                ('ffield_files', 'flowfield_files')]),
                (dartelTemplate, outputnode, [('dartel_flow_fields','dartel_flow_fields'),
                                              ('final_template_file', 'final_template_file'),
                                              ('template_files', 'template_files')]),
                (dartel2mni, outputnode, [('normalization_parameter_file', 'normalization_parameter_file'),
                                          ('normalized_files', 'normalized_files')])
                ])
    return wf


def t1_spm_full_pipeline(working_directory=None,
                         name='segmentation_wf',
                         tissue_classes=[1, 2, 3],
                         dartel_tissues=[1],
                         save_warped_unmodulated=False,
                         save_warped_modulated=False,
                         affine_regularization=None,
                         channel_info=None,
                         sampling_distance=None,
                         warping_regularization=None,
                         write_deformation_fields=None,
                         iteration_parameters=None,
                         optimization_parameters=None,
                         regularization_form=None,
                         template_prefix=None,
                         bounding_box=None,
                         fwhm=0,
                         modulate=True,
                         voxel_size=None):
    """
    Creates a pipeline that performs SPM preprocessing for T1 images.
    First unified segmentation is done and then a DARTEL registration.

    It takes a series of MRI T1 images and returns gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images registered
    into a DARTEL template obtained for the group.


    new_segment node inputs:

    [Mandatory]
    segmentation_wf.new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented and registered.


    outputnode node outputs:

    bias_corrected_images: (a list of items which are an existing file name)
        bias corrected images
    bias_field_images: (a list of items which are an existing file name)
            bias field images
    dartel_input_images: (a list of items which are a list of items which are an existing file name)
            dartel imported class images
    forward_deformation_field: (a list of items which are an existing file name)
    inverse_deformation_field: (a list of items which are an existing file name)
    modulated_class_images: (a list of items which are a list of items which are an existing file name)
            modulated+normalized class images
    native_class_images: (a list of items which are a list of items which are an existing file name)
            native space probability maps
    normalized_class_images: (a list of items which are a list of items which are an existing file name)
            normalized class images
    transformation_mat: (a list of items which are an existing file name)
            Normalization transformation
    dartel_flow_fields: (a list of items which are an existing file name)
        DARTEL flow fields
    final_template_file: (an existing file name)
            final DARTEL template
    template_files: (a list of items which are an existing file name)
            Templates from different stages of iteration
    normalization_parameter_file: (an existing file name)
        Transform parameters to MNI space
    normalized_files: (a list of items which are an existing file name)
        Normalized files in MNI space


    For more NewSegment, DARTEL and DARTELNorm2MNI outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartel
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartelnorm2mni

    :param output_directory: Directory to save the resulting images of segmentation and registration processes.
    :param working_directory: Directory to run the workflow.
    :param name: Workflow name
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'

    NewSegment parameters
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :param affine_regularization: ('mni' or 'eastern' or 'subj' or 'none')
    :param channel_info: a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
            A tuple with the following fields:
             - bias reguralisation (0-10)
             - FWHM of Gaussian smoothness of bias
             - which maps to save (Corrected, Field) - a tuple of two boolean values
    :param sampling_distance: (a float) Sampling distance on data for parameter estimation
    :param warping_regularization: (a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)
    :param write_deformation_fields: Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]


    DARTEL parameters
    :param iteration_parameters: (a list of from 3 to 12 items which are a tuple
         of the form: (1 <= an integer <= 10, a tuple of the form: (a float,
         a float, a float), 1 or 2 or 4 or 8 or 16 or 32 or 64 or 128 or 256
         or 512, 0 or 0.5 or 1 or 2 or 4 or 8 or 16 or 32))
        List of tuples for each iteration
         - Inner iterations
         - Regularization parameters
         - Time points for deformation model
         - smoothing parameter
    :param optimization_parameters: (a tuple of the form: (a float, 1 <= an
         integer <= 8, 1 <= an integer <= 8))
         Optimization settings a tuple
         - LM regularization
         - cycles of multigrid solver
         - relaxation iterations
    :param regularization_form: ('Linear' or 'Membrane' or 'Bending')
        Form of regularization energy term
    :param template_prefix: (a string, nipype default value: Template)
        Prefix for template

    DARTELNorm2MNI parameters
    :param bounding_box: (a tuple of the form: (a float, a float, a float, a
         float, a float, a float))
        Voxel sizes for output file
    :param fwhm: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :param modulate: (a boolean)
        Modulate out images - no modulation preserves concentrations
    :param voxel_size: (a tuple of the form: (a float, a float, a float))
        Voxel sizes for output file
    :return: Segmentation and registration workflow
    """

    segmentation_wf = segmentation_pipeline(working_directory=working_directory,
                                            name='segmentation_wf',
                                            tissue_classes=tissue_classes,
                                            dartel_tissues=dartel_tissues,
                                            save_warped_unmodulated=save_warped_unmodulated,
                                            save_warped_modulated=save_warped_modulated,
                                            affine_regularization=affine_regularization,
                                            channel_info=channel_info,
                                            sampling_distance=sampling_distance,
                                            warping_regularization=warping_regularization,
                                            write_deformation_fields=write_deformation_fields)

    dartel_wf = dartel_pipeline(working_directory=working_directory,
                                name='dartel_wf',
                                iteration_parameters=iteration_parameters,
                                optimization_parameters=optimization_parameters,
                                regularization_form=regularization_form,
                                template_prefix=template_prefix,
                                bounding_box=bounding_box,
                                fwhm=fwhm,
                                modulate=modulate,
                                voxel_size=voxel_size)

    outputnode = npe.Node(niu.IdentityInterface(fields=['bias_corrected_images', 'bias_field_images',
                                                       'dartel_input_images', 'forward_deformation_field',
                                                       'inverse_deformation_field', 'modulated_class_images',
                                                       'native_class_images', 'normalized_class_images',
                                                       'transformation_mat', 'dartel_flow_fields',
                                                       'final_template_file', 'template_files',
                                                       'normalization_parameter_file', 'normalized_files']),
                         name='outputnode')

    wf = npe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory

    wf.connect([
        (segmentation_wf, dartel_wf, [(('new_segment.native_class_images', get_class_images, tissue_classes), 'dartel2mni_input.native_space_images'),
                                      (('new_segment.dartel_input_images', get_class_images, dartel_tissues), 'dartelTemplate.image_files')]),
        (segmentation_wf, outputnode, [('outputnode.bias_corrected_images', 'bias_corrected_images'),
                                       ('outputnode.bias_field_images', 'bias_field_images'),
                                       ('outputnode.dartel_input_images', 'dartel_input_images'),
                                       ('outputnode.forward_deformation_field', 'forward_deformation_field'),
                                       ('outputnode.inverse_deformation_field', 'inverse_deformation_field'),
                                       ('outputnode.modulated_class_images', 'modulated_class_images'),
                                       ('outputnode.native_class_images', 'native_class_images'),
                                       ('outputnode.normalized_class_images', 'normalized_class_images'),
                                       ('outputnode.transformation_mat', 'transformation_mat')]),
        (dartel_wf, outputnode, [('outputnode.dartel_flow_fields','dartel_flow_fields'),
                                 ('outputnode.final_template_file', 'final_template_file'),
                                 ('outputnode.template_files', 'template_files'),
                                 ('outputnode.normalization_parameter_file', 'normalization_parameter_file'),
                                 ('outputnode.normalized_files', 'normalized_files')])
    ])

    return wf
