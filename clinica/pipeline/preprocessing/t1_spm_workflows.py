import os
import os.path as op
import nipype.interfaces.spm as spm
import nipype.interfaces.matlab as mlab
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from nipype.interfaces.utility import Function
from t1_spm_utils import get_class_images


def segmentation_pipeline(experiment_dir, datasink_directory, name='segmentation_wf'):
    """
    Creates a pipeline that performs SPM unified segmentation.
    It takes a series of MRI T1 images and extracts gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images.

    new_segment node inputs:
    [Mandatory]
    new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented.

    new_segment node outputs:
    new_segment.dartel_input_images: (a list of items which are a list of items which are an existing file name)
        dartel imported class images
    new_segment.native_class_images: (a list of items which are a list of items which are an existing file name)
        native space probability maps

    For more optional new_segment parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment

    :param experiment_dir: Directory to run the workflow.
    :param datasink_directory: Directory to save the resulting images of the segmentation process.
    :param name: Workflow name
    :return: Segmentation workflow
    """
    spm_home = os.getenv("SPM_HOME")
    mlab_home = os.getenv("MATLABCMD")
    mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
    mlab.MatlabCommand.set_default_paths(spm_home)
    
    version = spm.Info.version()
    print 'VERSION DE SPM'
    print version
    
    new_segment = pe.MapNode(spm.NewSegment(), name='new_segment', iterfield=['channel_files'])

    # Tissue images ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated))

    if version:
        spm_path = version['path']
        if version['name'] == 'SPM8':
            tissue_map = op.join(spm_path,'toolbox/Seg/TPM.nii')
        elif version['name'] == 'SPM12':
            tissue_map = op.join(spm_path,'tpm/TPM.nii')

    tissue1 = ((tissue_map, 1), 2, (True, True), (False, False))
    tissue2 = ((tissue_map, 2), 2, (True, False), (False, False))
    tissue3 = ((tissue_map, 3), 2, (True, False), (False, False))
    tissue4 = ((tissue_map, 4), 3, (False, False), (False, False))
    tissue5 = ((tissue_map, 5), 4, (False, False), (False, False))
    tissue6 = ((tissue_map, 6), 2, (False, False), (False, False))

    new_segment.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.parameterization = False
    datasink.inputs.base_directory = op.join(datasink_directory, 'segmentation')

    wf = pe.Workflow(name=name)
    wf.base_dir = experiment_dir
    wf.connect([
        (new_segment, datasink, [(('native_class_images', get_class_images, [1, 2, 3]), 'native_space'),
                                  (('dartel_input_images', get_class_images, [1, 2, 3]), 'dartel_input')])
    ])
    return wf


def dartel_pipeline(experiment_dir, datasink_directory, name='dartel_wf', modulate=True, smooth=0):
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

    For more dartelTemplate parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartel

    For more dartel2mni parameters and outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#dartelnorm2mni

    :param experiment_dir: Directory to run the workflow.
    :param datasink_directory: Directory to save the resulting images of the registration process.
    :param name: Workflow name
    :param modulate: (a boolean) Modulate out images - no modulation preserves concentrations
    :param smooth: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :return: Registration workflow
    """

    # DARTEL Template creation node
    dartelTemplate = pe.Node(spm.DARTEL(), name='dartelTemplate')

    def prepare_dartel2mni_input(native_class_images, flowfield_files):
        """
        Function receives a list of lists of native class images and a list of flow fields.
        It returns two lists of the same length to give as input to DARTEL2MNI
        This will allow to process each pair of files in parallel.
        :param native_class_images: list of lists of native class images
        :param flowfield_files: list of flow fields files
        :return: expanded list of native images,list of flow fields files of the same length
        """

        native_files = [image for nat_class in native_class_images for image in nat_class]

        if len(native_files)%len(flowfield_files) != 0:
            raise ValueError('Length of the list of native class images is not a multiple of the length of the list of flow fields images')

        ffield_files = flowfield_files * (len(native_files)/len(flowfield_files))

        return native_files, ffield_files

    dartel2mni_input = pe.Node(Function(input_names=['native_class_images', 'flowfield_files'],
                         output_names=['native_files', 'ffield_files'],
                         function=prepare_dartel2mni_input),
                name='dartel2mni_input')

    # DARTEL2MNI
    dartel2mni = pe.MapNode(spm.DARTELNorm2MNI(), name='dartel2MNI', iterfield=['apply_to_files', 'flowfield_files'])
    #Modulation
    dartel2mni.inputs.modulate = modulate
    #Smoothing
    dartel2mni.inputs.fwhm = smooth

    datasink_dartel = pe.Node(nio.DataSink(), name='datasink_dartel')
    datasink_dartel.inputs.parameterization = False
    datasink_dartel.inputs.base_directory = op.join(datasink_directory, 'dartel')

    wf = pe.Workflow(name=name)
    wf.base_dir = experiment_dir

    wf.connect([(dartelTemplate, dartel2mni_input,[('dartel_flow_fields', 'flowfield_files')]),
                        (dartelTemplate, dartel2mni, [('final_template_file', 'template_file')]),
                        (dartel2mni_input, dartel2mni, [('native_files', 'apply_to_files'),
                                                        ('ffield_files', 'flowfield_files')]),
                        (dartelTemplate, datasink_dartel, [('dartel_flow_fields', 'flow_fields'),
                                                           ('final_template_file', 'template')]),
                        (dartel2mni, datasink_dartel, [('normalized_files', 'final')])
                        ])
    return wf


def t1_spm_prep_pipeline(experiment_dir, datasink_directory, name='T1_SPM_prep_wf', modulate=True, smooth=0):
    """
    Creates a pipeline that performs SPM preprocessing for T1 images.
    First unified segmentation is done and then a DARTEL registration.

    It takes a series of MRI T1 images and returns gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images registered
    into a DARTEL template obtained for the group.




    :param experiment_dir: Directory to run the workflow.
    :param datasink_directory: Directory to save the resulting images of segmentation and registration processes.
    :param name: Workflow name
    :param modulate: (a boolean) Modulate out images - no modulation preserves concentrations
    :param smooth: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
    :return: Segmentation and registration workflow
    """

    """
    Outputs:
    new_segment.dartel_input_images: (a list of items which are a list of items which are an existing file name)
        dartel imported class images
    new_segment.native_class_images: (a list of items which are a list of items which are an existing file name)
        native space probability maps
    For more new_segment outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment
    """

    segmentation_wf = segmentation_pipeline(experiment_dir, datasink_directory)
    dartel_wf = dartel_pipeline(experiment_dir, datasink_directory)

    wf = pe.Workflow(name=name)
    wf.base_dir = experiment_dir
    wf.connect([
        (segmentation_wf, dartel_wf, [(('new_segment.native_class_images', get_class_images, [1, 2, 3]), 'dartel2mni_input.native_class_images'),
                                  (('new_segment.dartel_input_images', get_class_images, [1, 2, 3]), 'dartelTemplate.image_files')])
    ])

    return wf
