import os.path as op
import nipype.interfaces.spm as spm
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from nipype.interfaces.utility import Function

from T1_SPM_utils import get_class_images


def segmentation_pipeline(experiment_dir, datasink_directory, tissue_map, name='segmentation_wf'):
    """
    Creates a pipeline that performs SPM unified segmentation.
    It takes a series of MRI T1 images and extracts gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images.

    Inputs:
    [Mandatory]
    datasink_directory: Directory to save the resulting images of the segmentation process.
    new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented.
    [Optional]
    new_segment.tissues: (a list of items which are a tuple of the form: (a tuple of
         the form: (an existing file name, an integer), an integer, a tuple
         of the form: (a boolean, a boolean), a tuple of the form: (a
         boolean, a boolean)))
        A list of tuples (one per tissue) with the following fields:
         - tissue probability map (4D), 1-based index to frame
         - number of gaussians
         - which maps to save [Native, DARTEL] - a tuple of two boolean
        values
         - which maps to save [Unmodulated, Modulated] - a tuple of two
        boolean values
        By default GM, WM and CSF tissue maps are saved only in native space.
    For more optional new_segment parameters check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment

    Outputs:
    new_segment.dartel_input_images: (a list of items which are a list of items which are an existing file name)
        dartel imported class images
    new_segment.native_class_images: (a list of items which are a list of items which are an existing file name)
        native space probability maps
    For more new_segment outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment
    """

    new_segment = pe.MapNode(spm.NewSegment(), name='new_segment', iterfield=['channel_files'])

    # Tissue images ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated))
    GM = ((True, True), (False, False))
    WM = ((True, False), (False, False))
    CSF = ((True, False), (False, False))

    tissue1 = ((tissue_map, 1), 2, GM[0], GM[1])
    tissue2 = ((tissue_map, 2), 2, WM[0], WM[1])
    tissue3 = ((tissue_map, 3), 2, CSF[0], CSF[1])

    new_segment.inputs.tissues = [tissue1, tissue2, tissue3]

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

    Inputs:
    [Mandatory]
    datasink_directory: Directory to save the resulting images of the segmentation process.
    new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented.
    [Optional]
    new_segment.tissues: (a list of items which are a tuple of the form: (a tuple of
         the form: (an existing file name, an integer), an integer, a tuple
         of the form: (a boolean, a boolean), a tuple of the form: (a
         boolean, a boolean)))
        A list of tuples (one per tissue) with the following fields:
         - tissue probability map (4D), 1-based index to frame
         - number of gaussians
         - which maps to save [Native, DARTEL] - a tuple of two boolean
        values
         - which maps to save [Unmodulated, Modulated] - a tuple of two
        boolean values
        By default GM, WM and CSF tissue maps are saved only in native space.
    For more optional new_segment parameters check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment

    Outputs:
    new_segment.dartel_input_images: (a list of items which are a list of items which are an existing file name)
        dartel imported class images
    new_segment.native_class_images: (a list of items which are a list of items which are an existing file name)
        native space probability maps
    For more new_segment outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment
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


def T1_SPM_prep_pipeline(experiment_dir, datasink_directory, tissue_map, name='T1_SPM_prep_wf', modulate=True, smooth=0):
    """
    Creates a pipeline that performs SPM preprocessing for T1 images.
    First unified segmentation is done and then a DARTEL registration.

    It takes a series of MRI T1 images and returns gray matter(GM),
    white matter(WM) and cerebrospinal fluid(CSF) tissue images registered
    into a DARTEL template obtained for the group.

    Inputs:
    [Mandatory]
    datasink_directory: Directory to save the resulting images of the segmentation process.
    new_segment.channel_files: (a list of items which are an existing file name)
        A list of files to be segmented.
    [Optional]
    new_segment.tissues: (a list of items which are a tuple of the form: (a tuple of
         the form: (an existing file name, an integer), an integer, a tuple
         of the form: (a boolean, a boolean), a tuple of the form: (a
         boolean, a boolean)))
        A list of tuples (one per tissue) with the following fields:
         - tissue probability map (4D), 1-based index to frame
         - number of gaussians
         - which maps to save [Native, DARTEL] - a tuple of two boolean
        values
         - which maps to save [Unmodulated, Modulated] - a tuple of two
        boolean values
        By default GM, WM and CSF tissue maps are saved only in native space.
    For more optional new_segment parameters check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment

    Outputs:
    new_segment.dartel_input_images: (a list of items which are a list of items which are an existing file name)
        dartel imported class images
    new_segment.native_class_images: (a list of items which are a list of items which are an existing file name)
        native space probability maps
    For more new_segment outputs check:
    http://www.mit.edu/~satra/nipype-nightly/interfaces/generated/nipype.interfaces.spm.preprocess.html#newsegment
    """

    segmentation_wf = segmentation_pipeline(experiment_dir, datasink_directory, tissue_map)
    dartel_wf = dartel_pipeline(experiment_dir, datasink_directory)

    wf = pe.Workflow(name=name)
    wf.base_dir = experiment_dir
    wf.connect([
        (segmentation_wf, dartel_wf, [(('new_segment.native_class_images', get_class_images, [1, 2, 3]), 'dartel2mni_input.native_class_images'),
                                  (('new_segment.dartel_input_images', get_class_images, [1, 2, 3]), 'dartelTemplate.image_files')])
    ])

    return wf
