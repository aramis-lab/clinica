"""Define the command line parser for each pipeline

Please, add your command line parser at the end of this file.

Just take a look to the following example:
______________________________________________________________
class CmdParserT1(CmdParser):

    def define_name(self):
        self._name = 'T1'

    def define_options(self):
        self._args.add_argument("-s", "--source", dest='source')

    def run_pipeline(self, args):
        print "run pipeline %s" % args.source
______________________________________________________________

"""

import abc
from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_full_pipeline, datagrabber_t1_spm_segment_pipeline
from clinica.pipeline.t1.t1_freesurfer import recon_all_pipeline
from argparse import ArgumentParser
from os.path import join
from os import getcwd

class CmdParser:
    """Abstract class to exend in order to create your command line parser

    Take a look to the CmdParserT1 example and write your own
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.reset()
        self.build()

    def build(self):
        self.define_name()
        self.define_options()

    def reset(self):
        self._args = ArgumentParser()
        self._name = None

    @property
    def options(self): return self._args

    @options.setter
    def options(self, x): self._args = x

    @property
    def name(self): return self._name

    @name.setter
    def name(self, x): self._name = x

    @abc.abstractmethod
    def define_name(self): pass

    @abc.abstractmethod
    def define_options(self): pass

    @abc.abstractmethod
    def run_pipeline(self, args): pass

    def absolute_path(self, arg):
        return join(getcwd(), arg)





def get_cmdparser_objects():
    """UTIL: read all derived classes of CmdParser

    :return: all derived instance of CmdParser present to this file
    """

    import inspect
    import clinica.engine.cmdparser
    for name, obj in inspect.getmembers(clinica.engine.cmdparser):
        if name != 'CmdParser' and inspect.isclass(obj):
            x = obj()
            if isinstance(x, clinica.engine.cmdparser.CmdParser):
                yield x

def init_cmdparser_objects(parser,objects=None):
    """UTIL:Init all derived CmdParser instances with specific data

    :param parser: the ArgParser node
    :param objects: all CmpParser instances of this file
    :return: all CmpParser instances of this file
    """
    def init(x):
        x.options = parser.add_parser(x.name)
        x.options.set_defaults(func=x.run_pipeline)
        x.build()
    if objects is None: objects = get_cmdparser_objects()
    [init(x) for x in objects]

def get_cmdparser_names(objects=None):
    """Return the names of all pipelines

    :param objects: all CmpParser instances of this file
    :return: the names of all pipelines
    """
    if objects is None: objects = get_cmdparser_objects()
    for x in objects: yield x.name


#______ _   _ _____   _   _  ___________ _____  __   _______ _   _______   _____  _       ___   _____ _____
#| ___ \ | | |_   _| | | | ||  ___| ___ \  ___| \ \ / /  _  | | | | ___ \ /  __ \| |     / _ \ /  ___/  ___|
#| |_/ / | | | | |   | |_| || |__ | |_/ / |__    \ V /| | | | | | | |_/ / | /  \/| |    / /_\ \\ `--.\ `--.
#|  __/| | | | | |   |  _  ||  __||    /|  __|    \ / | | | | | | |    /  | |    | |    |  _  | `--. \`--. \
#| |   | |_| | | |   | | | || |___| |\ \| |___    | | \ \_/ / |_| | |\ \  | \__/\| |____| | | |/\__/ /\__/ /
#\_|    \___/  \_/   \_| |_/\____/\_| \_\____/    \_/  \___/ \___/\_| \_|  \____/\_____/\_| |_/\____/\____/

class CmdParserT1SPMFullPrep(CmdParser):

    def define_name(self):
        self._name = 't1-spm-full-prep'

    def define_options(self):
        self._args.add_argument("input_dir", help='Directory where the NIFTI images are stored')
        self._args.add_argument("experiment_dir", help='Directory to run the workflow')
        self._args.add_argument("datasink_dir",
                                help='Directory to save the resulting images of segmentation and registration processes')
        self._args.add_argument("-np", "--n_procs", type=int, default=4, help='Number of parallel process to run')
        self._args.add_argument("-ci", "--class_images", type=list, default=[1, 2, 3],
                                help='Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF')
        self._args.add_argument("-dci", "--dartel_class_images", type=list, default=[1],
                                help='Classes of images to use for DARTEL template calculation. Ex: [1] is only GM')
        self._args.add_argument("-af", "--affine_regularization",
                                help="('mni' or 'eastern' or 'subj' or 'none')")
        self._args.add_argument("-chi", "--channel_info", type=tuple,
                                help='''a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
        A tuple with the following fields:
         - bias reguralisation (0-10)
         - FWHM of Gaussian smoothness of bias
         - which maps to save (Corrected, Field) - a tuple of two boolean values''')
        self._args.add_argument("-sd", "--sampling_distance", type=float,
                                help="Sampling distance on data for parameter estimation")
        self._args.add_argument("-ts", "--tissues_to_save", type=list,
                                help='''A list of tuples (one per tissue) of the form:
        ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated)) with the boolen value for the type of images to save
        Ex.: [((True, True), (False, False)),
              ((True, False), (False, False))]''')
        self._args.add_argument("-wr", "--warping_regularization", type=list,
                                help='''(a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)''')
        self._args.add_argument("-wdf", "--write_deformation_fields", type=list,
                                help='''(a list of from 2 to 2 items which are a boolean)
        Which deformation fields to write:[Inverse, Forward]''')
        self._args.add_argument("-ip", "--iteration_parameters", type=list,
                                help='''(a list of from 3 to 12 items which are a tuple
         of the form: (1 <= an integer <= 10, a tuple of the form: (a float,
         a float, a float), 1 or 2 or 4 or 8 or 16 or 32 or 64 or 128 or 256
         or 512, 0 or 0.5 or 1 or 2 or 4 or 8 or 16 or 32))
        List of tuples for each iteration
         - Inner iterations
         - Regularization parameters
         - Time points for deformation model
         - smoothing parameter''')
        self._args.add_argument("-op", "--optimization_parameters", type=tuple,
                                help='''a tuple of the form: (a float, 1 <= an integer <= 8, 1 <= an integer <= 8))
         Optimization settings a tuple
         - LM regularization
         - cycles of multigrid solver
         - relaxation iterations''')
        self._args.add_argument("-rf", "--regularization_form",
                                help='''('Linear' or 'Membrane' or 'Bending')
        Form of regularization energy term''')
        self._args.add_argument("-tp", "--template_prefix",
                                help='''(a string, nipype default value: Template)
        Prefix for template''')
        self._args.add_argument("-bb", "--bounding_box", type=tuple,
                                help='''(a tuple of the form: (a float, a float, a float, a float, a float, a float))
        Voxel sizes for output file''')
        self._args.add_argument("-fwhm", "--fwhm", type=list, default=[0, 0, 0],
                                help='''(a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension''')
        self._args.add_argument("-m", "--modulate", type=bool, default=True,
                                help='Modulate out images - no modulation preserves concentrations')
        self._args.add_argument("-vs", "--voxel_size", type=tuple,
                                help='''(a tuple of the form: (a float, a float, a float))
        Voxel sizes for output file''')

    def run_pipeline(self, args):

        preproc_wf = datagrabber_t1_spm_full_pipeline(args.input_dir, args.experiment_dir, args.datasink_dir,
                                          class_images=args.class_images, dartel_class_images=args.dartel_class_images,
                                          in_affine_regularization=args.affine_regularization, in_channel_info=args.channel_info,
                                          in_sampling_distance=args.sampling_distance, in_tissues_to_save=args.tissues_to_save,
                                          in_warping_regularization=args.warping_regularization, in_write_deformation_fields=args.write_deformation_fields,
                                          in_iteration_parameters=args.iteration_parameters, in_optimization_parameters=args.optimization_parameters,
                                          in_regularization_form=args.regularization_form, in_template_prefix=args.template_prefix,
                                          in_bounding_box=args.bounding_box, in_fwhm=args.fwhm, in_modulate=args.modulate,
                                          in_voxel_size=args.voxel_size)

        preproc_wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})



class CmdParserT1SPMSegment(CmdParser):

    def define_name(self):
        self._name = 't1-spm-segment'

    def define_options(self):
        self._args.add_argument("input_dir", help='Directory where the NIFTI images are stored')
        self._args.add_argument("experiment_dir", help='Directory to run the workflow')
        self._args.add_argument("datasink_dir",
                                help='Directory to save the resulting images of segmentation and registration processes')
        self._args.add_argument("-np", "--n_procs", type=int, default=4, help='Number of parallel process to run')
        self._args.add_argument("-ci", "--class_images", type=list, default=[1, 2, 3],
                                help='Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF')
        self._args.add_argument("-af", "--affine_regularization",
                                help="('mni' or 'eastern' or 'subj' or 'none')")
        self._args.add_argument("-chi", "--channel_info", type=tuple,
                                help='''a tuple of the form: (a float, a float, a tuple of the form: (a boolean, a boolean)))
        A tuple with the following fields:
         - bias reguralisation (0-10)
         - FWHM of Gaussian smoothness of bias
         - which maps to save (Corrected, Field) - a tuple of two boolean values''')
        self._args.add_argument("-sd", "--sampling_distance", type=float,
                                help="Sampling distance on data for parameter estimation")
        self._args.add_argument("-ts", "--tissues_to_save", type=list,
                                help='''A list of tuples (one per tissue) of the form:
        ((Native space, DARTEL input),(Warped Unmodulated, Warped Modulated)) with the boolen value for the type of images to save
        Ex.: [((True, True), (False, False)),
              ((True, False), (False, False))]''')
        self._args.add_argument("-wr", "--warping_regularization", type=list,
                                help='''(a list of from 5 to 5 items which are a float or a float)
        Warping regularization parameter(s). Accepts float or list of floats (the latter is required by SPM12)''')
        self._args.add_argument("-wdf", "--write_deformation_fields", type=list,
                                help='''(a list of from 2 to 2 items which are a boolean)
        Which deformation fields to write:[Inverse, Forward]''')

    def run_pipeline(self, args):

        segment_wf = datagrabber_t1_spm_segment_pipeline(args.input_dir, args.experiment_dir, args.datasink_dir,
                                                         class_images=args.class_images, in_affine_regularization=args.affine_regularization,
                                                         in_channel_info=args.channel_info, in_sampling_distance=args.sampling_distance,
                                                         in_tissues_to_save=args.tissues_to_save, in_warping_regularization=args.warping_regularization,
                                                         in_write_deformation_fields=args.write_deformation_fields,)

        segment_wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})

class CmdParserT1ReconAll(CmdParser):

    def define_name(self):
        self._name = 't1-reconall'
    def define_options(self):
        self._args.add_argument("input_dir", help='Directory where the NIFTI images are stored')
        self._args.add_argument("output_dir", help='Directory to store the result of the pipeline')
        self._args.add_argument("field_template", help='A list to define the input structure')
        self._args.add_argument("template_args", type=str, help='A list of list to define the input structure, the name of the NIFTI images')
        self._args.add_argument("-ifs", "--intermediate_files", type=list, default=['orig', 'white'], help='The intermediate files stored in datasinker')
        self._args.add_argument("-ras", "--reconall_args", type=str, default='-qcache', help='additional flags for reconAll command line, default is -qcache')


    def run_pipeline(self, args):

        reconall_wf = recon_all_pipeline(self.absolute_path(args.input_dir), args.output_dir, args.field_template, args.template_args,
                                         datasink_para=args.intermediate_files, recon_all_args=args.reconall_args)

        reconall_wf.run("MultiProc", plugin_args={'n_procs':4})