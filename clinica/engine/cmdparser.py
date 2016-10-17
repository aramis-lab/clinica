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
from argparse import ArgumentParser
from os.path import join
from os import getcwd

class CmdParser:
    """Abstract class to extend in order to create your command line parser

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

def init_cmdparser_objects(rootparser, parser, objects=None):
    """UTIL:Init all derived CmdParser instances with specific data

    :param parser: the ArgParser node
    :param objects: all CmpParser instances of this file
    :return: all CmpParser instances of this file
    """

    def init(x):
        def silent_help(): pass
        def error_message(p):
            def error(x):
                p.print_help()
                rootparser.print_help = silent_help
                exit(-1)
            return error

        x.options = parser.add_parser(x.name)
        x.options.error = error_message(x.options)
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
        self._args.add_argument("input_directory",
                                help='Directory where the input NIFTI images are stored')
        self._args.add_argument("output_directory",
                                help='Directory to save the resulting images')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to run the workflow')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
                                help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
                                help="Save warped unmodulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
                                help="Save warped modulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
                                help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")
        # The default smoothing is 8mm isotropic, as in the Matlab version of SPM.
        # This is recommended when using modulated images.
        self._args.add_argument("-fwhm", "--fwhm", nargs=3, type=float, default=[8, 8, 8],
                                help="A list of 3 floats specifying the FWHM for each dimension")
        self._args.add_argument("-m", "--modulate", type=bool, default=True,
                                help='A boolean. Modulate output images - no modulation preserves concentrations')
        self._args.add_argument("-vs", "--voxel_size", nargs=3, type=float,
                                help="A list of 3 floats specifying voxel sizes for each dimension of output image")


    def run_pipeline(self, args):

        from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_full_pipeline

        working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None

        voxel_size = tuple(args.voxel_size) if args.voxel_size is not None else None

        preproc_wf = datagrabber_t1_spm_full_pipeline(self.absolute_path(args.input_directory),
                                                      self.absolute_path(args.output_directory),
                                                      working_directory=working_directory,
                                                      tissue_classes=args.tissue_classes,
                                                      dartel_tissues=args.dartel_tissues,
                                                      save_warped_unmodulated=args.save_warped_unmodulated,
                                                      save_warped_modulated=args.save_warped_modulated,
                                                      in_write_deformation_fields=args.write_deformation_fields,
                                                      in_fwhm=args.fwhm,
                                                      in_modulate=args.modulate,
                                                      in_voxel_size=voxel_size)

        print 'Workflow set'
        preproc_wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})


class CmdParserT1SPMSegment(CmdParser):

    def define_name(self):
        self._name = 't1-spm-segment'

    def define_options(self):
        self._args.add_argument("input_directory",
                                help='Directory where the input NIFTI images are stored')
        self._args.add_argument("output_directory",
                                help='Directory to save the resulting images')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to run the workflow')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
                                help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
                                help="Save warped unmodulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
                                help="Save warped modulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
                                help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")

    def run_pipeline(self, args):

        from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_segment_pipeline

        working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
        segment_wf = datagrabber_t1_spm_segment_pipeline(self.absolute_path(args.input_directory),
                                                         self.absolute_path(args.output_directory),
                                                         working_directory=working_directory,
                                                         tissue_classes=args.tissue_classes,
                                                         dartel_tissues=args.dartel_tissues,
                                                         save_warped_unmodulated=args.save_warped_unmodulated,
                                                         save_warped_modulated=args.save_warped_modulated,
                                                         in_write_deformation_fields=args.write_deformation_fields)

        print 'Workflow set'
        segment_wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})


class CmdParserT1SPMSegmentBIDS(CmdParser):

    def define_name(self):
        self._name = 't1-spm-segment-bids'

    def define_options(self):
        self._args.add_argument("input_directory",
                                help='Directory where the input NIFTI images are stored')
        self._args.add_argument("output_directory",
                                help='Directory to save the resulting images')
        self._args.add_argument("subjects_visits_tsv",
                                help='TSV file with subjects and sessions to be processed')
        self._args.add_argument("analysis_series_id",
                                help='Current analysis series name')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to run the workflow')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
                                help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
                                help="Save warped unmodulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
                                help="Save warped modulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
                                help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")

    def run_pipeline(self, args):

        from clinica.pipeline.t1.t1_spm_bids import datagrabber_t1_spm_segment_pipeline_bids

        working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
        segment_wf = datagrabber_t1_spm_segment_pipeline_bids(self.absolute_path(args.input_directory),
                                                         self.absolute_path(args.output_directory),
                                                         self.absolute_path(args.subjects_visits_tsv),
                                                         args.analysis_series_id,
                                                         working_directory=working_directory,
                                                         tissue_classes=args.tissue_classes,
                                                         dartel_tissues=args.dartel_tissues,
                                                         save_warped_unmodulated=args.save_warped_unmodulated,
                                                         save_warped_modulated=args.save_warped_modulated,
                                                         in_write_deformation_fields=args.write_deformation_fields)

        print 'Workflow set'
        segment_wf.run('MultiProc', plugin_args={'n_procs': args.n_procs})


class CmdParserT1ReconAll(CmdParser):

    def define_name(self):
        self._name = 't1-reconall'

    def define_options(self):
        self._args.add_argument("input_directory", help='Directory where the NIFTI images are stored')
        self._args.add_argument("output_dir", help='Directory to store the result of the pipeline')
        self._args.add_argument("field_template", help='A list to define the input structure')
        self._args.add_argument("template_args", help='A list of list to define the input structure, the name of the NIFTI images')
        self._args.add_argument("-ifs", "--intermediate_files", type=list, default=['orig', 'white'], help='The intermediate files stored in datasinker')
        self._args.add_argument("-ras", "--reconall_args", type=str, default='-qcache', help='additional flags for reconAll command line, default is -qcache')

    def run_pipeline(self, args):

        from clinica.pipeline.t1.t1_freesurfer import recon_all_pipeline

        reconall_wf = recon_all_pipeline(self.absolute_path(args.input_directory), self.absolute_path(args.output_dir), args.field_template, args.template_args,
                                         datasink_para=args.intermediate_files, recon_all_args=args.reconall_args)

        reconall_wf.run("MultiProc", plugin_args={'n_procs':4})

class CmdParserStatisticsSurfStat(CmdParser):

    def define_name(self):
        self._name = 't1-surfstat'

    def define_options(self):
        self._args.add_argument("input_directory", help='Directory where the input files(output of reconAll pipeline) are stored')
        self._args.add_argument("output_dir", help='Directory to store the result images of the pipeline')
        self._args.add_argument("linear_model", help='A list to define the model that fits into GLM')
        self._args.add_argument("contrast", help='A list to define the contrast matrix for GLM')
        self._args.add_argument("csv_file", help='Directory where the csv files are stored')
        self._args.add_argument("str_format", help='A list to define the format string for the csv files')
        self._args.add_argument("-sof", "--size_of_fwhm", type=int, default=20, help='FWHM for the surface smoothing')
        self._args.add_argument("-tup", "--threshold_uncorrected_pvalue", type=float, default='0.001', help='Threshold to display the uncorrected Pvalue')
        self._args.add_argument("-tcp", "--threshold_corrected_pvalue", type=float, default=0.05, help='Threshold to display the corrected cluster')
        self._args.add_argument("-ct", "--cluster_threshold", type=float, default=0.001, help='Threshold to define a cluster in the process of cluster-wise correction')

    def run_pipeline(self, args):

        from clinica.pipeline.statistics.t1_surfstat import clinica_surfstat
        
        surfstat_wf = clinica_surfstat(self.absolute_path(args.input_directory), self.absolute_path(args.output_dir), args.linear_model, args.contrast,
                                         self.absolute_path(args.csv_file), args.str_format,
                                         size_of_fwhm=args.size_of_fwhm, threshold_uncorrected_pvalue=args.threshold_uncorrected_pvalue,
                                         threshold_corrected_pvalue=args.threshold_corrected_pvalue, cluster_threshold=args.cluster_threshold)

        surfstat_wf.run()
