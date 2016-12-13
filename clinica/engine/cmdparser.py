
# -*- coding: utf-8 -*-

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
from os.path import expanduser

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
        if arg is None:
            return None
        elif arg[:1] == '~':
            return expanduser(arg)
        elif arg[:1] == '.':
            return getcwd()
        else:
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
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("group_id",
                                help='Current group name')
        self._args.add_argument("-as", "--analysis_series_id",
                                help='Current analysis series name', default='default')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_threads", type=int, default=4,
                                help='Number of threads to run in parallel')
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

        preproc_wf = datagrabber_t1_spm_full_pipeline(self.absolute_path(args.bids_directory),
                                                      self.absolute_path(args.caps_directory),
                                                      self.absolute_path(args.subjects_sessions_tsv),
                                                      args.group_id,
                                                      analysis_series_id=args.analysis_series_id,
                                                      working_directory=self.absolute_path(working_directory),
                                                      tissue_classes=args.tissue_classes,
                                                      dartel_tissues=args.dartel_tissues,
                                                      save_warped_unmodulated=args.save_warped_unmodulated,
                                                      save_warped_modulated=args.save_warped_modulated,
                                                      in_write_deformation_fields=args.write_deformation_fields,
                                                      in_fwhm=args.fwhm,
                                                      in_modulate=args.modulate,
                                                      in_voxel_size=voxel_size)

        # print 'Workflow set'
        preproc_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserT1SPMSegment(CmdParser):

    def define_name(self):
        self._name = 't1-spm-segment'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-as", "--analysis_series_id",
                                help='Current analysis series name', default='default')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_threads", type=int, default=4,
                                help='Number of threads to run in parallel')
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
        segment_wf = datagrabber_t1_spm_segment_pipeline(self.absolute_path(args.bids_directory),
                                                         self.absolute_path(args.caps_directory),
                                                         self.absolute_path(args.subjects_sessions_tsv),
                                                         analysis_series_id=args.analysis_series_id,
                                                         working_directory=self.absolute_path(working_directory),
                                                         tissue_classes=args.tissue_classes,
                                                         dartel_tissues=args.dartel_tissues,
                                                         save_warped_unmodulated=args.save_warped_unmodulated,
                                                         save_warped_modulated=args.save_warped_modulated,
                                                         in_write_deformation_fields=args.write_deformation_fields)

        # print 'Workflow set'
        segment_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserT1FreeSurfer(CmdParser):

    def define_name(self):
        self._name = 't1-freesurfer'

    def define_options(self):
        self._args.add_argument("bids_dir",
                                help='Path to the BIDS directory')
        self._args.add_argument("caps_dir",
                                help='Path to the CAPS output directory')
        self._args.add_argument("subjects_sessions",
                                help='TSV file with subjects and sessions to be processed')
        self._args.add_argument("-asi", "--analysis_series_id", type=str, default='default',
                                help='Current analysis series name')
        self._args.add_argument("-wd", "--working_directory", type=str, default=None,
                                help='Temporary directory to run the workflow')
        self._args.add_argument("-ras", "--reconall_args", type=str, default='-qcache',
                                help='additional flags for recon-all command line, default is -qcache')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')

    def run_pipeline(self, args):

        from clinica.pipeline.t1.t1_freesurfer import t1_freesurfer_pipeline

        reconall_wf = t1_freesurfer_pipeline(self.absolute_path(args.bids_dir),
                                         self.absolute_path(args.caps_dir),
                                         self.absolute_path(args.subjects_sessions),
                                         analysis_series_id=args.analysis_series_id,
                                         working_directory=self.absolute_path(args.working_directory),
                                         recon_all_args=args.reconall_args)

        reconall_wf.run("MultiProc", plugin_args={'n_procs': args.n_procs})



class CmdParserStatisticsSurfStat(CmdParser):

    def define_name(self):
        self._name = 'statistics-surfstat'

    def define_options(self):
        self._args.add_argument("caps_dir",
                                help='Directory where the input files(output of FreeSurfer pipeline) are stored')
        self._args.add_argument("subjects_visits_tsv",
                                help='Directory where the tsv files are stored')
        self._args.add_argument("linear_model",
                                help='A list to define the model that fits into GLM, eg, 1 + group_label + sex + age')
        self._args.add_argument("contrast",
                                help='A list to define the contrast matrix for GLM, eg, group_label')
        self._args.add_argument("str_format",
                                help='A list to define the format string for the tsv column , eg, %%s %%s %%s %%f')
        self._args.add_argument("group_label",
                                help='Current group name')
        self._args.add_argument("analysis_series_id",
                                help='an existed recon-alled series, correspondent to your former t1_freesurfer pipeline')
        self._args.add_argument("-sof", "--size_of_fwhm", type=int, default=20, help='FWHM for the surface smoothing')
        self._args.add_argument("-tup", "--threshold_uncorrected_p_value", type=float, default='0.001',
                                help='Threshold to display the uncorrected Pvalue')
        self._args.add_argument("-tcp", "--threshold_corrected_p_value", type=float, default=0.05,
                                help='Threshold to display the corrected cluster')
        self._args.add_argument("-ct", "--cluster_threshold", type=float, default=0.001,
                                help='Threshold to define a cluster in the process of cluster-wise correction')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-wd", "--working_directory", type=str, default=None,
                                help='Temporary directory to run the workflow')

    def run_pipeline(self, args):

        from clinica.pipeline.statistics.surfstat import clinica_surfstat
        # working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
        surfstat_wf = clinica_surfstat(self.absolute_path(args.caps_dir),
                                       self.absolute_path(args.subjects_visits_tsv),
                                       args.linear_model,
                                       args.contrast,
                                       args.str_format,
                                       args.group_label,
                                       args.analysis_series_id,
                                       size_of_fwhm=args.size_of_fwhm,
                                       threshold_uncorrected_pvalue=args.threshold_uncorrected_p_value,
                                       threshold_corrected_pvalue=args.threshold_corrected_p_value,
                                       cluster_threshold=args.cluster_threshold,
                                       working_directory=self.absolute_path(args.working_directory))

        surfstat_wf.run("MultiProc", plugin_args={'n_procs': args.n_procs})


class CmdParserMachineLearningVBLinearSVM(CmdParser):

    def define_name(self):
        self._name = 'ml-vb-linear-svm'

    def define_options(self):
        self._args.add_argument("caps_directory",
                                help='Directory where the input NIFTI images are stored')
        self._args.add_argument("subjects_visits_tsv",
                                help='TSV file with subjects and sessions to be processed')
        self._args.add_argument("analysis_series_id",
                                help='Current analysis series name')
        self._args.add_argument("group_id",
                                help='Current group name')
        self._args.add_argument("diagnoses_tsv",
                                help='TSV file with subjects diagnoses')
        self._args.add_argument("-p", "--prefix",
                                help='Images prefix')
        self._args.add_argument("-mz", "--mask_zeros", type=bool, default=True,
                                help='Use a mask to remove zero valued voxels across images')
        self._args.add_argument("-b", "--balanced", type=bool, default=True,
                                help='Balance the weights of subjects for the SVM proportionally to their number in each class')
        self._args.add_argument("-cv", "--cv_folds", type=int, default=10,
                                help='Number of folds to use in the cross validation')
        self._args.add_argument("-fc", "--folds_c", type=int, default=10,
                                help='Number of folds to use in the cross validation to determine parameter C')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-crl", "--c_range_logspace", nargs=3, type=int, default=[-6, 2, 17],
                                help="numpy logspace function arguments defining the range of search for SVM parameter C. Ex: -6 2 17")
        self._args.add_argument("-sgm", "--save_gram_matrix", action='store_true',
                                help="Save feature weights for each classification as a matrix")
        self._args.add_argument("-sc", "--save_subject_classification", action='store_true',
                                help="Save list of classification results for each subject for each classification")
        self._args.add_argument("-sw", "--save_original_weights", action='store_true',
                                help="Save feature weights for each classification as a matrix")
        # self._args.add_argument("-sf", "--save_features_image", action='store_true',
        #                         help="Save feature weights for each classification as an image")

    def run_pipeline(self, args):

        from clinica.pipeline.machine_learning.voxel_based_svm import linear_svm_binary_classification_caps
        from numpy import logspace

        c_range = logspace(args.c_range_logspace[0], args.c_range_logspace[1], args.c_range_logspace[2])

        linear_svm_binary_classification_caps(self.absolute_path(args.caps_directory),
                                              self.absolute_path(args.subjects_visits_tsv),
                                              args.analysis_series_id,
                                              args.group_id,
                                              self.absolute_path(args.diagnoses_tsv),
                                              prefix=args.prefix,
                                              mask_zeros=args.mask_zeros,
                                              balanced=args.balanced,
                                              outer_folds=args.cv_folds,
                                              inner_folds=args.folds_c,
                                              n_threads=args.n_procs,
                                              c_range=c_range,
                                              save_gram_matrix=args.save_gram_matrix,
                                              save_subject_classification=args.save_subject_classification,
                                              save_original_weights=args.save_original_weights,
                                              save_features_image=args.save_features_image)


class CmdParserT1DWIRegistration(CmdParser):

    def define_name(self):
        self._name = 't1-dwi-registration'


    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-analysis_series_id", default='default',
                                help='Label for analysis series id (default name is \'default\')')
        self._args.add_argument("-n_threads", type=int, default=0,
                                help='Number of threads (default=0, which disables multi-threading).')
        group = self._args.add_mutually_exclusive_group(required=True)
        group.add_argument('-is_bias_corrected', action='store_true',
                           help='Set this flag if your images are bias corrected (mutually exclusive with \'-is_not_bias_corrected\').')
        group.add_argument('-is_not_bias_corrected', action='store_false',
                           help='Set this flag if your images are not bias corrected (mutually exclusive with \'-is_bias_corrected\').')

    def run_pipeline(self, args):
        import csv
        import os.path
        from clinica.pipeline.registration.mri_registration import t1_b0_registration_pipeline

        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            for row in tsv_reader:
                bids_path_to_t1 = os.path.join(self.absolute_path(args.bids_directory), row[0], row[1], 'anat',
                                               row[0] + '_' + row[1] + '_T1w.nii.gz')
                caps_path_to_bet_t1 = os.path.join(self.absolute_path(args.caps_directory), row[0], row[1], 't1',
                                                   'fsl', row[0] + '_' + row[1] + '_brain-extracted-T1w.nii.gz')
                caps_path_to_desikan_parcellation = os.path.join(self.absolute_path(args.caps_directory),
                                                                 row[0], row[1], 'registration', row[0] + '_' + row[1] + '_desikan-parcellation.nii.gz')
                caps_path_to_destrieux_parcellation = os.path.join(self.absolute_path(args.caps_directory),
                                                                 row[0], row[1], 'registration', row[0] + '_' + row[1] + '_destrieux-parcellation.nii.gz')
                caps_path_to_wm_mask = os.path.join(self.absolute_path(args.caps_directory),
                                                    row[0], row[1], 'dwi', row[0] + '_' + row[1] + '_binary-white-matter-mask.nii.gz')
                assert(os.path.isfile(caps_path_to_bet_t1))
                assert(os.path.isfile(caps_path_to_desikan_parcellation))
                assert(os.path.isfile(caps_path_to_destrieux_parcellation))
                assert(os.path.isfile(caps_path_to_wm_mask))

                t1_b0_registration = t1_b0_registration_pipeline(subject_id=row[0],
                                                                 session_id=row[1],
                                                                 analysis_series_id=args.analysis_series_id,
                                                                 caps_directory=self.absolute_path(args.caps_directory),
                                                                 working_directory=self.absolute_path(args.working_directory),
                                                                 )
                t1_b0_registration.inputs.inputnode.in_bias_corrected_bet_t1 = caps_path_to_bet_t1
                t1_b0_registration.inputs.inputnode.in_preprocessed_dwi = 'subject_id_preprocessed_dwi.nii'
                t1_b0_registration.inputs.inputnode.in_b0_mask = 'subject_id_b0_mask.nii'
                t1_b0_registration.inputs.inputnode.in_white_matter_binary_mask = 'subject_id_wm_mask.nii'
                t1_b0_registration.inputs.inputnode.in_desikan_parcellation = caps_path_to_desikan_parcellation
                t1_b0_registration.inputs.inputnode.in_destrieux_parcellation = caps_path_to_destrieux_parcellation
                t1_b0_registration.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserT1FSL(CmdParser):

    def define_name(self):
        self._name = 't1-fsl'


    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-analysis_series_id", default='default',
                                help='Label for analysis series id (default name is \'default\')')
        self._args.add_argument("-n_threads", type=int, default=0,
                                help='Number of threads (default=0, which disables multi-threading).')
        group = self._args.add_mutually_exclusive_group(required=True)
        group.add_argument('-is_bias_corrected', action='store_true',
                           help='Set this flag if your images are bias corrected (mutually exclusive with \'-is_not_bias_corrected\').')
        group.add_argument('-is_not_bias_corrected', action='store_false',
                           help='Set this flag if your images are not bias corrected (mutually exclusive with \'-is_bias_corrected\').')

    def run_pipeline(self, args):
        import csv
        import os.path
        from clinica.pipeline.t1.t1_fsl import t1_fsl_segmentation_pipeline

        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            # Check inputs:
            for row in tsv_reader:
                bids_path_to_t1 = os.path.join(self.absolute_path(args.bids_directory), row[0], row[1], 'anat',
                                               row[0] + '_' + row[1] + '_T1w.nii.gz')
                assert(os.path.isfile(bids_path_to_t1))
                t1_fsl_segmentation = t1_fsl_segmentation_pipeline(
                    subject_id=row[0],
                    session_id=row[1],
                    analysis_series_id=args.analysis_series_id,
                    caps_directory=self.absolute_path(args.caps_directory),
                    working_directory=self.absolute_path(args.working_directory),
                    is_bias_corrected=args.is_bias_corrected
                    )
                t1_fsl_segmentation.inputs.inputnode.in_t1 = bids_path_to_t1
                t1_fsl_segmentation.run('MultiProc', plugin_args={'n_procs': args.n_threads})



class CmdParserDWIPreprocessingFieldmapBased(CmdParser):

    def define_name(self):
        self._name = 'dwi-preprocessing-fieldmap-based'


    def define_options(self):
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-analysis_series_id", default='default',
                                help='Label for analysis series id (default name is default)')
        self._args.add_argument("-n_threads", type=int, default=0,
                                help='Number of threads (default=0, which disables multi-threading).')

    def run_pipeline(self, args):
        import csv
        import os.path
        from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_fieldmap_based

        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            # Check inputs:
            for row in tsv_reader:
                bids_path_to_dwi = os.path.join(self.absolute_path(args.caps_directory),
                                               'sub-' + row[0], 'ses-' + row[1], 'anat',
                                               'sub-' + row[0] + '_ses-' + row[1] + '_T1w.nii.gz')
                assert(os.path.isfile(bids_path_to_dwi))
                dwi_preprocessing = diffusion_preprocessing_fieldmap_based(

                )
                dwi_preprocessing.inputs.inputnode.in_t1 = bids_path_to_dwi
                dwi_preprocessing.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserDWIProcessing(CmdParser):

    def define_name(self):
        self._name = 'dwi-dti-and-tracto-for-the-moment'


    def define_options(self):
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-analysis_series_id", default='default',
                                help='Label for analysis series id (default name is default)')
        self._args.add_argument("-n_threads", type=int, default=0,
                                help='Number of threads (default=0, which disables multi-threading).')


    def run_pipeline(self, args):
        import csv
        import os.path
        from clinica.pipeline.dwi.dwi_processing import tractography_and_dti_pipeline

        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            # Check inputs:
            for row in tsv_reader:
                caps_path_to_dwi = os.path.join(
                    self.absolute_path(args.caps_directory), 'sub-' + row[0], 'ses-' + row[1], 'dwi', 'preprocessing',
                    'sub-' + row[0] + '_ses-' + row[1] + '_dwi')
                caps_path_to_b0_mask = os.path.join(
                    self.absolute_path(args.caps_directory), 'sub-' + row[0], 'ses-' + row[1], 'dwi',
                    'sub-' + row[0] + '_ses-' + row[1] + '_dwi')
                caps_path_to_white_matter_binary_mask = os.path.join(self.absolute_path(args.caps_directory),
                                                   'sub-' + row[0], 'ses-' + row[1], 'dwi',
                                                   'sub-' + row[0] + '_ses-' + row[1] + '_binary-white-matter-mask.nii.gz')

                assert(os.path.isfile(caps_path_to_dwi + '.bval'))
                assert(os.path.isfile(caps_path_to_dwi + '.bvec'))
                assert(os.path.isfile(caps_path_to_dwi + '.nii.gz'))
                assert(os.path.isfile(caps_path_to_b0_mask))

                dwi_processing = tractography_and_dti_pipeline()
                dwi_processing.inputs.inputnode.in_t1 = caps_path_to_dwi
                dwi_processing.run('MultiProc', plugin_args={'n_procs': args.n_threads})

                from clinica.pipeline.dwi.dwi_processing import tractography_and_dti_pipeline
                tractography_and_dti = tractography_and_dti_pipeline(
                    subject_id=row[0], session_id=row[1], analysis_series_id=args.analysis_series_id,
                    caps_directory=self.absolute_path(args.caps_directory),
                    working_directory=self.absolute_path(args.working_directory),
                    max_harmonic_order=None,
                    tractography_algorithm='iFOD2',
                    tractography_nb_of_tracks="100K",
                    tractography_fod_threshold=None,
                    tractography_step_size=None,
                    tractography_angle=None,
                    nthreads=2
                )
                tractography_and_dti.inputs.inputnode.in_dwi = caps_path_to_dwi + '.nii.gz'
                tractography_and_dti.inputs.inputnode.in_bvecs = caps_path_to_dwi + '.bvec'
                tractography_and_dti.inputs.inputnode.in_bvals = caps_path_to_dwi + '.bval'
                tractography_and_dti.inputs.inputnode.in_b0_mask = caps_path_to_b0_mask
                tractography_and_dti.inputs.inputnode.in_white_matter_binary_mask = caps_path_to_white_matter_binary_mask
                tractography_and_dti.run()


class CmdParserCappToBids(CmdParser):

    def define_name(self):
        self._name = 'capp-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized CAPP directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("-co", type=bool, default=False,
                                help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')

    def run_pipeline(self, args):
        from clinica.bids import capp_to_bids
        capp_to_bids.convert(args.dataset_directory, args.bids_directory)


class CmdParserInsightToBids(CmdParser):

    def define_name(self):
        self._name = 'insight-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized INSIGHT directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("-co", type=bool, default=False,
                                help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')

    def run_pipeline(self, args):
        from clinica.bids import insight_to_bids
        insight_to_bids.convert(args.dataset_directory, args.bids_directory)


class CmdParserPrevDemAlsToBids(CmdParser):

    def define_name(self):
        self._name = 'prevdemals-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized INSIGHT directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("-co", type=bool, default=False,
                                help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')

    def run_pipeline(self, args):
        from clinica.bids import prevdemals_to_bids
        prevdemals_to_bids.convert(args.dataset_directory, args.bids_directory)

class CmdParserSubsSess(CmdParser):

    def define_name(self):
        self._name = 'create-subjects-visits'

    def define_options(self):
        self._args.add_argument("bids_dir",
                               help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')


    def run_pipeline(self, args):
        from clinica.bids.utils import data_handling as dt
        dt.create_subs_sess_list(args.bids_dir, args.out_directory)


class CmdParserMergeTsv(CmdParser):

    def define_name(self):
        self._name = 'merge-tsv'

    def define_options(self):
        self._args.add_argument("bids_dir",
                               help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')
        self._args.add_argument("-tf", '--true_false_mode', type= bool, default = False,
                                help = '(Optional) Conver all the field with binary meaning into True and False values.')


    def run_pipeline(self, args):
        from clinica.bids.utils import data_handling as dt
        dt.create_merge_file(args.bids_dir, args.out_directory, args.true_false_mode)


class CmdParserMissingModalities(CmdParser):

    def define_name(self):
        self._name = 'compute-missing-mods'

    def define_options(self):
        self._args.add_argument("bids_dir",
                               help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')


    def run_pipeline(self, args):
        from clinica.bids.utils import data_handling as dt
        dt.compute_missing_mods(args.bids_dir, args.out_directory)