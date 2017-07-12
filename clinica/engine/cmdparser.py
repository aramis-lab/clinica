
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

def init_cmdparser_objects(rootparser, parser, objects):
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
    for x in objects:
        try:
            init(x)
        except:pass


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

# class CmdParserT1SPMFullPrep(CmdParser):
#
#     def define_name(self):
#         self._name = 't1-spm-full-prep'
#
#     def define_options(self):
#         self._args.add_argument("bids_directory",
#                                 help='Path to the BIDS directory.')
#         self._args.add_argument("caps_directory",
#                                 help='Path to the CAPS directory.')
#         self._args.add_argument("subjects_sessions_tsv",
#                                 help='TSV file containing the subjects with their sessions.')
#         self._args.add_argument("group_id",
#                                 help='Current group name')
#         self._args.add_argument("-wd", "--working_directory",
#                                 help='Temporary directory to store pipeline intermediate results')
#         self._args.add_argument("-np", "--n_threads", type=int, default=4,
#                                 help='Number of threads to run in parallel')
#         self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
#                                 help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
#         self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
#                                 help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
#         self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
#                                 help="Save warped unmodulated images for tissues specified in --tissue_classes")
#         self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
#                                 help="Save warped modulated images for tissues specified in --tissue_classes")
#         self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
#                                 help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")
#         # The default smoothing is 8mm isotropic, as in the Matlab version of SPM.
#         # This is recommended when using modulated images.
#         self._args.add_argument("-fwhm", "--fwhm", nargs=3, type=float, default=[8, 8, 8],
#                                 help="A list of 3 floats specifying the FWHM for each dimension")
#         self._args.add_argument("-m", "--modulate", type=bool, default=True,
#                                 help='A boolean. Modulate output images - no modulation preserves concentrations')
#         self._args.add_argument("-vs", "--voxel_size", nargs=3, type=float,
#                                 help="A list of 3 floats specifying voxel sizes for each dimension of output image")
#
#
#     def run_pipeline(self, args):
#
#         from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_full_pipeline
#
#         working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
#
#         voxel_size = tuple(args.voxel_size) if args.voxel_size is not None else None
#
#         preproc_wf = datagrabber_t1_spm_full_pipeline(self.absolute_path(args.bids_directory),
#                                                       self.absolute_path(args.caps_directory),
#                                                       self.absolute_path(args.subjects_sessions_tsv),
#                                                       args.group_id,
#                                                       working_directory=self.absolute_path(working_directory),
#                                                       tissue_classes=args.tissue_classes,
#                                                       dartel_tissues=args.dartel_tissues,
#                                                       save_warped_unmodulated=args.save_warped_unmodulated,
#                                                       save_warped_modulated=args.save_warped_modulated,
#                                                       in_write_deformation_fields=args.write_deformation_fields,
#                                                       in_fwhm=args.fwhm,
#                                                       in_modulate=args.modulate,
#                                                       in_voxel_size=voxel_size)
#
#         # print 'Workflow set'
#         preproc_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})

#
# class CmdParserT1SPMSegment(CmdParser):
#
#     def define_name(self):
#         self._name = 't1-spm-segment'
#
#     def define_options(self):
#         self._args.add_argument("bids_directory",
#                                 help='Path to the BIDS directory.')
#         self._args.add_argument("caps_directory",
#                                 help='Path to the CAPS directory.')
#         self._args.add_argument("subjects_sessions_tsv",
#                                 help='TSV file containing the subjects with their sessions.')
#         self._args.add_argument("-wd", "--working_directory",
#                                 help='Temporary directory to store pipeline intermediate results')
#         self._args.add_argument("-np", "--n_threads", type=int, default=4,
#                                 help='Number of threads to run in parallel')
#         self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
#                                 help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
#         self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
#                                 help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
#         self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
#                                 help="Save warped unmodulated images for tissues specified in --tissue_classes")
#         self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
#                                 help="Save warped modulated images for tissues specified in --tissue_classes")
#         self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
#                                 help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")
#
#     def run_pipeline(self, args):
#
#         from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_segment_pipeline
#
#         working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
#         segment_wf = datagrabber_t1_spm_segment_pipeline(self.absolute_path(args.bids_directory),
#                                                          self.absolute_path(args.caps_directory),
#                                                          self.absolute_path(args.subjects_sessions_tsv),
#                                                          working_directory=self.absolute_path(working_directory),
#                                                          tissue_classes=args.tissue_classes,
#                                                          dartel_tissues=args.dartel_tissues,
#                                                          save_warped_unmodulated=args.save_warped_unmodulated,
#                                                          save_warped_modulated=args.save_warped_modulated,
#                                                          in_write_deformation_fields=args.write_deformation_fields)
#
#         # print 'Workflow set'
#         segment_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})
#
#
# class CmdParserT1SPMDartelTemplate(CmdParser):
#
#     def define_name(self):
#         self._name = 't1-spm-dartel-template'
#
#     def define_options(self):
#         self._args.add_argument("caps_directory",
#                                 help='Path to the CAPS directory.')
#         self._args.add_argument("subjects_sessions_tsv",
#                                 help='TSV file containing the subjects with their sessions.')
#         self._args.add_argument("group_id",
#                                 help='Current group name')
#         self._args.add_argument("-wd", "--working_directory",
#                                 help='Temporary directory to store pipeline intermediate results')
#         self._args.add_argument("-np", "--n_threads", type=int, default=4,
#                                 help='Number of threads to run in parallel')
#         self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1], choices=range(1, 7),
#                                 help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
#
#     def run_pipeline(self, args):
#
#         from clinica.pipeline.t1.t1_spm import datagrabber_t1_spm_dartel_template
#
#         working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
#
#         dartel_template_wf = datagrabber_t1_spm_dartel_template(self.absolute_path(args.caps_directory),
#                                                                 self.absolute_path(args.subjects_sessions_tsv),
#                                                                 args.group_id,
#                                                                 working_directory=self.absolute_path(working_directory),
#                                                                 dartel_tissues=args.dartel_tissues)
#
#         # print 'Workflow set'
#         dartel_template_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserPETPreprocessing(CmdParser):

    def define_name(self):
        self._name = 'pet-preprocessing'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("group_id",
                                help='Current group name. Used to obtain t1 images template and transformation.')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_threads", type=int, default=4,
                                help='Number of threads to run in parallel')
        self._args.add_argument("-pt", "--pet_type", default='FDG', choices=['FDG', 'AV45'],
                                help="Type of PET scan. Possible values are FDG or AV45.")
        self._args.add_argument("-pvc", "--pvc", action='store_true',
                                help="To apply or not partial value correction to the scan. If yes, FWHM is required.")
        self._args.add_argument("-fwhm", "--pvc_fwhm", nargs=3, type=float, default=[6, 6, 6],
                                help="A list of 3 floats specifying the FWHM for each dimension X, Y, Z")
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes to compose the final mask to apply (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...). Default is GM, WM, CSF. Up to 6 tissue classes previously obtained can be used. Ex: 1 2 3 is GM, WM and CSF")

    def run_pipeline(self, args):

        from clinica.pipeline.pet.pet import bids_caps_pet_pipeline

        working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None

        pet_wf = bids_caps_pet_pipeline(self.absolute_path(args.bids_directory),
                                        self.absolute_path(args.caps_directory),
                                        self.absolute_path(args.subjects_sessions_tsv),
                                        args.group_id,
                                        pet_type=args.pet_type,
                                        working_directory=self.absolute_path(working_directory),
                                        pvc=args.pvc,
                                        fwhm_x=args.pvc_fwhm[0],
                                        fwhm_y=args.pvc_fwhm[1],
                                        fwhm_z=args.pvc_fwhm[2],
                                        mask_tissues=args.tissue_classes)

        pet_wf.run('MultiProc', plugin_args={'n_procs': args.n_threads})


# class CmdParserT1FreeSurfer(CmdParser):
#
#     def define_name(self):
#         self._name = 't1-freesurfer'
#
#     def define_options(self):
#         self._args.add_argument("bids_dir",
#                                 help='Path to the BIDS directory')
#         self._args.add_argument("caps_dir",
#                                 help='Path to the CAPS output directory')
#         self._args.add_argument("-tsv", "--subjects_visits_tsv", type=str, default=None,
#                                 help='The path to the tsv file, which by default contains all the subjects in your BIDS dataset')
#         self._args.add_argument("-wd", "--working_directory", type=str, default=None,
#                                 help='Path to contain the information about your workflow')
#         self._args.add_argument("-ras", "--reconall_args", type=str, default='-qcache',
#                                 help='additional flags for recon-all command line, default is -qcache')
#         self._args.add_argument("-np", "--n_procs", type=int, default=4,
#                                 help='Number of parallel processes to run')
#
#     def run_pipeline(self, args):
#
#         from clinica.pipeline.t1.t1_freesurfer import datagrabber_t1_freesurfer_pipeline
#
#         reconall_wf = datagrabber_t1_freesurfer_pipeline(self.absolute_path(args.bids_dir),
#                                          self.absolute_path(args.caps_dir),
#                                          subjects_visits_tsv=self.absolute_path(args.subjects_visits_tsv),
#                                          working_directory=self.absolute_path(args.working_directory),
#                                          recon_all_args=args.reconall_args)
#
#         reconall_wf.run("MultiProc", plugin_args={'n_procs': args.n_procs})
#
#
# class CmdParserStatisticsSurfStat(CmdParser):
#
#     def define_name(self):
#         self._name = 'statistics-surfstat'
#
#     def define_options(self):
#         self._args.add_argument("caps_dir",
#                                 help='Directory where the input files(output of FreeSurfer pipeline) are stored')
#         self._args.add_argument("subjects_visits_tsv",
#                                 help='Directory where the tsv files are stored, this is based on your GLM')
#         self._args.add_argument("design_matrix",
#                                 help='A str to define the design matrix that fits into GLM, eg, 1 + group + sex + age')
#         self._args.add_argument("contrast",
#                                 help='A str to define the contrast matrix for GLM, eg, group_label')
#         self._args.add_argument("str_format",
#                                 help='A str to define the format string for the tsv column , eg, %%s %%s %%s %%f')
#         self._args.add_argument("group_label",
#                                 help='A str for current group name')
#         self._args.add_argument("glm_type",
#                                 help='A str based on glm type for the hypothesis, one between group_comparison and correlation' )
#         self._args.add_argument("-fwhm", "--full_width_at_half_maximum", type=int, default=20, help='FWHM for the surface smoothing (default=20)')
#         self._args.add_argument("-tup", "--threshold_uncorrected_p_value", type=float, default=0.001,
#                                 help='Threshold to display the uncorrected Pvalue (default=0.001)')
#         self._args.add_argument("-tcp", "--threshold_corrected_p_value", type=float, default=0.05,
#                                 help='Threshold to display the corrected cluster (default=0.05)')
#         self._args.add_argument("-ct", "--cluster_threshold", type=float, default=0.001,
#                                 help='Threshold to define a cluster in the process of cluster-wise correction (default=0.001)')
#         self._args.add_argument("-np", "--n_procs", type=int, default=4,
#                                 help='Number of parallel processes to run (default=4)')
#         self._args.add_argument("-wd", "--working_directory", type=str, default=None,
#                                 help='Temporary directory to run the workflow')
#
#     def run_pipeline(self, args):
#
#         from clinica.pipeline.statistics.surfstat import clinica_surfstat
#         # working_directory = self.absolute_path(args.working_directory) if (args.working_directory is not None) else None
#         surfstat_wf = clinica_surfstat(self.absolute_path(args.caps_dir),
#                                        self.absolute_path(args.subjects_visits_tsv),
#                                        args.design_matrix,
#                                        args.contrast,
#                                        args.str_format,
#                                        args.group_label,
#                                        args.glm_type,
#                                        full_width_at_half_maximum=args.full_width_at_half_maximum,
#                                        threshold_uncorrected_pvalue=args.threshold_uncorrected_p_value,
#                                        threshold_corrected_pvalue=args.threshold_corrected_p_value,
#                                        cluster_threshold=args.cluster_threshold,
#                                        working_directory=self.absolute_path(args.working_directory))
#
#         surfstat_wf.run("MultiProc", plugin_args={'n_procs': args.n_procs})


class CmdParserMachineLearningVBLinearSVM(CmdParser):

    def define_name(self):
        self._name = 'machinelearning-svm-voxel'

    def define_options(self):
        self._args.add_argument("image_type",
                                help='it can assume two values: pet/t1, according to the images used')
        self._args.add_argument("caps_directory",
                                help='Directory where the input NIFTI images are stored')
        self._args.add_argument("group_id",
                                help='Current group name')
        self._args.add_argument("diagnoses_tsv",
                                help='TSV file with subjects diagnoses')
        self._args.add_argument("-p", "--prefix",
                                help='Images prefix')
        self._args.add_argument("-t", "--tissue",
                                help='')
        self._args.add_argument("-svt", "--subjects_visits_tsv", type=str, default=None,
                                    help='TSV file with subjects and sessions to be processed')
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
        self._args.add_argument("-sf", "--save_features_image", action='store_true',
                                help="Save feature weights for each classification as an image")
        self._args.add_argument("-sdc", "--save_dual_coefficients", action='store_true',
                                help="Save ")
    def run_pipeline(self, args):

        from clinica.pipeline.machine_learning.voxel_based_svm import linear_svm_binary_classification_caps
        from numpy import logspace

        if args.subjects_visits_tsv is None:
            subjects_visits_tsv=() #where it's saved for t1 and pet
        else:
            subjects_visits_tsv = pandas.io.parsers.read_csv(self.absolute_path(args.participants_sessions_tsv),
                                                             sep='\t')

        c_range = logspace(args.c_range_logspace[0], args.c_range_logspace[1], args.c_range_logspace[2])

        linear_svm_binary_classification_caps(self.absolute_path(args.caps_directory),
                                              self.absolute_path(args.subjects_visits_tsv),
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

class CmdParserMachineLearningSVMRB(CmdParser):

    def define_name(self):
        self._name = 'machinelearning-svm-region'

    def define_options(self):
        self._args.add_argument("image_type",
                                   help='it can assume two values: pet/t1, according to the images used')
        self._args.add_argument("caps_directory",
                                    help='Directory where the input NIFTI images are stored')
        self._args.add_argument("group_id",
                                    help='Current group name')
        self._args.add_argument("diagnosis_tsv",
                                    help='TSV file with subjects diagnosis')
        self._args.add_argument("atlas_id",
                                    help='Name of the atlas used to extract features')
        self._args.add_argument("-svt", "--subjects_visits_tsv", type=str, default=None,
                                    help='TSV file with subjects and sessions to be processed')
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
                                    help="Save gram matrix for each classification as a matrix")
        self._args.add_argument("-sdc", "--save_dual_coefficients", action='store_true',
                                    help="Save ")
        self._args.add_argument("-sc", "--save_subject_classification", action='store_true',
                                    help="Save list of classification results for each subject for each classification")
        self._args.add_argument("-sw", "--save_original_weights", action='store_true',
                                    help="Save feature weights for each classification as a matrix")
        self._args.add_argument("-sf", "--save_features_image", action='store_true',
                                    help="Save feature weights for each classification as an image")

    def run_pipeline(self, args):
        from clinica.pipeline.machine_learning.region_based_svm import svm_binary_classification
        from clinica.pipeline.machine_learning.region_based_io import get_caps_pet_list, get_caps_t1_list, load_data
        from clinica.pipeline.machine_learning.svm_utils import gram_matrix_linear
        from numpy import logspace
        import pandas
        import os.path

        output_directory = os.path.join(self.absolute_path(args.caps_directory),
                                            'group-' +args.group_id + '/machine_learning/region_based_svm/',
                                            'space' + args.atlas_id, args.image_type)

        if args.subjects_visits_tsv is None:
            subjects_visits_tsv=() #where it's saved for t1 and pet
        else:
            subjects_visits_tsv = pandas.io.parsers.read_csv(self.absolute_path(args.participants_sessions_tsv),
                                                             sep='\t')

        if args.image_type == 't1':

            image_list = get_caps_t1_list(self.absolute_path(args.caps_directory),
                                              subjects_visits_tsv,
                                              args.group_id,
                                              args.atlas_id)
        else:

            image_list = get_caps_pet_list(self.absolute_path(args.caps_directory),
                                               subjects_visits_tsv,
                                               args.group_id,
                                               args.atlas_id)

        data = load_data(image_list, subjects_visits_tsv)
        input_image_atlas = os.path.join(os.path.expandvars('$CLINICA_HOME'),clinica, resources,atlases_spm, args.atlas_id+'.nii')
        subjects_diagnosis = pandas.io.parsers.read_csv(args.diagnosis_tsv, sep='\t')
        if list(subjects_diagnosis.columns.values) != ['participant_id', 'diagnosis']:
            raise Exception('Subjects and visits file is not in the correct format.')
        diagnosis_list = list(subjects_diagnosis.diagnosis)
        gram_matrix = gram_matrix_linear(data)
        c_range = logspace(args.c_range_logspace[0], args.c_range_logspace[1], args.c_range_logspace[2])



        svm_binary_classification(input_image_atlas,image_list,diagnosis_list,output_directory, kernel_function=None, existing_gram_matrix=gram_matrix, mask_zeros=True,
                                  scale_data=False, balanced=False,
                                  outer_folds=args.cv_folds,
                                  inner_folds=args.folds_c,
                                  n_threads=args.n_procs,
                                  c_range=c_range,
                                  save_gram_matrix=args.save_gram_matrix,
                                  save_subject_classification=args.save_subject_classification,
                                  save_dual_coefficients=args.save_dual_coefficients,
                                  scaler=None, data_mask=None,
                                  save_original_weights=args.save_original_weights,
                                  save_features_image=args.save_features_image)




class CmdParserT1FSL(CmdParser):

    def define_name(self):
        self._name = 't1-fsl'


    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("participants_sessions_tsv",
                                help='TSV file containing the participants with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-n_threads", type=int, default=1,
                                help='Number of threads (default=1, which disables multi-threading).')
        group = self._args.add_mutually_exclusive_group(required=True)
        group.add_argument('-is_bias_corrected', action='store_true',
                           help='Set this flag if your images are bias corrected (mutually exclusive with \'-is_not_bias_corrected\').')
        group.add_argument('-is_not_bias_corrected', action='store_false',
                           help='Set this flag if your images are not bias corrected (mutually exclusive with \'-is_bias_corrected\').')

    def run_pipeline(self, args):
        from clinica.pipeline.t1.t1_fsl import t1_fsl_segmentation_pipeline
        import pandas
        import os.path

        if not os.path.isfile(self.absolute_path(args.participants_sessions_tsv)):
            raise Exception('The TSV file does not exist.')
        subjects_visits = pandas.io.parsers.read_csv(self.absolute_path(args.participants_sessions_tsv), sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
            raise Exception('Subjects and visits file is not in the correct format.')
        subjects = list(subjects_visits.participant_id)
        sessions = list(subjects_visits.session_id)

#        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
        for index in xrange(len(subjects)):
            participant_id=subjects[index]
            session_id=sessions[index]
            bids_path_to_t1 = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'anat',
                                           participant_id + '_' + session_id + '_T1w.nii.gz')
            assert(os.path.isfile(bids_path_to_t1))
            t1_fsl_segmentation = t1_fsl_segmentation_pipeline(
                participant_id=participant_id,
                session_id=session_id,
                caps_directory=self.absolute_path(args.caps_directory),
                working_directory=self.absolute_path(args.working_directory),
                is_bias_corrected=args.is_bias_corrected
                )
            t1_fsl_segmentation.inputs.inputnode.in_t1 = bids_path_to_t1
            t1_fsl_segmentation.run('MultiProc', plugin_args={'n_procs': args.n_threads})



class CmdParserDWIPreprocessingPhaseDifferenceFieldmap(CmdParser):

    def define_name(self):
        self._name = 'dwi-preprocessing-phase-difference-fieldmap'


    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-n_threads", type=int, default=1,
                                help='Number of threads (default=1, which disables multi-threading).')
        registration = self._args.add_mutually_exclusive_group(required=True)
        registration.add_argument('-register_fmap_on_b0', action='store_true',
                           help='Choose to register fmap on b0. (mutually exclusive with \'-do_not_register_fmap_on_b0\').')
        registration.add_argument('-do_not_register_fmap_on_b0', action='store_false',
                           help='Choose not to register fmap on b0. (mutually exclusive with \'-register_fmap_on_b0\').')



    def run_pipeline(self, args):
        import os.path
        from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_phasediff_fieldmap
        from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
        import pandas

        subjects_visits = pandas.io.parsers.read_csv(self.absolute_path(args.subjects_sessions_tsv), sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id', 'dwi_effective_echo_spacing',
                                                    'delta_echo_time', 'dwi_phase_encoding_direction']:
            raise Exception(
                'The TSV file should contain the following columns: participant_id, session_id, dwi_effective_echo_spacing (in seconds), dwi_phase_encoding_direction (x/x-/y/y-/z/z-), delta_echo_time (in seconds).')
        subjects = list(subjects_visits.participant_id)
        sessions = list(subjects_visits.session_id)
        delta_echo_times = list(subjects_visits.delta_echo_time)
        echo_spacings = list(subjects_visits.dwi_effective_echo_spacing)
        phase_encoding_directions = list(subjects_visits.dwi_phase_encoding_direction)

        if args.register_fmap_on_b0:
            print("The fieldmap will be registered on the B0")
        else:
            print("The fieldmap will NOT be registered on the B0")


#        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
        for index in xrange(len(subjects)):
            participant_id = subjects[index]
            session_id = sessions[index]
            delta_echo_time = delta_echo_times[index]
            echo_spacing = echo_spacings[index]
            phase_encoding_direction = phase_encoding_directions[index]

            bids_path_to_dwi = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.nii.gz')
            bids_path_to_bval = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.bval')
            bids_path_to_bvec = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.bvec')
            bids_path_to_fmap_magnitude = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'fmap',
                                           participant_id + '_' + session_id + '_magnitude1.nii.gz')
            bids_path_to_fmap_phasediff = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'fmap',
                                           participant_id + '_' + session_id + '_phasediff.nii.gz')
            assert(os.path.isfile(bids_path_to_dwi))
            assert(os.path.isfile(bids_path_to_bval))
            assert(os.path.isfile(bids_path_to_bvec))
            assert(os.path.isfile(bids_path_to_fmap_magnitude))
            assert(os.path.isfile(bids_path_to_fmap_phasediff))
            preprocessing = diffusion_preprocessing_phasediff_fieldmap(
                participant_id=participant_id, session_id=session_id,
                caps_directory=self.absolute_path(args.caps_directory),
                delta_echo_time=delta_echo_time,
                effective_echo_spacing=echo_spacing,
                phase_encoding_direction=phase_encoding_direction,
                num_b0s=count_b0s(bids_path_to_bval),
                register_fmap_on_b0=args.register_fmap_on_b0,
                working_directory=self.absolute_path(args.working_directory)
            )
            preprocessing.inputs.inputnode.in_dwi = bids_path_to_dwi
            preprocessing.inputs.inputnode.in_bvals = bids_path_to_bval
            preprocessing.inputs.inputnode.in_bvecs = bids_path_to_bvec
            preprocessing.inputs.inputnode.in_fmap_magnitude = bids_path_to_fmap_magnitude
            preprocessing.inputs.inputnode.in_fmap_phasediff = bids_path_to_fmap_phasediff

            preprocessing.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserDWIPreprocessingTwoPhaseImagesFieldmap(CmdParser):
    def define_name(self):
        self._name = 'dwi-preprocessing-two-phase-images-fieldmap'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-n_threads", type=int, default=1,
                                help='Number of threads (default=1, which disables multi-threading).')

    def run_pipeline(self, args):
        import os.path
        from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_twophase_fieldmap
        from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
        import pandas

        subjects_visits = pandas.io.parsers.read_csv(self.absolute_path(args.subjects_sessions_tsv), sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id', 'dwi_effective_echo_spacing',
                                                    'delta_echo_time', 'dwi_phase_encoding_direction']:
            raise Exception(
                'The TSV file should contain the following columns: participant_id, session_id, dwi_effective_echo_spacing (in seconds), dwi_phase_encoding_direction (x/x-/y/y-/z/z-), delta_echo_time (in seconds).')
        subjects = list(subjects_visits.participant_id)
        sessions = list(subjects_visits.session_id)
        delta_echo_times = list(subjects_visits.delta_echo_time)
        echo_spacings = list(subjects_visits.dwi_effective_echo_spacing)
        phase_encoding_directions = list(subjects_visits.dwi_phase_encoding_direction)

        #        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
        for index in xrange(len(subjects)):
            participant_id = subjects[index]
            session_id = sessions[index]
            delta_echo_time = delta_echo_times[index]
            phase_encoding_direction = phase_encoding_directions[index]
            echo_spacing = echo_spacings[index]

            bids_path_to_dwi = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                            participant_id + '_' + session_id + '_dwi.nii.gz')
            bids_path_to_bval = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                             participant_id + '_' + session_id + '_dwi.bval')
            bids_path_to_bvec = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                             participant_id + '_' + session_id + '_dwi.bvec')
            bids_path_to_fmap_magnitude1 = os.path.join(self.absolute_path(args.bids_directory), participant_id,
                                                       session_id, 'fmap',
                                                       participant_id + '_' + session_id + '_magnitude1.nii.gz')
            bids_path_to_fmap_magnitude2 = os.path.join(self.absolute_path(args.bids_directory), participant_id,
                                                       session_id, 'fmap',
                                                       participant_id + '_' + session_id + '_magnitude2.nii.gz')
            bids_path_to_fmap_phase1 = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id,
                                                   'fmap', participant_id + '_' + session_id + '_phase1.nii.gz')
            bids_path_to_fmap_phase2 = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id,
                                                   'fmap', participant_id + '_' + session_id + '_phase2.nii.gz')
            assert (os.path.isfile(bids_path_to_dwi))
            assert (os.path.isfile(bids_path_to_bval))
            assert (os.path.isfile(bids_path_to_bvec))
            assert (os.path.isfile(bids_path_to_fmap_magnitude1))
            assert (os.path.isfile(bids_path_to_fmap_magnitude2))
            assert (os.path.isfile(bids_path_to_fmap_phase1))
            assert (os.path.isfile(bids_path_to_fmap_phase2))
            preprocessing = diffusion_preprocessing_twophase_fieldmap(
                participant_id=participant_id, session_id=session_id,
                caps_directory=self.absolute_path(args.caps_directory),
                delta_echo_time=delta_echo_time,
                effective_echo_spacing=echo_spacing,
                phase_encoding_direction=phase_encoding_direction,
                num_b0s=count_b0s(bids_path_to_bval),
                working_directory=self.absolute_path(args.working_directory)
            )
            preprocessing.inputs.inputnode.in_dwi = bids_path_to_dwi
            preprocessing.inputs.inputnode.in_bvals = bids_path_to_bval
            preprocessing.inputs.inputnode.in_bvecs = bids_path_to_bvec
            preprocessing.inputs.inputnode.in_fmap_magnitude1 = bids_path_to_fmap_magnitude1
            preprocessing.inputs.inputnode.in_fmap_magnitude2 = bids_path_to_fmap_magnitude2
            preprocessing.inputs.inputnode.in_fmap_phase1 = bids_path_to_fmap_phase1
            preprocessing.inputs.inputnode.in_fmap_phase2 = bids_path_to_fmap_phase2

            preprocessing.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserDWIPreprocessingT1Based(CmdParser):

    def define_name(self):
        self._name = 'dwi-preprocessing-t1-based'


    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-n_threads", type=int, default=1,
                                help='Number of threads (default=1, which disables multi-threading).')

    def run_pipeline(self, args):
        import os.path
        from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_t1_based
        from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
        import pandas

        subjects_visits = pandas.io.parsers.read_csv(self.absolute_path(args.subjects_sessions_tsv), sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
            raise Exception('Subjects and visits file is not in the correct format.')
        subjects = list(subjects_visits.participant_id)
        sessions = list(subjects_visits.session_id)

#        with open(self.absolute_path(args.subjects_sessions_tsv), 'rb') as tsv_file:
        for index in xrange(len(subjects)):
            participant_id=subjects[index]
            session_id=sessions[index]

            bids_path_to_t1 = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'anat',
                                           participant_id + '_' + session_id + '_T1w.nii.gz')
            bids_path_to_dwi = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.nii.gz')
            bids_path_to_bval = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.bval')
            bids_path_to_bvec = os.path.join(self.absolute_path(args.bids_directory), participant_id, session_id, 'dwi',
                                           participant_id + '_' + session_id + '_dwi.bvec')
            assert(os.path.isfile(bids_path_to_t1))
            assert(os.path.isfile(bids_path_to_dwi))
            assert(os.path.isfile(bids_path_to_bval))
            assert(os.path.isfile(bids_path_to_bvec))
            preprocessing = diffusion_preprocessing_t1_based(
                participant_id=participant_id, session_id=session_id,
                caps_directory=self.absolute_path(args.caps_directory),
                num_b0s=count_b0s(bids_path_to_bval),
                working_directory=self.absolute_path(args.working_directory)
            )
            preprocessing.inputs.inputnode.in_t1= bids_path_to_t1
            preprocessing.inputs.inputnode.in_dwi = bids_path_to_dwi
            preprocessing.inputs.inputnode.in_bvals = bids_path_to_bval
            preprocessing.inputs.inputnode.in_bvecs = bids_path_to_bvec

            preprocessing.run('MultiProc', plugin_args={'n_procs': args.n_threads})


class CmdParserDWIProcessing(CmdParser):

    def define_name(self):
        self._name = 'dwi-processing'


    def define_options(self):
        self._args.add_argument("caps_directory",
                                help='Path to the output/input directory in a CAPS format.')
        self._args.add_argument("participants_sessions_tsv",
                                help='TSV file containing the subjects/sessions list to be processed.')
        self._args.add_argument("-working_directory", default=None,
                                help='Temporary directory to store intermediate results')
        self._args.add_argument("-n_threads", type=int, default=1,
                                help='Number of threads (default=1, which disables multi-threading).')
        self._args.add_argument("number_of_tracks",
                                help="Set the desired number of tracks (e.g. 100K, 1M, ...).")
        self._args.add_argument("-max_harmonic_order", default=None,
                                help="@TODO")
        self._args.add_argument("-tractography_algorithm", default='iFOD2',
                                help="@TODO")
        self._args.add_argument("-tractography_fod_threshold", default=None,
                                help="@TODO")
        self._args.add_argument("-tractography_step_size", default=None,
                                help="@TODO")
        self._args.add_argument("-tractography_angle", default=None,
                                help="@TODO")


    def run_pipeline(self, args):
        import pandas
        import os.path
        from clinica.pipeline.dwi.dwi_processing import dwi_processing_pipeline

        if not os.path.isfile(self.absolute_path(args.participants_sessions_tsv)):
            raise Exception('The TSV file does not exist.')
        subjects_visits = pandas.io.parsers.read_csv(self.absolute_path(args.participants_sessions_tsv), sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
            raise Exception('Subjects and visits file is not in the correct format.')
        subjects = list(subjects_visits.participant_id)
        sessions = list(subjects_visits.session_id)

        for index in xrange(len(subjects)):
            participant_id = subjects[index]
            session_id = sessions[index]

            # DWI-Preprocessing:
            caps_path_to_preprocessed_dwi = os.path.join(self.absolute_path(args.caps_directory), 'subjects',
                                                         participant_id, session_id, 'dwi', 'preprocessing',
                                                         participant_id + '_' + session_id + '_dwi.nii.gz')
            caps_path_to_bvecs = os.path.join(self.absolute_path(args.caps_directory), 'subjects', participant_id,
                                              session_id, 'dwi', 'preprocessing',
                                              participant_id + '_' + session_id + '_dwi.bvec')
            caps_path_to_bvals = os.path.join(self.absolute_path(args.caps_directory), 'subjects', participant_id,
                                              session_id, 'dwi', 'preprocessing',
                                              participant_id + '_' + session_id + '_dwi.bval')
            caps_path_to_b0_mask = os.path.join(self.absolute_path(args.caps_directory), 'subjects', participant_id,
                                                session_id, 'dwi', 'preprocessing',
                                                participant_id + '_' + session_id + '_b0Mask.nii.gz')
            # T1-FSL:
            caps_path_to_bias_corrected_bet_t1 = os.path.join(self.absolute_path(args.caps_directory), 'subjects',
                                                              participant_id, session_id, 't1', 'fsl',
                                                              participant_id + '_' + session_id + '_brainExtractedT1w.nii.gz')
            caps_path_to_white_matter_binary_mask = os.path.join(self.absolute_path(args.caps_directory), 'subjects',
                                                                 participant_id, session_id, 't1', 'fsl',
                                                                 participant_id + '_' + session_id + '_tissue-whitematter_binaryMask.nii.gz')
            # T1-FreeSurfer:
            caps_path_to_desikan_parcellation = os.path.join(self.absolute_path(args.caps_directory), 'subjects',
                                                             participant_id, session_id, 't1',
                                                             'freesurfer-cross-sectional',
                                                             participant_id + '_' + session_id, 'mri', 'aparc+aseg.mgz')
            print(caps_path_to_desikan_parcellation)
            caps_path_to_destrieux_parcellation = os.path.join(self.absolute_path(args.caps_directory), 'subjects',
                                                               participant_id, session_id, 't1',
                                                               'freesurfer-cross-sectional',
                                                               participant_id + '_' + session_id, 'mri',
                                                               'aparc.a2009s+aseg.mgz')
            print(caps_path_to_destrieux_parcellation)

            assert (os.path.isfile(caps_path_to_preprocessed_dwi))
            assert (os.path.isfile(caps_path_to_bvecs))
            assert (os.path.isfile(caps_path_to_bvals))
            assert (os.path.isfile(caps_path_to_bias_corrected_bet_t1))
            assert (os.path.isfile(caps_path_to_b0_mask))
            assert (os.path.isfile(caps_path_to_white_matter_binary_mask))
            assert (os.path.isfile(caps_path_to_desikan_parcellation))
            assert (os.path.isfile(caps_path_to_destrieux_parcellation))

            dwi_processing = dwi_processing_pipeline(
                participant_id=participant_id,
                session_id=session_id,
                caps_directory=self.absolute_path(args.caps_directory),
                working_directory=self.absolute_path(args.working_directory),
                atlas_name='JHU-ICBM-tracts-maxprob-thr25',
                max_harmonic_order=args.max_harmonic_order,
                tractography_algorithm=args.tractography_algorithm,
                tractography_nb_of_tracks=args.number_of_tracks,
                tractography_fod_threshold=args.tractography_fod_threshold,
                tractography_step_size=args.tractography_step_size,
                tractography_angle=args.tractography_angle,
                nthreads=args.n_threads,
                zero_diagonal=True
            )
            dwi_processing.inputs.inputnode.in_preprocessed_dwi = caps_path_to_preprocessed_dwi
            dwi_processing.inputs.inputnode.in_bvecs = caps_path_to_bvecs
            dwi_processing.inputs.inputnode.in_bvals = caps_path_to_bvals
            dwi_processing.inputs.inputnode.in_bias_corrected_bet_t1 = caps_path_to_bias_corrected_bet_t1
            dwi_processing.inputs.inputnode.in_b0_mask = caps_path_to_b0_mask
            dwi_processing.inputs.inputnode.in_white_matter_binary_mask = caps_path_to_white_matter_binary_mask
            dwi_processing.inputs.inputnode.in_desikan_parcellation = caps_path_to_desikan_parcellation
            dwi_processing.inputs.inputnode.in_destrieux_parcellation = caps_path_to_destrieux_parcellation

            dwi_processing.run('MultiProc', plugin_args={'n_procs': args.n_threads})




# class CmdParserInsightToBids(CmdParser):
#
#     def define_name(self):
#         self._name = 'insight-to-bids'
#
#     def define_options(self):
#         self._args.add_argument("dataset_directory",
#                                help='Path of the unorganized INSIGHT directory.')
#         self._args.add_argument("bids_directory",
#                                 help='Path to the BIDS directory.')
#         self._args.add_argument("-co", type=bool, default=False,
#                                 help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')
#
#     def run_pipeline(self, args):
#         from clinica.bids import insight_to_bids
#         insight_to_bids.convert(args.dataset_directory, args.bids_directory)


# class CmdParserPrevDemAlsToBids(CmdParser):
#
#     def define_name(self):
#         self._name = 'prevdemals-to-bids'
#
#     def define_options(self):
#         self._args.add_argument("dataset_directory",
#                                help='Path of the unorganized INSIGHT directory.')
#         self._args.add_argument("bids_directory",
#                                 help='Path to the BIDS directory.')
#         self._args.add_argument("-co", type=bool, default=False,
#                                 help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')
#
#     def run_pipeline(self, args):
#         from clinica.bids import prevdemals_to_bids
#         prevdemals_to_bids.convert(args.dataset_directory, args.bids_directory)


class CmdParserHmtcToBids(CmdParser):

    def define_name(self):
        self._name = 'hmtc-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized HMTC directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_pipeline(self, args):
        from clinica.iotools import hmtc_to_bids
        hmtc_to_bids.convert(args.dataset_directory, args.bids_directory)


class CmdParserSubsSess(CmdParser):

    def define_name(self):
        self._name = 'create-subjects-visits'

    def define_options(self):
        self._args.add_argument("bids_dir",
                               help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')
        self._args.add_argument("-on", '--output_name', type=str, default='',
                            help='(Optional) Name of the output file.')

    def run_pipeline(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.create_subs_sess_list(args.bids_dir, args.out_directory, args.output_name)

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
        from clinica.iotools.utils import data_handling as dt
        dt.create_merge_file(args.bids_dir, args.out_directory, args.true_false_mode)

class CmdParserMissingModalities(CmdParser):

    def define_name(self):
        self._name = 'check-missing-mods'

    def define_options(self):
        self._args.add_argument("bids_dir",
                               help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')
        self._args.add_argument("-op", '--output_prefix', type= str, default= '',
                                help='Prefix for the name of output files.')



    def run_pipeline(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.compute_missing_mods(args.bids_dir, args.out_directory, args.output_prefix)
