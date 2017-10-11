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
        from os.path import join, split, realpath

        output_directory = join(self.absolute_path(args.caps_directory),
                                            'group-' +args.group_id + '/machine_learning/region_based_svm/',
                                            'space' + args.atlas_id, args.image_type)

        if args.subjects_visits_tsv is None:
            subjects_visits_tsv = () # TODO where it's saved for t1 and pet
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
        input_image_atlas = join(split(realpath(__file__))[0], '../resources/atlases_spm', args.atlas_id + '.nii')
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
