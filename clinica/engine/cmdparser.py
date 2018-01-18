# -*- coding: utf-8 -*-

"""Define the command line parser for each utility

Please, add your command line parser at the end of this file.

Some examples can be found at the end of this file.
"""

import abc
from argparse import ArgumentParser
from os.path import join
from os import getcwd
from os.path import expanduser


class CmdParser:
    """Abstract class to extend in order to create your command line parser

    For pipelines, please the 'clinica generate template' command.
    For utilities, see examples at the end of this file.
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
    """
    Init all derived CmdParser instances with specific data.

    Args:
        rootparser: The root parser
        parser: The ArgParser node (e.g. 'run' or 'convert')
        objects: All CmdParser instances of this file
    """

    def silent_help():
        pass

    def error_message(p):
        def error(x):
            p.print_help()
            rootparser.print_help = silent_help
            exit(-1)
        return error

    def init(x):
        x.options = parser.add_parser(x.name)
        x.options.error = error_message(x.options)
        x.options.set_defaults(func=x.run_pipeline)
        x.build()

    for x in objects:
        try:
            init(x)
        except:
            pass


def get_cmdparser_names(objects=None):
    """
    Return the names of all pipelines

    Args:
        objects: All CmdParser instances of this file

    Returns:
        The names of all pipelines
    """
    if objects is None:
        objects = get_cmdparser_objects()
    for x in objects:
        yield x.name


#
# Insert here pipelines which do not fit in standard pipelines
#

class CmdParserMachineLearningVBLinearSVM(CmdParser):

    def define_name(self):
        self._name = 'machinelearning-svm-voxel'

    def define_options(self):
        self._args.add_argument("image_type",
                                help='it can assume two values: pet/t1, according to the images used')  # noqa
        self._args.add_argument("caps_directory",
                                help='Directory where the input NIFTI images are stored')  # noqa
        self._args.add_argument("group_id",
                                help='Current group name')  # noqa
        self._args.add_argument("diagnoses_tsv",
                                help='TSV file with subjects diagnoses')  # noqa
        self._args.add_argument("-p", "--prefix",
                                help='Images prefix')  # noqa
        self._args.add_argument("-t", "--tissue",
                                help='')
        self._args.add_argument("-svt", "--subjects_visits_tsv",
                                type=str, default=None,
                                help='TSV file with subjects and sessions to be processed')  # noqa
        self._args.add_argument("-mz", "--mask_zeros",
                                type=bool, default=True,
                                help='Use a mask to remove zero valued voxels across images')  # noqa
        self._args.add_argument("-b", "--balanced",
                                type=bool, default=True,
                                help='Balance the weights of subjects for the SVM proportionally to their number in each class')  # noqa
        self._args.add_argument("-cv", "--cv_folds",
                                type=int, default=10,
                                help='Number of folds to use in the cross validation')  # noqa
        self._args.add_argument("-fc", "--folds_c",
                                type=int, default=10,
                                help='Number of folds to use in the cross validation to determine parameter C')   # noqa
        self._args.add_argument("-np", "--n_procs",
                                type=int, default=4,
                                help='Number of parallel processes to run')  # noqa
        self._args.add_argument("-crl", "--c_range_logspace",
                                nargs=3, type=int, default=[-6, 2, 17],
                                help="numpy logspace function arguments defining the range of search for SVM parameter C. Ex: -6 2 17")  # noqa
        self._args.add_argument("-sgm", "--save_gram_matrix",
                                action='store_true',
                                help="Save feature weights for each classification as a matrix")  # noqa
        self._args.add_argument("-sc", "--save_subject_classification",
                                action='store_true',
                                help="Save list of classification results for each subject for each classification")  # noqa
        self._args.add_argument("-sw", "--save_original_weights",
                                action='store_true',
                                help="Save feature weights for each classification as a matrix")  # noqa
        self._args.add_argument("-sf", "--save_features_image",
                                action='store_true',
                                help="Save feature weights for each classification as an image")  # noqa
        self._args.add_argument("-sdc", "--save_dual_coefficients",
                                action='store_true',
                                help="Save ")  # noqa

    def run_pipeline(self, args):

        from clinica.pipelines.machine_learning.voxel_based_svm import linear_svm_binary_classification_caps
        from numpy import logspace

        if args.subjects_visits_tsv is None:
            subjects_visits_tsv = ()  # where it's saved for t1 and pet
        else:
            subjects_visits_tsv = pandas.io.parsers.read_csv(self.absolute_path(args.participants_sessions_tsv),
                                                             sep='\t')

        c_range = logspace(args.c_range_logspace[0], args.c_range_logspace[1], args.c_range_logspace[2])

        linear_svm_binary_classification_caps(
            self.absolute_path(args.caps_directory),
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
            save_features_image=args.save_features_image
        )


class CmdParserMachineLearningSVMRB(CmdParser):

    def define_name(self):
        self._name = 'machinelearning-svm-region'

    def define_options(self):
        self._args.add_argument("image_type",
                                help='it can assume two values: pet/t1, according to the images used')  # noqa
        self._args.add_argument("caps_directory",
                                help='Directory where the input NIFTI images are stored')  # noqa
        self._args.add_argument("group_id",
                                help='Current group name')  # noqa
        self._args.add_argument("diagnosis_tsv",
                                help='TSV file with subjects diagnosis')  # noqa
        self._args.add_argument("atlas_id",
                                help='Name of the atlas used to extract features')  # noqa
        self._args.add_argument("-svt", "--subjects_visits_tsv",
                                type=str, default=None,
                                help='TSV file with subjects and sessions to be processed')  # noqa
        self._args.add_argument("-b", "--balanced",
                                type=bool, default=True,
                                help='Balance the weights of subjects for the SVM proportionally to their number in each class')  # noqa
        self._args.add_argument("-cv", "--cv_folds",
                                type=int, default=10,
                                help='Number of folds to use in the cross validation')  # noqa
        self._args.add_argument("-fc", "--folds_c",
                                type=int, default=10,
                                help='Number of folds to use in the cross validation to determine parameter C')  # noqa
        self._args.add_argument("-np", "--n_procs",
                                type=int, default=4,
                                help='Number of parallel processes to run')
        self._args.add_argument("-crl", "--c_range_logspace",
                                nargs=3, type=int, default=[-6, 2, 17],
                                help="numpy logspace function arguments defining the range of search for SVM parameter C. Ex: -6 2 17")  # noqa
        self._args.add_argument("-sgm", "--save_gram_matrix",
                                action='store_true',
                                help="Save gram matrix for each classification as a matrix")  # noqa
        self._args.add_argument("-sdc", "--save_dual_coefficients",
                                action='store_true',
                                help="Save ")
        self._args.add_argument("-sc", "--save_subject_classification",
                                action='store_true',
                                help="Save list of classification results for each subject for each classification")  # noqa
        self._args.add_argument("-sw", "--save_original_weights",
                                action='store_true',
                                help="Save feature weights for each classification as a matrix")  # noqa
        self._args.add_argument("-sf", "--save_features_image",
                                action='store_true',
                                help="Save feature weights for each classification as an image")  # noqa

    def run_pipeline(self, args):
        from clinica.pipelines.machine_learning.region_based_svm import svm_binary_classification
        from clinica.pipelines.machine_learning.region_based_io import get_caps_pet_list, get_caps_t1_list, load_data
        from clinica.pipelines.machine_learning.svm_utils import gram_matrix_linear
        from numpy import logspace
        import pandas
        from os.path import join, split, realpath

        output_directory = join(
            self.absolute_path(args.caps_directory),
            'group-' + args.group_id + '/machine_learning/region_based_svm/',
            'space' + args.atlas_id,
            args.image_type)

        if args.subjects_visits_tsv is None:
            subjects_visits_tsv = ()  # TODO where it's saved for t1 and pet
        else:
            subjects_visits_tsv = pandas.io.parsers.read_csv(
                self.absolute_path(args.participants_sessions_tsv),
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
        c_range = logspace(args.c_range_logspace[0],
                           args.c_range_logspace[1],
                           args.c_range_logspace[2])

        svm_binary_classification(
            input_image_atlas,
            image_list,
            diagnosis_list,
            output_directory,
            kernel_function=None, existing_gram_matrix=gram_matrix,
            mask_zeros=True,
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
            save_features_image=args.save_features_image
        )

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


# class CmdParserHmtcToBids(CmdParser):
#
#     def define_name(self):
#         self._name = 'hmtc-to-bids'
#
#     def define_options(self):
#         self._args.add_argument("dataset_directory",
#                                 help='Path of the unorganized HMTC directory.')
#         self._args.add_argument("bids_directory",
#                                 help='Path to the BIDS directory.')
#
#     def run_pipeline(self, args):
#         from clinica.iotools import hmtc_to_bids
#         hmtc_to_bids.convert(args.dataset_directory, args.bids_directory)


#
# Insert here iotools parsers
#

class CmdParserSubsSess(CmdParser):

    def define_name(self):
        self._name = 'create-subjects-visits'

    def define_options(self):
        self._args.add_argument("bids_dir",
                                help='Path to the BIDS dataset directory.')  # noqa
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')  # noqa
        self._args.add_argument("-on", '--output_name',
                                type=str, default='',
                                help='(Optional) Name of the output file.')  # noqa

    def run_pipeline(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.create_subs_sess_list(
            args.bids_dir,
            args.out_directory,
            args.output_name
        )


class CmdParserMergeTsv(CmdParser):

    def define_name(self):
        self._name = 'merge-tsv'

    def define_options(self):
        self._args.add_argument("bids_dir",
                                help='Path to the BIDS dataset directory.')  # noqa
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')  # noqa
        self._args.add_argument("-tf", '--true_false_mode',
                                type=bool, default=False,
                                help='(Optional) Convert all the field with binary meaning into True and False values.')  # noqa

    def run_pipeline(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.create_merge_file(
            args.bids_dir,
            args.out_directory,
            args.true_false_mode
        )


class CmdParserMissingModalities(CmdParser):

    def define_name(self):
        self._name = 'check-missing-modalities'

    def define_options(self):
        self._args.add_argument("bids_dir",
                                help='Path to the BIDS dataset directory.')  # noqa
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')  # noqa
        self._args.add_argument("-op", '--output_prefix',
                                type=str, default='',
                                help='Prefix for the name of output files.')  # noqa

    def run_pipeline(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.compute_missing_mods(
            args.bids_dir,
            args.out_directory,
            args.output_prefix
        )
