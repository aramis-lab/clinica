# coding: utf8

# coding: utf-8

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.3.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@icm-institute.org"
__status__ = "Development"

import clinica.engine as ce


class StatisticsVolumeCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'statistics-volume'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Brief description:\n'
                             'http://clinica.run/doc/Pipelines/Statistics_Volume/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])

        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        clinica_comp.add_argument("subject_visits_with_covariates_tsv",
                                  help='TSV file containing a list of subjects with their sessions and all '
                                       'the covariates and factors needed for the GLM.')

        clinica_comp.add_argument("contrast",
                                  help='Defines the contrast. Must be one of the column names form the TSV file.')

        clinica_comp.add_argument("feature_type",
                                  help='Define what type of file are grabbed for the analysis. Use "fdg" or "graymatter". Use your own if you want to use the --custom_file flag. ')

        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')

        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])

        optional.add_argument("--custom_files", "-cf", type=str, default=None,
                              help=('Custom file string. Specify filename using * when the subject or session name appear. '
                                    + 'Example : \'*_task-rest_acq-fdg_pet_space-Ixi549Space_pet.nii.gz\' will grab the'
                                    + ' corresponding file in all the subjects/sessions'))

        optional.add_argument("-tup", "--threshold_uncorrected_pvalue",
                              type=float, default=0.001,
                              help='Threshold to display the uncorrected p-value '
                                   + '(--threshold_uncorrected_pvalue 0.001).')

        optional.add_argument("-tcp", "--threshold_corrected_pvalue",
                              type=float, default=0.05,
                              help='Threshold to display the corrected p-value '
                                   '(default: --threshold_corrected_pvalue 0.05)')

        optional.add_argument("-gic", "--group_id_caps", type=str, default=None,
                              help='Name of the group that Clinica needs to use to grab input file')

        optional.add_argument("-fwhm", "--smoothing", type=int, default=8,
                              help='Full Width at Half Maximum (FWHM) of the smoothing used in your input file.')

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        import os
        from networkx import Graph
        from colorama import Fore
        from .statistics_volume_pipeline import StatisticsVolume
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and/or CAPS directory as inputs. If the BIDS directory is not needed
        # for your pipeline, simply remove:
        # bids_directory=self.absolute_path(args.bids_directory),
        pipeline = StatisticsVolume(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subject_visits_with_covariates_tsv),
            base_dir=self.absolute_path(args.working_directory),
            name=self.name
        )
        pipeline.parameters = {'contrast': args.contrast,
                               'feature_type': args.feature_type,
                               'group_id': args.group_id,
                               'custom_files': args.custom_files,
                               'threshold_uncorrected_pvalue': args.threshold_uncorrected_pvalue,
                               'threshold_corrected_pvalue': args.threshold_corrected_pvalue,
                               'group_id_caps': args.group_id_caps,
                               'smoothing': args.smoothing}

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, os.path.join(pipeline.base_dir, self.name), pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
