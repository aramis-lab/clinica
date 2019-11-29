# coding: utf8

# coding: utf-8

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.3.0"
__maintainer__ = "Alexandre Routier"
__email__ = "alexandre.routier@icm-institute.org"
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

        clinica_comp.add_argument("file_id",
                                  help='Define what type of file are grabbed for the analysis')

        clinica_comp.add_argument("covariables",
                                  help='Defines subject list and list of associated covariables')

        clinica_comp.add_argument("contrast",
                                  help='Defines the contrast')

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
            tsv_file=self.absolute_path(args.covariables),
            base_dir=self.absolute_path(args.working_directory),
            name=self.name
        )
        pipeline.parameters = {'contrast': args.contrast,
                               'file_id': args.file_id}

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, os.path.join(pipeline.base_dir, self.name), pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)