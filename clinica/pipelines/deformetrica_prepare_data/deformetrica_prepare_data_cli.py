# -*- coding: utf-8 -*-

import clinica.engine as ce

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"


class deformetrica_prepare_dataCLI(ce.CmdParser):
    """CLI class
    """

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'deformetrica-prepare-data'

    def define_description(self):
        """Provide a description of this pipeline.
        """
        self._description = (
            'Generation of subcortical meshes input for Deformetrica:\n'
            'http://clinica.run/doc/Pipelines/Deformetrica_prepare_data/')

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument(
            "group_id",
            help='User-defined identifier for the provided group of subjects.')
        # Optional arguments
        optional = self._args.add_argument_group(
            PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-st", "--structure",
                              default=None,
                              help='Structure name.')
        optional.add_argument("-stf", "--structure_file",
                              default=None,
                              help='Path to structure file.')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(
            PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument(
            "-tsv", "--subjects_sessions_tsv",
            help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument(
            "-wd", "--working_directory",
            help='Temporary directory where the pipeline’s intermediate results are stored.')
        clinica_opt.add_argument(
            "-np", "--n_procs",
            metavar=('N'), type=int,
            help='Number of cores used to run in parallel')

    def run_command(self, args):
        """
        Run the pipelines with defined args
        """

        from tempfile import mkdtemp
        from ./deformetrica_prepare_data_pipeline import DeformetricaPrepareData

        pipeline = DeformetricaPrepareData(
            # pass these args by the class attribute itself
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            group_id=args.group_id)
        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'structure': args.structure or '-qcache',
            'structure_file': args.structure_file or '-qcache',
        }
        # make sure if working_directory is not defined
        # using a temp folder as the working directory.
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        # run the pipelines with [n_procs] cores
        if args.n_procs:
            pipeline.run(
                plugin='MultiProc',
                plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
