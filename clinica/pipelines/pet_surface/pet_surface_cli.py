# coding: utf-8

import clinica.engine as ce


__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2018 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux", "Michael Bacci"]
__license__ = "See LICENSE.txt file"
__version__ = "1.0.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@inria.fr"
__status__ = "Development"


class PetSurfaceCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'pet-surface'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Surface-based processing of PET images:\nhttp://clinica.run/doc/Pipelines/PET_Surface/'

    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory. (Filled with results from t1-freesurfer-cross-sectional')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing a list of subjects with their sessions.')
        self._args.add_argument("-pt", "--pet_tracer", type=str, default='fdg',
                                help='PET tracer type. Can be fdg or av45. Default value : fdg')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')

    def run_command(self, args):

        from tempfile import mkdtemp
        from pet_surface_pipeline import PetSurface
        if args.working_directory is None:
            args.working_directory = mkdtemp()

        pipeline = PetSurface(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'pet_type': args.pet_tracer,
            'wd': self.absolute_path(args.working_directory),
            'n_procs': args.n_procs
        }
        pipeline.base_dir = self.absolute_path(mkdtemp())
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
