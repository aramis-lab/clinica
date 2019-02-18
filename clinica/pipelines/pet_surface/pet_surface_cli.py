# coding: utf-8

import clinica.engine as ce


__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
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
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory. (Filled with results from t1-freesurfer pipeline')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-pt", "--pet_tracer", type=str, default='fdg',
                              help='PET tracer type. Can be fdg or av45 (default: --pet_tracer fdg)')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores to use when running the pipeline in parallel (default: --n_procs 1).')

    def run_command(self, args):
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

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
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
