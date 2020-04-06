# coding: utf-8

import clinica.engine as ce


class PetSurfaceLongitudinalCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'pet-surface-longitudinal'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Longitudinal surface-based processing of PET images:\n'
                             'http://clinica.run/doc/Pipelines/PET_Surface_Longitudinal/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory. (Filled with results from t1-freesurfer-longitudinal pipeline')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-pt", "--pet_tracer", type=str, default='fdg',
                              help='PET tracer type. Can be fdg or av45 (default: --pet_tracer %(default)s)')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .pet_surface_pipeline import PetSurface
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'pet_tracer': args.pet_tracer,
            'longitudinal': True
        }
        pipeline = PetSurface(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name
        )

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
