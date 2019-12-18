# coding: utf8

import clinica.engine as ce


class T1VolumeParcellationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 't1-volume-parcellation'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Computation of mean GM concentration for a set of regions:\n'
                             'http://clinica.run/doc/Pipelines/T1_Volume/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        default_atlases = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
        advanced.add_argument("-atlases", "--atlases",
                              nargs='+', type=str, metavar='', default=default_atlases, choices=default_atlases,
                              help='A list of atlases used to calculate the regional mean GM concentrations (default: '
                                   'all atlases i.e. --atlases %s).' % self.list_to_string(default_atlases))

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .t1_volume_parcellation_pipeline import T1VolumeParcellation
        from clinica.utils.check_dependency import verify_cat12_atlases
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        # If the user wants to use any of the atlases of cat12 and has not installed it, we just remove it from the list
        # of the computed atlases
        args.atlases = verify_cat12_atlases(args.atlases)

        parameters = {
            'group_id': args.group_id,
            'atlases': args.atlases,
        }
        pipeline = T1VolumeParcellation(
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
