# coding: utf8

import clinica.engine as ce


class T1VolumeCreateDartelCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 't1-volume-create-dartel'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Inter-subject registration using Dartel (creating a new Dartel template):\n'
                             'http://clinica.run/doc/Pipelines/T1_Volume/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-dt", "--dartel_tissues",
                              metavar='', nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help='Tissues to use for DARTEL template calculation '
                                   '(default: GM, WM and CSF i.e. --dartel_tissues 1 2 3).')

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .t1_volume_create_dartel_pipeline import T1VolumeCreateDartel
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'group_id': args.group_id,
            'dartel_tissues': args.dartel_tissues
        }
        pipeline = T1VolumeCreateDartel(
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
