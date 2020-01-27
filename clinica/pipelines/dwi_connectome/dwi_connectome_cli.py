# coding: utf8

import clinica.engine as ce


class DwiConnectomeCli(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'dwi-connectome'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Connectome-based processing of DWI datasets:\n'
                             'http://clinica.run/doc/DWI_Connectome')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-nt", "--n_tracks",
                              metavar='N', type=int, default=1000000,
                              help=('Set the desired number of streamlines to generate the tractography and connectome '
                                    '(default: --n_tracks %(default)s).'))
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .dwi_connectome_pipeline import DwiConnectome
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'n_tracks': args.n_tracks
        }
        pipeline = DwiConnectome(
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
