# coding: utf8

"""Tractography - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""

import clinica.engine as ce


class TractographyCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'tractography'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Run tractography algorithm using preprocessed ' \
                            'DWI ' \
                            'images:\nhttp://clinica.run/doc/Pipelines' \
                            '/Tractography/'

    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with '
                                     'their sessions.')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline '
                                     'intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        self._args.add_argument("-nt", "--n_tracks", type=int,
                                help='Number of tracts gene')
        self._args.add_argument("-sl", "--slurm", action='store_true',
                                help='Run the pipeline using SLURM')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from tractography_pipeline import Tractography

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = Tractography(
                bids_directory=self.absolute_path(args.bids_directory),
                caps_directory=self.absolute_path(args.caps_directory),
                tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'n_tracks'        : args.n_tracks or 100000,
        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        pipeline.write_graph()
        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
