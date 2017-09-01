# coding: utf8

import clinica.engine as ce


class DWIProcessingCLI(ce.CmdParser):

    def __init__(self):
        super(DWIProcessingCLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-processing'

    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')  # noqa

        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')  # noqa
        self._args.add_argument("-np", "--n_procs",
                                type=int,
                                help='Number of cores used to run in parallel')
        self._args.add_argument("-sl", "--slurm",
                                action='store_true',
                                help='Run the pipeline using SLURM')

    def run_pipeline(self, args):
        """
        """
        from tempfile import mkdtemp
        from dwi_processing_pipeline import DWIProcessing

        from clinica.utils.stream import cprint

        cprint('TSV file = ' + args.subjects_sessions_tsv + '.')

        pipeline = DWIProcessing(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
        )
#        pipeline.parameters = {
#            # Add your own pipeline parameters here to use them inside your
#            # pipeline. See the file `dwi_processing_pipeline.py` to
#            # see an example of use.
#            'hello_word'        : args.hello_word_arg or 'Hello'
#        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
