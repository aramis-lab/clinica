# coding: utf8

import clinica.engine as ce


class DWIProcessingDTICLI(ce.CmdParser):

    def __init__(self):
        super(DWIProcessingDTICLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        from colorama import Fore, init
        init()

        self._name = 'dwi-processing-dti'
        self._args.epilog = 'Example: clinica run dwi-processing-dti BIDS CAPS'
        self._args.description = "%sDTI-based pipeline: http:://clinica.run/doc/DWIProcessing%s" % (Fore.GREEN, Fore.RESET)
        self._args._positionals.title = '%sCompulsory arguments:%s' % (Fore.BLUE, Fore.RESET)
        self._args._optionals.title = '%sOptional arguments:%s' % (Fore.BLUE, Fore.RESET)

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')  # noqa

        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipelines intermediate results')  # noqa
        self._args.add_argument("-np", "--n_procs",
                                type=int,
                                help='Number of cores used to run in parallel')  # noqa
        self._args.add_argument("-sl", "--slurm",
                                action='store_true',
                                help='Run the pipelines using SLURM')  # noqa

    def run_pipeline(self, args):
        """
        """
        from tempfile import mkdtemp
        from dwi_processing_dti_pipeline import DWIProcessingDTI

        pipeline = DWIProcessingDTI(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
        )
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
