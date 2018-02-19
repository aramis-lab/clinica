# coding: utf8

import clinica.engine as ce


class DWIPreprocessingUsingT1CLI(ce.CmdParser):

    def __init__(self):
        super(DWIPreprocessingUsingT1CLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-preprocessing-using-t1'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Preprocessing of raw DWI datasets using a T1w image:\nhttp://clinica.run/doc/Pipelines/DWI_Preprocessing/'

    def define_options(self):
        """Define the sub-command arguments.
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')  # noqa
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')  # noqa
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')  # noqa

        self._args.add_argument("--low_bval",
                                type=int, default=5,
                                help='Define the b0 volumes as all volume bval <= lowbval. (Default: --low_bval 5)')  # noqa

        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')  # noqa
        self._args.add_argument("-np", "--n_procs",
                                type=int,
                                help='Number of cores used to run in parallel')  # noqa
        self._args.add_argument("-sl", "--slurm",
                                action='store_true',
                                help='Run the pipeline using SLURM')  # noqa

    def run_command(self, args):
        """
        Run the DWIPreprocessingUsingT1 Pipeline from command line.
        """
        from tempfile import mkdtemp
        from dwi_preprocessing_using_t1_pipeline import DWIPreprocessingUsingT1

        pipeline = DWIPreprocessingUsingT1(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            low_bval=args.low_bval
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
