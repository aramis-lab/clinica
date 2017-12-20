# coding: utf8

import clinica.engine as ce


class DWIPreprocessingUsingPhaseDiffFieldmapCLI(ce.CmdParser):

    def __init__(self):
        super(DWIPreprocessingUsingPhaseDiffFieldmapCLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-preprocessing-using-phasediff-fieldmap'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')  # noqa

        self._args.add_argument("--low_bval",
                                type=int, default=5,
                                help='Define the b0 volumes as all volume bval <= lowbval. (Default: --low_bval 5)')  # noqa

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
        Run the DWIPreprocessingUsingPhaseDiffFieldmap pipeline from command line.
        """
        from tempfile import mkdtemp
        from dwi_preprocessing_using_phasediff_fieldmap_pipeline import DWIPreprocessingUsingPhaseDiffFieldmap

        pipeline = DWIPreprocessingUsingPhaseDiffFieldmap(
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
