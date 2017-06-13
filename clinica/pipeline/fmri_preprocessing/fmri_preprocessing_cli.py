"""fMRI Preprocessing - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

class fMRIPreprocessingCLI(ce.CmdParser):
    """
    """

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'fmri-preprocessing'


    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of processors to run in parallel')
        self._args.add_argument("-ns", "--num_slices", type=int,
                                help="Number of slices")
        self._args.add_argument("-tr", "--time_repetition", type=float,
                                help='TR in seconds')
        self._args.add_argument("-et", "--echo_times", nargs=2, type=float,
                                help="Echo times in seconds (ex.: '-et 5.19 7.65')")
        self._args.add_argument("-bd", "--blip_direction", type=int,
                                help="Blip direction (1 or -1)")
        self._args.add_argument("-trt", "--total_readout_time", type=float,
                                help="Total readout time (TRT) in seconds")


    def run_pipeline(self, args):
        """
        """

        from fmri_preprocessing_pipeline import fMRIPreprocessing

        pipeline = fMRIPreprocessing(bids_directory=self.absolute_path(args.bids_directory),
                                     caps_directory=self.absolute_path(args.caps_directory),
                                     tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'num_slices'        : args.num_slices,
            'time_repetition'   : args.time_repetition,
            'echo_times'        : args.echo_times,
            'blip_direction'    : args.blip_direction,
            'total_readout_time': args.total_readout_time
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()