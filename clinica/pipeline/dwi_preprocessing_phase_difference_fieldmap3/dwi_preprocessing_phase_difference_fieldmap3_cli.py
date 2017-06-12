"""Dwi Preprocessing Phase Difference Fieldmap3 - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce


class DwiPreprocessingPhaseDifferenceFieldmap3CLI(ce.CmdParser):


    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'dwi-preprocessing-phase-difference-fieldmap3'


    def define_options(self):
        """Define the sub-command arguments
        """


        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        # registration = self._args.add_mutually_exclusive_group(required=True)
        # registration.add_argument('-register_fmap_on_b0', action='store_true',
        #                    help='Choose to register fmap on b0. (mutually exclusive with \'-do_not_register_fmap_on_b0\').')
        # registration.add_argument('-do_not_register_fmap_on_b0', action='store_false',
        #                    help='Choose not to register fmap on b0. (mutually exclusive with \'-register_fmap_on_b0\').')
        self._args.add_argument("-ees",  "--echospacing", type=float, default=0.36e-3,
                                help='Effective echo spacing time of epi sequence')
        self._args.add_argument("-det",  "--delta_te", type=float, default=2.46e-3,
                                help='Echo time difference between two phase images')
        self._args.add_argument("-ped", "--enc_dir", type=str, default='y-',
                                help='Phase encoding direction of epi sequence')
        ## TODO, exclued with other epi options, add the option to read from json file
        # self._args.add_argument("-epijson", "--epi_read_from_json", type=str, default=False,
        #                         help='Reading epi paramters from json file if True')
        self._args.add_argument("-rfmb0", "--register_fmap_on_b0", type=str, default=True,
                                help='Register fmap on b0, default is Ture')
        self._args.add_argument("-rb0", "--register_b0_on_b0", type=str, default=True,
                                help='Register all b0s on first b0, default is Ture')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of cores used to run in parallel')


    def run_pipeline(self, args):
        """
        """

        from tempfile import mkdtemp
        from dwi_preprocessing_phase_difference_fieldmap3_pipeline import DwiPreprocessingPhaseDifferenceFieldmap3
        from nipype import config
        cfg = dict(execution={'parameterize_dirs': False})
        config.update_config(cfg)

        pipeline = DwiPreprocessingPhaseDifferenceFieldmap3(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'register_fmap_on_b0': args.register_fmap_on_b0,
            'register_b0_on_b0': args.register_b0_on_b0,
            'echospacing': args.echospacing,
            'delta_te': args.delta_te,
            'enc_dir': args.enc_dir,
            'n_procs': args.n_procs
        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()