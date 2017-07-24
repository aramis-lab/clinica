"""T1 SPM Dartel Existing Template - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce


class T1SPMDartelExistingTemplateCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 't1-spm-dartel-existing-template'

    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("group_id",
                                help='Current group name')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-t", "--tissues", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help='Tissues to create flow fields to DARTEL template. Ex: 1 is only GM')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        self._args.add_argument("-sl", "--slurm", action='store_true',
                                help='Run the pipeline using SLURM')

    def run_pipeline(self, args):
        """
        """

        from t1_spm_dartel_existing_template_pipeline import T1SPMDartelExistingTemplate

        pipeline = T1SPMDartelExistingTemplate(bids_directory=self.absolute_path(args.bids_directory),
                                               caps_directory=self.absolute_path(args.caps_directory),
                                               tsv_file=self.absolute_path(args.subjects_sessions_tsv),
                                               group_id=args.group_id
                                               )

        pipeline.parameters.update({'tissues': args.tissues})

        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
