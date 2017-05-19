"""T1 FreeSurfer2 - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

class T1FreeSurfer2CLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 't1-freesurfer2'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        # Add your own pipeline command line arguments here to be used in the
        # method below. Example below:
        self._args.add_argument("-ras", "--recon_all_args",
                                help='additional flags for recon-all command line, default will be -qcache')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        # END OF EXAMPLE

    def run_pipeline(self, args):
        """
        """
        from t1_freesurfer2_pipeline import T1FreeSurfer2

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = T1FreeSurfer2(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        # pipeline = T1FreeSurfer2()
        pipeline.parameters = {
            # Add your own pipeline parameters here to use them inside your
            # pipeline. See the file `t1_freesurfer2_pipeline.py` to
            # see an example of use.
            'recon_all_args'        : args.recon_all_args or '-qcache'
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()