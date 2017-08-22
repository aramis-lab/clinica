"""DWI Preprocessing using T1 - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

__author__ = ["Thomas Jacquemont", "Alexandre Routier"]
__copyright__ = "Copyright 2016,2017 The Aramis Lab Team"
__credits__ = ["Nipype", "Thomas Jacquemont", "Alexandre Routier"]
__license__ = "See LICENSE.txt file"
__version__ = "1.0.0"
__status__ = "Development"


class DWIPreprocessingUsingT1CLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'dwi-preprocessing-using-t1'

    def define_options(self):
        """Define the sub-command arguments
        """

        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')

        self._args.add_argument("--low_bval", type=int, default=5,
                                help='Define the b0 volumes as all volume bval <= lowbval. (Default=5)')
        self._args.add_argument("--save_intermediate_files", type=bool, default=False,
                                help='Will save some intermediate results (for debugging purposes).')

        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        self._args.add_argument("-sl", "--slurm", action='store_true',
                                help='Run the pipeline using SLURM')

    def run_pipeline(self, args):
        """
        Run the DWIPreprocessingUsingT1 Pipeline from command line.
        """

        from tempfile import mkdtemp
        from dwi_preprocessing_using_t1_pipeline import DWIPreprocessingUsingT1

        pipeline = DWIPreprocessingUsingT1(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            low_bval=args.low_bval,
            save_intermediate_files=args.save_intermediate_files
        )

#        pipeline.parameters.update({
#            'low_bval': args.low_bval,
#            'save_intermediate_files': args.save_intermediate_files
#        })

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
