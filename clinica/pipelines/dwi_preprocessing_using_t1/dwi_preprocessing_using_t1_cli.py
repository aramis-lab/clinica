# coding: utf8

import clinica.engine as ce


class DwiPreprocessingUsingT1Cli(ce.CmdParser):

    def __init__(self):
        super(DwiPreprocessingUsingT1Cli, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'dwi-preprocessing-using-t1'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Preprocessing of raw DWI datasets using a T1w image:\n'
                             'http://clinica.run/doc/Pipelines/DWI_Preprocessing/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("--low_bval",
                              metavar=('N'), type=int, default=5,
                              help='Define the b0 volumes as all volume bval <= low_bval '
                                   '(default: --low_bval %(default)s).')

        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results.')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel.')

        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        cuda_action = advanced.add_mutually_exclusive_group(required=False)
        cuda_action.add_argument('--use_cuda_8_0', action='store_true', default=False,
                                 help=('Use CUDA 8.0 implementation of FSL eddy. Please note that '
                                       '--use_cuda_8_0 and --use_cuda_9_1 flags are mutually exclusive.'))
        cuda_action.add_argument('--use_cuda_9_1', action='store_true', default=False,
                                 help=('Use CUDA 9.1 implementation of FSL eddy. Please note that '
                                       '--use_cuda_8_0 and --use_cuda_9_1 flags are mutually exclusive.'))
        advanced.add_argument("--initrand",
                              metavar=('N'), type=int,
                              help="Initialize the random number generator for FSL eddy.")

    def run_command(self, args):
        """Run the  Pipeline from command line."""
        import os
        import datetime
        from tempfile import mkdtemp
        from colorama import Fore
        from clinica.utils.stream import cprint
        from .dwi_preprocessing_using_t1_pipeline import DwiPreprocessingUsingT1

        pipeline = DwiPreprocessingUsingT1(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            low_bval=args.low_bval,
            use_cuda_8_0=args.use_cuda_8_0,
            use_cuda_9_1=args.use_cuda_9_1,
            seed_fsl_eddy=args.initrand,
        )

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        now = datetime.datetime.now().strftime('%H:%M:%S')
        cprint('%s[%s]%s The %s pipeline has completed. You can now delete the working directory (%s).' %
               (Fore.GREEN, now, Fore.RESET, self._name,
                os.path.join(os.path.abspath(args.working_directory), pipeline.__class__.__name__)))
