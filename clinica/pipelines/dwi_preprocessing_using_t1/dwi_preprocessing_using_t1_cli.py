# coding: utf8

import clinica.engine as ce


class DwiPreprocessingUsingT1Cli(ce.CmdParser):

    def __init__(self):
        super(DwiPreprocessingUsingT1Cli, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-preprocessing-using-t1'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = ('Preprocessing of raw DWI datasets using a T1w image:\n'
                             'http://clinica.run/doc/Pipelines/DWI_Preprocessing/')

    def define_options(self):
        """Define the sub-command arguments.
        """
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
                              help='Define the b0 volumes as all volume bval <= lowbval. (default: --low_bval 5)')  # noqa
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """
        Run the DWIPreprocessingUsingT1 Pipeline from command line.
        """
        from clinica.utils.stream import cprint
        from .dwi_preprocessing_using_t1_pipeline import DwiPreprocessingUsingT1

        pipeline = DwiPreprocessingUsingT1(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            low_bval=args.low_bval,
        )

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
