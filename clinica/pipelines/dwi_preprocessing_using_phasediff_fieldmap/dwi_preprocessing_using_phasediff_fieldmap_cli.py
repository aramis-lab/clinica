# coding: utf8

import clinica.engine as ce


class DwiPreprocessingUsingPhaseDiffFieldmapCli(ce.CmdParser):

    def __init__(self):
        super(DwiPreprocessingUsingPhaseDiffFieldmapCli, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-preprocessing-using-fieldmap'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Preprocessing of raw DWI datasets using a phase difference image:\nhttp://clinica.run/doc/Pipelines/DWI_Preprocessing/'

    def define_options(self):
        """Define the sub-command arguments
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
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel')

    def run_command(self, args):
        """
        Run the DWIPreprocessingUsingPhaseDiffFieldmap pipeline from command line.
        """
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from .dwi_preprocessing_using_phasediff_fieldmap_pipeline import DwiPreprocessingUsingPhaseDiffFieldmap

        pipeline = DwiPreprocessingUsingPhaseDiffFieldmap(
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
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
