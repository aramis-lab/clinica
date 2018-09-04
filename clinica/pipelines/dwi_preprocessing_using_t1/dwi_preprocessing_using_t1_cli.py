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
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("phase_encoding_direction", type=str,
                                  help='The phase encoding direction (e.g. For ADNI data, the phase_encoding_direction is y(j).')
        clinica_comp.add_argument("total_readout_time", type=str,
                                  help='The total readout time (see https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq for details)')

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

        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'epi_param': dict([('readout_time', args.total_readout_time),  ('enc_dir', args.phase_encoding_direction)]),
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
