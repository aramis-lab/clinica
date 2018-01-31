# coding: utf8

import clinica.engine as ce


class DWIProcessingDTICLI(ce.CmdParser):

    def __init__(self):
        super(DWIProcessingDTICLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-processing-dti'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'DTI-based processing of DWI datasets: http://clinica.run/doc/DWIProcessing'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 metavar=('participants.tsv'),
                                 help='TSV file containing the subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 metavar=('Working_Directory'),
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'),
                                 type=int,
                                 help='Number of cores used to run in parallel')
        clinica_opt.add_argument("-sl", "--slurm",
                                 action='store_true',
                                 help='Run the pipelines using SLURM')

    def run_pipeline(self, args):
        """
        """
        from tempfile import mkdtemp
        from dwi_processing_dti_pipeline import DWIProcessingDTI

        pipeline = DWIProcessingDTI(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
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
