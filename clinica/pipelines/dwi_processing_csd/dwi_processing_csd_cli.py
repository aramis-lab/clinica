# coding: utf8

"""DWIProcessingCSD - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""

import clinica.engine as ce


class DWIProcessingCSDCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'dwi-processing-csd'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Run tractography algorithm using preprocessed ' \
                            'DWI ' \
                            'images:\nhttp://clinica.run/doc/Pipelines' \
                            '/DWIProcessingCSD/'

    def define_options(self):
        """Define the sub-command arguments
        """

        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        clinica_comp = self._args.add_argument_group(
                PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        clinica_opt = self._args.add_argument_group(
                PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipelines intermediate results.')
        clinica_opt.add_argument("-np", "--n_procs", type=int,
                                help='Number of processors to run in parallel.')
        clinica_opt.add_argument("-sl", "--slurm", action='store_true',
                                help='Run the pipelines using SLURM.')
        clinica_opt.add_argument("-sa", "--sbatch_args",
                                help='SLURM\'s sbatch tool arguments.')

        optional = self._args.add_argument_group(
                PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-nt", "--n_tracks", type=int,
                              help='Number of tracks generated (default: 100000).')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from .dwi_processing_csd_pipeline import DWIProcessingCSD

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = DWIProcessingCSD(
                bids_directory=self.absolute_path(args.bids_directory),
                caps_directory=self.absolute_path(args.caps_directory),
                tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'n_tracks'        : args.n_tracks or 100000,
        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        pipeline.write_graph()
        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURMGraph', plugin_args = {
                'dont_resubmit_completed_jobs': True, 'sbatch_args':
                    args.sbatch_args})
        else:
            pipeline.run()
