# coding: utf8

import clinica.engine as ce


class DwiConnectomeCli(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-connectome'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Connectome-based processing of DWI datasets:\nhttp://clinica.run/doc/DWI_Connectome'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-nt", "--n_tracks",
                              metavar=('N'), type=int,
                              help='Set the desired number of streamlines to generate the tractography and connectome (default: 1M --n_tracks 1000000).')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results.')
        clinica_opt.add_argument("-np", "--n_procs",
                                 type=int,
                                 help='Number of cores used to run in parallel.')
        clinica_opt.add_argument("-overwrite", "--overwrite_outputs",
                                 action='store_true', default=False,
                                 help='Force overwrite of output files (can not be used together with --skip_if_outputs_present flag).')
        clinica_opt.add_argument("-skip", "--skip_if_outputs_present",
                                 action='store_true', default=False,
                                 help='Skip input image if its outputs are present (can not be used together with --overwrite_outputs flag).')

    def run_command(self, args):
        """
        """
        import os
        import datetime
        import sys
        from colorama import Fore
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from .dwi_connectome_pipeline import DwiConnectome

        if args.overwrite_outputs and args.skip_if_outputs_present:
            cprint("%s\n[Error] You can not use the --skip_if_outputs_present flag and --overwrite_outputs flag at the same time.%s\n"
                   "\n%sExplanations on the flags:%s\n"
                   "\t%s--skip_if_outputs_present%s flag will skip input image if its outputs are present.\n"
                   "\t%s--overwrite_outputs%s flag will force overwrite of output files.\n"
                   "\nEither remove these two flags or chose a single one. The program will now exit." %
                   (Fore.RED, Fore.RESET, Fore.BLUE, Fore.RESET, Fore.BLUE, Fore.RESET, Fore.BLUE, Fore.RESET))
            sys.exit(1)

        pipeline = DwiConnectome(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
        )
        pipeline.parameters = {
            'n_tracks'                  : args.n_tracks or 1000000,
            'overwrite_outputs'         : args.overwrite_outputs,
            'skip_if_outputs_present'   : args.skip_if_outputs_present,
        }
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
