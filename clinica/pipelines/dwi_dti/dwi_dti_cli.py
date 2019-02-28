# coding: utf8

import clinica.engine as ce


class DwiDtiCli(ce.CmdParser):

    def __init__(self):
        super(DwiDtiCli, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'dwi-dti'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'DTI-based processing of DWI datasets:\nhttp://clinica.run/doc/DWI_DTI'

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
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 type=int,
                                 help='Number of cores used to run in parallel')

    def run_command(self, args):
        """
        """
        import os
        import datetime
        from colorama import Fore
        from tempfile import mkdtemp
        from .dwi_dti_pipeline import DwiDti
        from clinica.utils.stream import cprint

        pipeline = DwiDti(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv)
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
