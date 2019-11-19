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
        self._description = ('Connectome-based processing of DWI datasets:\n'
                             'http://clinica.run/doc/DWI_Connectome')

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
                              help=('Set the desired number of streamlines to generate the tractography and connectome '
                                    '(default: --n_tracks 1000000).'))
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """
        """
        import os
        import datetime
        from colorama import Fore
        from clinica.utils.stream import cprint
        from .dwi_connectome_pipeline import DwiConnectome

        pipeline = DwiConnectome(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory)
        )
        pipeline.parameters = {
            'n_tracks': args.n_tracks or 1000000,
        }
        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        now = datetime.datetime.now().strftime('%H:%M:%S')
        cprint('%s[%s]%s The %s pipeline has completed. You can now delete the working directory (%s).' %
               (Fore.GREEN, now, Fore.RESET, self._name,
                os.path.join(os.path.abspath(args.working_directory), pipeline.__class__.__name__)))
