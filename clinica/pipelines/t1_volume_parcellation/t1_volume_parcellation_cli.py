# coding: utf8

import clinica.engine as ce


class T1VolumeParcellationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 't1-volume-parcellation'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Computation of mean GM concentration for a set of regions:\n' \
                            + 'http://clinica.run/doc/Pipelines/T1_Volume/'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 type=int,
                                 help='Number of cores used to run in parallel')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-m", "--modulation",
                              metavar="[on|off]", default='on',
                              help='Specify if modulation must be enabled (default: --modulation on')
        list_atlases = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
        advanced.add_argument("-atlases", "--atlases",
                              nargs='+', type=str, metavar='',
                              default=list_atlases, choices=list_atlases,
                              help='A list of atlases used to calculate the regional mean GM concentrations (default: all atlases i.e. --atlases AAL2 AICHA Hammers LPBA40 Neuromorphometrics).')

    def run_command(self, args):
        """
        """
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import T1VolumeParcellation

        pipeline = T1VolumeParcellation(
             caps_directory=self.absolute_path(args.caps_directory),
             tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        assert args.modulation in ['on', 'off']
        pipeline.parameters = {
            'group_id': args.group_id,
            'atlases': args.atlases,
            'wd': self.absolute_path(args.working_directory),
            'n_procs': args.n_procs,
            'modulate': args.modulation
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
