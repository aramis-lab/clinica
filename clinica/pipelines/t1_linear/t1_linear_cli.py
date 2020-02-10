# coding: utf8

"""T1 Linear - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


import clinica.engine as ce


class T1LinearCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 't1-linear'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = ('Brief description:\n'
                             'http://clinica.run/doc/Pipelines/T1Linear/')

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id...)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        # group_id can be used by certain pipelines when some operations are performed at the group level
        # (for example, generation of a template in pipeline t1-volume)
        # clinica_comp.add_argument("group_id",
        #                          help='User-defined identifier for the provided group of subjects.')

        # Clinica standard arguments
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results.')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel.')

        # Add your own pipeline command line arguments here to be used in the
        # method below. Example below:
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-rt", "--ref_template",
                              help='Reference template for registration.')

        # Add advanced arguments
        # advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        # advanced.add_argument("-aa", "--advanced_arg",
        #                       help='Your advanced argument.')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from .t1_linear_pipeline import T1Linear

        parameters = {
            'ref_template'        : args.ref_template or 'Reference Template'
        }

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = T1Linear(
             bids_directory=self.absolute_path(args.bids_directory),
             caps_directory=self.absolute_path(args.caps_directory),
             tsv_file=self.absolute_path(args.subjects_sessions_tsv),
             parameters=parameters,
             name=self.name,
             overwrite_caps=args.overwrite_outputs)
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
