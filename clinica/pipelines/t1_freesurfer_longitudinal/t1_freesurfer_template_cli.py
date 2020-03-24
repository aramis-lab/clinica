# coding: utf8

import clinica.engine as ce


class T1FreeSurferTemplateCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 't1-freesurfer-template'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Creation of unbiased template with FreeSurfer:\n'
                             'http://clinica.run/doc/Pipelines/T1_FreeSurfer_Longitudinal/')

    def define_options(self):
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_overwrite_flag=True)

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .t1_freesurfer_template_pipeline import T1FreeSurferTemplate
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        pipeline = T1FreeSurferTemplate(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            name='t1-freesurfer-template',
            overwrite_caps=args.overwrite_outputs
        )
        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
