# coding: utf8

import clinica.engine as ce


class T1LinearCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 't1-linear'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = ('Affine registration of T1w images to the MNI standard space:\n'
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

        # Clinica optional arguments
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])

        optional.add_argument("-cp", "--crop_image",
                              help='Crop the image using a template (suggested for using with DL models)',
                              action='store_true',
                              default=False)

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .t1_linear_pipeline import T1Linear
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
                'crop_image': args.crop_image
        }

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = T1Linear(
             bids_directory=self.absolute_path(args.bids_directory),
             caps_directory=self.absolute_path(args.caps_directory),
             tsv_file=self.absolute_path(args.subjects_sessions_tsv),
             base_dir=self.absolute_path(args.working_directory),
             parameters=parameters,
             name=self.name
             )

        if args.n_procs:
            exec_pipeline = pipeline.run(
                    plugin='MultiProc',
                    plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
