# coding: utf8

import clinica.engine as ce


class T1LinearCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""

        self._name = "t1-linear"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Affine registration of T1w images to the MNI standard space:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Linear/"
        )

    def define_options(self):
        """Define the sub-command arguments"""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label...)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("bids_directory", help="Path to the BIDS directory.")
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")

        # Clinica optional arguments
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])

        optional.add_argument(
            "-ui",
            "--uncropped_image",
            help="""Do not crop the image with template
                              (cropped image are suggested for using with DL
                              models)""",
            action="store_true",
            default=False,
        )

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .t1_linear_pipeline import T1Linear

        parameters = {"uncropped_image": args.uncropped_image}

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = T1Linear(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name,
        )

        if args.n_procs:
            exec_pipeline = pipeline.run(
                plugin="MultiProc", plugin_args={"n_procs": args.n_procs}
            )
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(
                self.name, pipeline.base_dir, pipeline.base_dir_was_specified
            )
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
