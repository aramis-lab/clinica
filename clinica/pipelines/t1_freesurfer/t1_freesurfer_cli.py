# coding: utf8

import clinica.engine as ce


class T1FreeSurferCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "t1-freesurfer"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Cross-sectional pre-processing of T1w images with FreeSurfer:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("bids_directory", help="Path to the BIDS directory.")
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")

        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-raa",
            "--recon_all_args",
            metavar="flag(s)",
            type=str,
            default="-qcache",
            help="Additional flags for recon-all command line "
            '(default: --recon_all_args="%(default)s"). '
            "Please note that = is compulsory after --recon_all_args/-raa flag "
            "(this is not the case for other flags).",
        )

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_overwrite_flag=True)

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .t1_freesurfer_pipeline import T1FreeSurfer

        parameters = {
            "recon_all_args": args.recon_all_args,
        }
        pipeline = T1FreeSurfer(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name,
            overwrite_caps=args.overwrite_outputs,
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
