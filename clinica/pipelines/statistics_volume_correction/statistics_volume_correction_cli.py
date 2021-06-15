# coding: utf8

import clinica.engine as ce


class StatisticsVolumeCorrectionCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "statistics-volume-correction"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Statistical correction of statistics-volume pipeline:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Volume/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")
        clinica_comp.add_argument("t_map", help="t-statistics map")
        clinica_comp.add_argument(
            "height_threshold",
            type=float,
            help="T value corresponding to an uncorrected p-value of 0.001",
        )
        clinica_comp.add_argument(
            "FWEp",
            type=float,
            help="height threshold (i.e. voxel-level (= peak) threshold)",
        )
        clinica_comp.add_argument(
            "FDRp",
            type=float,
            help="height threshold (i.e. voxel-level (= peak) threshold)",
        )
        clinica_comp.add_argument(
            "FWEc", type=int, help="extent threshold (i.e. cluster size threshold)"
        )
        clinica_comp.add_argument(
            "FDRc", type=int, help="extent threshold (i.e. cluster size threshold)"
        )

        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-nc",
            "--n_cuts",
            default=8,
            type=int,
            help="Number of cuts along each direction",
        )

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .statistics_volume_correction_pipeline import StatisticsVolumeCorrection

        parameters = {
            "t_map": args.t_map,
            "height_threshold": args.height_threshold,
            "FWEp": args.FWEp,
            "FDRp": args.FDRp,
            "FWEc": args.FWEc,
            "FDRc": args.FDRc,
            "n_cuts": args.n_cuts,
        }

        pipeline = StatisticsVolumeCorrection(
            caps_directory=self.absolute_path(args.caps_directory),
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
