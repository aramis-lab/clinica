# coding: utf8

import clinica.engine as ce


class T1VolumeDartel2MNICLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "t1-volume-dartel2mni"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Register DARTEL template to MNI space:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/"
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
        clinica_comp.add_argument(
            "group_label",
            help="User-defined identifier for the provided group of subjects.",
        )
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-s",
            "--smooth",
            nargs="+",
            type=int,
            default=[8],
            help="A list of integers specifying the different isomorphic FWHM in millimeters "
            "to smooth the image (default: --smooth 8).",
        )
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES["ADVANCED"])
        advanced.add_argument(
            "-t",
            "--tissues",
            metavar="",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissues to create flow fields to DARTEL template "
            "(default: GM, WM and CSF i.e. --tissues 1 2 3).",
        )
        advanced.add_argument(
            "-m",
            "--modulate",
            type=bool,
            default=True,
            metavar=("True/False"),
            help="A boolean. Modulate output images - no modulation preserves concentrations "
            "(default: --modulate True).",
        )
        advanced.add_argument(
            "-vs",
            "--voxel_size",
            metavar=("float"),
            nargs=3,
            type=float,
            help="A list of 3 floats specifying the voxel size of the output image "
            "(default: --voxel_size 1.5 1.5 1.5).",
        )

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .t1_volume_dartel2mni_pipeline import T1VolumeDartel2MNI

        parameters = {
            "group_label": args.group_label,
            "tissues": args.tissues,
            "voxel_size": args.voxel_size,
            "modulate": args.modulate,
            "smooth": args.smooth,
        }
        pipeline = T1VolumeDartel2MNI(
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
