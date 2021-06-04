# coding: utf8

import clinica.engine as ce


class SpatialSVMCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "machinelearning-prepare-spatial-svm"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Prepare input data for SVM with spatial and anatomical regularization:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/MachineLearning_PrepareSVM/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")
        clinica_comp.add_argument(
            "group_label",
            help="User-defined identifier for the provided group of subjects.",
        )
        clinica_comp.add_argument(
            "orig_input_data",
            help="""Origin of input data. Type
            't1-volume' to use gray matter maps or
            'pet-volume' to use SUVr maps.""",
            choices=["t1-volume", "pet-volume"],
        )
        # Optional arguments for inputs from pet-volume pipeline
        optional_pet = self._args.add_argument_group(
            "Pipeline options if you use inputs from pet-volume pipeline"
        )
        optional_pet.add_argument(
            "-acq",
            "--acq_label",
            type=str,
            default=None,
            help="Name of the label given to the PET acquisition, specifying the tracer used (acq-<acq_label>).",
        )
        optional_pet.add_argument(
            "-suvr",
            "--suvr_reference_region",
            choices=LIST_SUVR_REFERENCE_REGIONS,
            default=None,
            help="Intensity normalization using the average PET uptake in reference regions "
            "resulting in a standardized uptake value ratio (SUVR) map. It can be "
            "cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers).",
        )
        optional_pet.add_argument(
            "-pvc",
            "--use_pvc_data",
            action="store_true",
            default=False,
            help="Use PET data with partial value correction (by default, PET data with no PVC are used)",
        )
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES["ADVANCED"])
        advanced.add_argument(
            "-fwhm",
            "--full_width_half_maximum",
            type=float,
            metavar="N",
            default=4,
            help="Amount of regularization (in mm). In practice, we found the default value "
            "(--full_width_half_maximum %(default)s) to be optimal. We therefore "
            "do not recommend to change it unless you have a specific reason to do so.",
        )

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .spatial_svm_pipeline import SpatialSVM

        if args.orig_input_data == "pet-volume":
            if args.acq_label is None:
                raise ClinicaException(
                    "You selected pet-volume pipeline without setting --acq_label flag. "
                    "Clinica will now exit."
                )
            if args.suvr_reference_region is None:
                raise ClinicaException(
                    "You selected pet-volume pipeline without setting --suvr_reference_region flag. "
                    "Clinica will now exit."
                )

        parameters = {
            # Clinica compulsory arguments
            "group_label": args.group_label,
            "orig_input_data": args.orig_input_data,
            # Optional arguments for inputs from pet-volume pipeline
            "acq_label": args.acq_label,
            "use_pvc_data": args.use_pvc_data,
            "suvr_reference_region": args.suvr_reference_region,
            # Advanced arguments
            "fwhm": args.full_width_half_maximum,
        }
        pipeline = SpatialSVM(
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
