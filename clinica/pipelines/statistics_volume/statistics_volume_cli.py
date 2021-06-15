# coding: utf-8

import clinica.engine as ce


class StatisticsVolumeCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "statistics-volume"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Volume-based mass-univariate analysis with SPM:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Volume/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from colorama import Fore

        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS

        # Clinica compulsory arguments
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
            help="Type of volume-based feature: type "
            "'t1-volume' to use gray matter maps, "
            "'pet-volume' to use PET data or "
            "'custom-pipeline' to use you own data in CAPS directory "
            "(see Wiki for details).",
            choices=["t1-volume", "pet-volume", "custom-pipeline"],
        )
        clinica_comp.add_argument(
            "subject_visits_with_covariates_tsv",
            help="TSV file containing a list of subjects with their sessions and all "
            "the covariates and factors needed for the GLM.",
        )
        clinica_comp.add_argument(
            "contrast",
            help="Defines the contrast. Must be one of the column names form the TSV file.",
        )

        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-dartel",
            "--group_label_dartel",
            type=str,
            default="*",
            help="Name of the DARTEL template that Clinica needs to use to grab the input files.",
        )
        optional.add_argument(
            "-fwhm",
            "--full_width_at_half_maximum",
            type=int,
            default=8,
            help="Full Width at Half Maximum (FWHM) of the smoothing used in your input files "
            "(default: --full_width_at_half_maximum %(default)s).",
        )

        # Optional arguments for inputs from pet-volume pipeline
        optional_pet = self._args.add_argument_group(
            f"{Fore.BLUE}Pipeline options if you use inputs from pet-volume pipeline{Fore.RESET}"
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

        # Optional arguments for custom pipeline
        opt_custom_input = self._args.add_argument_group(
            f"{Fore.BLUE}Pipeline options if you selected custom-pipeline{Fore.RESET}"
        )
        opt_custom_input.add_argument(
            "-cf",
            "--custom_file",
            type=str,
            default=None,
            help="Custom file string. Specify filename using * when the subject or session name "
            "appears e.g. '*_task-rest_acq-fdg_pet_space-Ixi549Space_pet.nii.gz' will grab "
            "the corresponding file in all the subjects/sessions. "
            "This flag must be specified with the --measure_label flag). "
            "See Wiki for an example.",
        )
        opt_custom_input.add_argument(
            "-ml",
            "--measure_label",
            type=str,
            default=None,
            help="Name of the feature type, it will be saved on the CAPS "
            "_measure-MEASURE_LABEL key-value association. "
            "This flag must be specified with the --custom_file flag). "
            "See Wiki for an example.",
        )

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_tsv_flag=False)

        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES["ADVANCED"])
        advanced.add_argument(
            "-ct",
            "--cluster_threshold",
            type=float,
            default=0.001,
            help="Threshold to define a cluster in the process of cluster-wise correction "
            "(default: --cluster_threshold %(default)s).",
        )

    def run_command(self, args):
        from colorama import Fore
        from networkx import Graph

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .statistics_volume_pipeline import StatisticsVolume

        # PET-Volume pipeline
        if args.orig_input_data == "pet-volume":
            if args.acq_label is None:
                raise ClinicaException(
                    f"{Fore.RED}You selected pet-volume pipeline without setting --acq_label flag. "
                    f"Clinica will now exit.{Fore.RESET}"
                )
            if args.suvr_reference_region is None:
                raise ClinicaException(
                    f"{Fore.RED}You selected pet-volume pipeline without setting --suvr_reference_region flag. "
                    f"Clinica will now exit.{Fore.RESET}"
                )

        # Custom pipeline
        if args.orig_input_data == "custom-pipeline":
            if (args.custom_file is None) or (args.measure_label is None):
                raise ClinicaException(
                    "You must set --measure_label and --custom_file flags."
                )

        parameters = {
            # Clinica compulsory arguments
            "group_label": args.group_label,
            "orig_input_data": args.orig_input_data,
            "contrast": args.contrast,
            # Optional arguments
            "group_label_dartel": args.group_label_dartel,
            "full_width_at_half_maximum": args.full_width_at_half_maximum,
            # Optional arguments for inputs from pet-volume pipeline
            "acq_label": args.acq_label,
            "use_pvc_data": args.use_pvc_data,
            "suvr_reference_region": args.suvr_reference_region,
            # Optional arguments for custom pipeline
            "measure_label": args.measure_label,
            "custom_file": args.custom_file,
            # Advanced arguments
            "cluster_threshold": args.cluster_threshold,
        }

        pipeline = StatisticsVolume(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subject_visits_with_covariates_tsv),
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
