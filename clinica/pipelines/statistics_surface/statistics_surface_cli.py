# coding: utf8

import clinica.engine as ce


class StatisticsSurfaceCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "statistics-surface"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Surface-based mass-univariate analysis with SurfStat:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Surface/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
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
            help="Type of surface-based feature: type "
            "'t1-freesurfer' to use cortical thickness, "
            "'pet-surface' to use projected PET data or "
            "'custom-pipeline' to use you own data in CAPS directory "
            "(see Wiki for details).",
            choices=["t1-freesurfer", "pet-surface", "custom-pipeline"],
        )
        clinica_comp.add_argument(
            "glm_type",
            help="String based on the GLM type for the hypothesis. "
            "You can choose between 'group_comparison' and 'correlation'.",
            choices=["group_comparison", "correlation"],
        )
        clinica_comp.add_argument(
            "subject_visits_with_covariates_tsv",
            help="TSV file containing a list of subjects with their sessions and all "
            "the covariates and factors needed for the GLM.",
        )
        clinica_comp.add_argument(
            "contrast",
            help="String to define the contrast matrix for the GLM, e.g. group. Please note "
            "that, when you want to perform negative correlation, the sign is ignored "
            "by the command line.",
        )
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-c",
            "--covariates",
            type=str,
            default=None,
            metavar="'covariate_1 ... covariate_n'",
            help="Covariates of the form --covariates 'covariate_1 ... covariate_n'. "
            "Each covariate must match the column name of the TSV file. "
            "On default, no covariates are taken.",
        )
        optional.add_argument(
            "-fwhm",
            "--full_width_at_half_maximum",
            type=int,
            default=20,
            help="FWHM for the surface smoothing "
            "(default: --full_width_at_half_maximum %(default)s).",
        )
        # Optional arguments for inputs from pet-surface pipeline
        opt_pet = self._args.add_argument_group("Pipeline options if you use inputs from pet-surface pipeline")
        opt_pet.add_argument(
            "-acq",
            "--acq_label",
            type=str,
            default=None,
            help="Name of the label given to the PET acquisition, specifying the tracer used (acq-<acq_label>).",
        )
        opt_pet.add_argument(
            "-suvr",
            "--suvr_reference_region",
            choices=LIST_SUVR_REFERENCE_REGIONS,
            default=None,
            help="Intensity normalization using the average PET uptake in reference regions "
            "resulting in a standardized uptake value ratio (SUVR) map. It can be "
            "cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers).",
        )
        # Optional arguments for custom pipeline
        opt_custom_input = self._args.add_argument_group("Pipeline options if you selected custom-pipeline")
        opt_custom_input.add_argument(
            "-cf",
            "--custom_file",
            type=str,
            default=None,
            help="Pattern of file inside CAPS directory using @subject, @session, "
            "@fwhm, @hemi. "
            "This flag must be specified with the --measure_label flag). "
            "See Wiki for an example.",
        )
        opt_custom_input.add_argument(
            "-ml",
            "--measure_label",
            type=str,
            default=None,
            help="Name of the feature type, it will be saved on the CAPS "
            "_measure-FEATURE_LABEL key-value association. "
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
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .statistics_surface_pipeline import StatisticsSurface
        from .statistics_surface_utils import (
            get_pet_surface_custom_file,
            get_t1_freesurfer_custom_file,
        )

        # PET-Surface pipeline
        if args.orig_input_data == "pet-surface":
            if args.acq_label is None:
                raise ClinicaException(
                    "You selected pet-surface pipeline without setting --acq_label flag. "
                    "Clinica will now exit."
                )
            if args.suvr_reference_region is None:
                raise ClinicaException(
                    "You selected pet-surface pipeline without setting --suvr_reference_region flag. "
                    "Clinica will now exit."
                )

        # FreeSurfer cortical thickness
        if args.orig_input_data == "t1-freesurfer":
            args.custom_file = get_t1_freesurfer_custom_file()
            args.measure_label = "ct"
        # PET cortical projection
        elif args.orig_input_data == "pet-surface":
            args.custom_file = get_pet_surface_custom_file(
                args.acq_label, args.suvr_reference_region
            )
            args.measure_label = args.acq_label
        else:
            if (args.custom_file is None) or (args.measure_label is None):
                raise ClinicaException(
                    "You must set --measure_label and --custom_file flags."
                )

        parameters = {
            # Clinica compulsory arguments
            "group_label": args.group_label,
            "orig_input_data": args.orig_input_data,
            "glm_type": args.glm_type,
            "contrast": args.contrast,
            # Optional arguments
            "covariates": args.covariates,
            "full_width_at_half_maximum": args.full_width_at_half_maximum,
            # Optional arguments for inputs from pet-surface pipeline
            "acq_label": args.acq_label,
            "suvr_reference_region": args.suvr_reference_region,
            # Optional arguments for custom pipeline
            "custom_file": args.custom_file,
            "measure_label": args.measure_label,
            # Advanced arguments (i.e. tricky parameters)
            "cluster_threshold": args.cluster_threshold,
        }
        pipeline = StatisticsSurface(
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
