import clinica.pipelines.engine as cpe


class StatisticsSurface(cpe.Pipeline):
    """StatisticsSurface - Surface-based mass-univariate analysis with SurfStat.

    See documentation at https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Surface/

    Note:
        The `tsv_file` attribute is overloaded for this pipeline. It must contain a list of subjects
        with their sessions and all the covariates and factors needed for the GLM.

        Pipeline parameters are explained in StatisticsSurfaceCLI.define_options()

    Returns:
        A clinica pipeline object containing the StatisticsSurface pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.group import check_group_label

        from ._inputs import _get_t1_freesurfer_custom_file_template

        # Clinica compulsory parameters
        self.parameters.setdefault("group_label", None)
        check_group_label(self.parameters["group_label"])

        if "orig_input_data" not in self.parameters.keys():
            raise KeyError(
                "Missing compulsory orig_input_data key in pipeline parameter."
            )

        self.parameters.setdefault("glm_type", None)
        if self.parameters["glm_type"] not in ["group_comparison", "correlation"]:
            raise ClinicaException(
                f"The glm_type you specified is wrong: it should be group_comparison or "
                f"correlation (given value: {self.parameters['glm_type']})."
            )

        if "contrast" not in self.parameters.keys():
            raise KeyError("Missing compulsory contrast key in pipeline parameter.")

        # Optional parameters
        self.parameters.setdefault("covariates", None)
        self.parameters.setdefault("full_width_at_half_maximum", 20)

        # Optional parameters for inputs from pet-surface pipeline
        self.parameters.setdefault("acq_label", None)
        self.parameters.setdefault("suvr_reference_region", None)

        # Optional parameters for custom pipeline
        self.parameters.setdefault(
            "custom_file", _get_t1_freesurfer_custom_file_template(self.base_dir)
        )
        self.parameters.setdefault("measure_label", "ct")

        # Advanced parameters
        self.parameters.setdefault("cluster_threshold", 0.001)
        if (
            self.parameters["cluster_threshold"] < 0
            or self.parameters["cluster_threshold"] > 1
        ):
            raise ClinicaException(
                f"Cluster threshold should be between 0 and 1 "
                f"(given value: {self.parameters['cluster_threshold']})."
            )

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return []

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["output_dir"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint

        # Check if already present in CAPS
        # ================================
        # Check if the group label has been existed, if yes, give an error to the users
        # Note(AR): if the user wants to compare Cortical Thickness measure with PET measure
        # using the group_id, Clinica won't allow it.
        # TODO: Modify this behaviour
        if os.path.exists(
            os.path.join(
                self.caps_directory, "groups", f"group-{self.parameters['group_label']}"
            )
        ):
            raise ClinicaException(
                f"Group label {self.parameters['group_label']} already exists, "
                f"please choose another one or delete the existing folder and "
                f"also the working directory and rerun the pipeline"
            )

        # Check input files before calling SurfStat with Matlab
        # =====================================================
        all_errors = []
        # clinica_files_reader expects regexp to start at subjects/ so sub-*/ses-*/ is removed here
        fwhm = str(self.parameters["full_width_at_half_maximum"])
        for direction, hemi in zip(["left", "right"], ["lh", "rh"]):
            cut_pattern = "sub-*/ses-*/"
            query = {"subject": "sub-*", "session": "ses-*", "hemi": hemi, "fwhm": fwhm}
            pattern_hemisphere = self.parameters["custom_file"] % query
            surface_based_info = {
                "pattern": pattern_hemisphere[
                    pattern_hemisphere.find(cut_pattern) + len(cut_pattern) :
                ],
                "description": f"surface-based features on {direction} hemisphere at FWHM = {fwhm}",
            }
            try:
                clinica_file_reader(
                    self.subjects,
                    self.sessions,
                    self.caps_directory,
                    surface_based_info,
                )
            except ClinicaException as e:
                all_errors.append(e)
        # Raise all errors if something happened
        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your CAPS directory.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise RuntimeError(error_message)

        # Give pipeline info
        # ==================
        cprint(
            "The pipeline will last a few minutes. Images generated by Matlab will popup during the pipeline."
        )

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .statistics_surface_utils import save_to_caps

        # Writing results into CAPS
        # =========================
        save_to_caps = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "source_dir",
                    "caps_dir",
                    "overwrite_caps",
                    "pipeline_parameters",
                ],
                function=save_to_caps,
            ),
            name="SaveToCaps",
        )
        save_to_caps.inputs.caps_dir = self.caps_directory
        save_to_caps.inputs.overwrite_caps = self.overwrite_caps
        save_to_caps.inputs.pipeline_parameters = self.parameters

        self.connect(
            [
                (self.output_node, save_to_caps, [("output_dir", "source_dir")]),
            ]
        )

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.pipelines.statistics_surface.statistics_surface_utils as utils

        from .statistics_surface_utils import init_input_node

        init_input = npe.Node(
            interface=nutil.Function(
                input_names=["parameters", "base_dir", "subjects_visits_tsv"],
                output_names=["group_label", "surfstat_results_dir"],
                function=init_input_node,
            ),
            name="0-InitPipeline",
        )
        init_input.inputs.parameters = self.parameters
        init_input.inputs.base_dir = os.path.join(self.base_dir, self.name)
        init_input.inputs.subjects_visits_tsv = self.tsv_file

        # Node to wrap the SurfStat matlab script
        surfstat = npe.Node(
            name="1-RunSurfStat",
            interface=nutil.Function(
                input_names=[
                    "caps_dir",
                    "output_dir",
                    "subjects_visits_tsv",
                    "pipeline_parameters",
                ],
                output_names=["output_dir"],
                function=utils.run_matlab,
            ),
        )
        surfstat.inputs.caps_dir = self.caps_directory
        surfstat.inputs.subjects_visits_tsv = self.tsv_file
        surfstat.inputs.pipeline_parameters = self.parameters

        # Connection
        # ==========
        self.connect(
            [
                (init_input, surfstat, [("surfstat_results_dir", "output_dir")]),
                (surfstat, self.output_node, [("output_dir", "output_dir")]),
            ]
        )
