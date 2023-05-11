import clinica.pipelines.engine as cpe


class SpatialSVM(cpe.Pipeline):
    """SpatialSVM - Prepare input data for SVM with spatial and anatomical regularization.

    Returns:
        A clinica pipeline object containing the SpatialSVM pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.group import check_group_label

        # Clinica compulsory parameters
        self.parameters.setdefault("group_label", None)
        check_group_label(self.parameters["group_label"])

        if "orig_input_data_ml" not in self.parameters.keys():
            raise KeyError(
                "Missing compulsory orig_input_data key in pipeline parameter."
            )

        # Optional parameters for inputs from pet-volume pipeline
        self.parameters.setdefault("acq_label", None)
        self.parameters.setdefault("suvr_reference_region", None)
        self.parameters.setdefault("use_pvc_data", False)

        # Advanced parameters
        self.parameters.setdefault("fwhm", 4)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["dartel_input", "input_image"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["regularized_image"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from clinica.utils.input_files import (
            pet_volume_normalized_suvr_pet,
            t1_volume_final_group_template,
        )
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        from clinica.utils.ux import print_groups_in_caps_directory

        # Check that group already exists
        if not os.path.exists(
            os.path.join(
                self.caps_directory, "groups", "group-" + self.parameters["group_label"]
            )
        ):
            print_groups_in_caps_directory(self.caps_directory)
            raise ClinicaException(
                f"Group {self.parameters['group_label']} does not exist. "
                "Did you run pet-volume, t1-volume or t1-volume-create-dartel pipeline?"
            )

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
        )
        all_errors = []

        if self.parameters["orig_input_data_ml"] == "t1-volume":
            caps_files_information = {
                "pattern": os.path.join(
                    "t1",
                    "spm",
                    "dartel",
                    "group-" + self.parameters["group_label"],
                    "*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability.nii.gz",
                ),
                "description": "graymatter tissue segmented in T1w MRI in Ixi549 space",
                "needed_pipeline": "t1-volume-tissue-segmentation",
            }
        elif self.parameters["orig_input_data_ml"] == "pet-volume":
            if not (
                self.parameters["acq_label"]
                and self.parameters["suvr_reference_region"]
            ):
                raise ValueError(
                    f"Missing value(s) in parameters from pet-volume pipeline. Given values:\n"
                    f"- acq_label: {self.parameters['acq_label']}\n"
                    f"- suvr_reference_region: {self.parameters['suvr_reference_region']}\n"
                    f"- use_pvc_data: {self.parameters['use_pvc_data']}\n"
                )
            caps_files_information = pet_volume_normalized_suvr_pet(
                acq_label=self.parameters["acq_label"],
                suvr_reference_region=self.parameters["suvr_reference_region"],
                use_brainmasked_image=False,
                use_pvc_data=self.parameters["use_pvc_data"],
                fwhm=0,
            )
        else:
            raise ValueError(
                f"Image type {self.parameters['orig_input_data_ml']} unknown."
            )

        try:
            input_image, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                caps_files_information,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            dartel_input = clinica_group_reader(
                self.caps_directory,
                t1_volume_final_group_template(self.parameters["group_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        # Raise all errors if some happened
        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your CAPS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaCAPSError(error_message)

        read_parameters_node.inputs.dartel_input = dartel_input
        read_parameters_node.inputs.input_image = input_image

        # fmt: off
        self.connect(
            [
                (read_parameters_node, self.input_node, [("dartel_input", "dartel_input")]),
                (read_parameters_node, self.input_node, [("input_image", "input_image")]),
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils

        fisher_tensor_generation = npe.Node(
            name="obtain_g_fisher_tensor",
            interface=nutil.Function(
                input_names=["dartel_input", "FWHM"],
                output_names=["fisher_tensor", "fisher_tensor_path"],
                function=utils.obtain_g_fisher_tensor,
            ),
        )
        fisher_tensor_generation.inputs.FWHM = self.parameters["fwhm"]

        time_step_generation = npe.Node(
            name="estimation_time_step",
            interface=nutil.Function(
                input_names=["dartel_input", "FWHM", "g"],
                output_names=["t_step", "json_file"],
                function=utils.obtain_time_step_estimation,
            ),
        )
        time_step_generation.inputs.FWHM = self.parameters["fwhm"]

        heat_solver_equation = npe.MapNode(
            name="heat_solver_equation",
            interface=nutil.Function(
                input_names=["input_image", "g", "FWHM", "t_step", "dartel_input"],
                output_names=["regularized_image"],
                function=utils.heat_solver_equation,
            ),
            iterfield=["input_image"],
        )
        heat_solver_equation.inputs.FWHM = self.parameters["fwhm"]

        datasink = npe.Node(nio.DataSink(), name="sinker")
        datasink.inputs.base_directory = self.caps_directory
        datasink.inputs.parameterization = True
        if self.parameters["orig_input_data_ml"] == "t1-volume":
            datasink.inputs.regexp_substitutions = [
                (
                    r"(.*)/regularized_image/.*/(.*(sub-(.*)_ses-(.*))_T1w(.*)_probability(.*))$",
                    r"\1/subjects/sub-\4/ses-\5/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"/\3_T1w\6_spatialregularization\7",
                ),
                (
                    r"(.*)json_file/(output_data.json)$",
                    r"\1/groups/group-"
                    + self.parameters["group_label"]
                    + r"/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"_space-Ixi549Space_parameters.json",
                ),
                (
                    r"(.*)fisher_tensor_path/(output_fisher_tensor.npy)$",
                    r"\1/groups/group-"
                    + self.parameters["group_label"]
                    + r"/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"_space-Ixi549Space_gram.npy",
                ),
            ]

        elif self.parameters["orig_input_data_ml"] == "pet-volume":
            datasink.inputs.regexp_substitutions = [
                (
                    r"(.*)/regularized_image/.*/(.*(sub-(.*)_ses-(.*))_(task.*)_pet(.*))$",
                    r"\1/subjects/sub-\4/ses-\5/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"/\3_\6_spatialregularization\7",
                ),
                (
                    r"(.*)json_file/(output_data.json)$",
                    r"\1/groups/group-"
                    + self.parameters["group_label"]
                    + r"/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"_space-Ixi549Space_parameters.json",
                ),
                (
                    r"(.*)fisher_tensor_path/(output_fisher_tensor.npy)$",
                    r"\1/groups/group-"
                    + self.parameters["group_label"]
                    + r"/machine_learning/input_spatial_svm/group-"
                    + self.parameters["group_label"]
                    + r"_space-Ixi549Space_gram.npy",
                ),
            ]
        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, fisher_tensor_generation, [("dartel_input", "dartel_input")]),
                (fisher_tensor_generation, time_step_generation, [("fisher_tensor", "g")]),
                (self.input_node, time_step_generation, [("dartel_input", "dartel_input")]),
                (self.input_node, heat_solver_equation, [("input_image", "input_image")]),
                (fisher_tensor_generation, heat_solver_equation, [("fisher_tensor", "g")]),
                (time_step_generation, heat_solver_equation, [("t_step", "t_step")]),
                (self.input_node, heat_solver_equation, [("dartel_input", "dartel_input")]),
                (fisher_tensor_generation, datasink, [("fisher_tensor_path", "fisher_tensor_path")]),
                (time_step_generation, datasink, [("json_file", "json_file")]),
                (heat_solver_equation, datasink, [("regularized_image", "regularized_image")]),
            ]
        )
        # fmt: on
