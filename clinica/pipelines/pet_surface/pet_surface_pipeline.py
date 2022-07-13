import clinica.pipelines.engine as cpe


class PetSurface(cpe.Pipeline):
    """PetSurface - Surface-based processing of PET images.

    Returns:
        A clinica pipeline object containing the PetSurface pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        if "acq_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory acq_label key in pipeline parameter.")
        if "suvr_reference_region" not in self.parameters.keys():
            raise KeyError(
                "Missing compulsory suvr_reference_region key in pipeline parameter."
            )
        if "pvc_psf_tsv" not in self.parameters.keys():
            raise KeyError("Missing compulsory pvc_psf_tsv key in pipeline parameter.")

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return [
            "orig_nu",
            "pet",
            "white_surface_left",
            "white_surface_right",
            "destrieux_left",
            "destrieux_right",
            "desikan_left",
            "desikan_right",
        ]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline."""
        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        if self.parameters["longitudinal"]:
            self.build_input_node_longitudinal()
        else:
            self.build_input_node_cross_sectional()

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_participants_sessions(
            self.subjects, self.sessions, folder_participants_tsv
        )

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                "List available in %s"
                % os.path.join(folder_participants_tsv, "participants.tsv")
            )
            cprint("The pipeline will last approximately a few hours per image.")

    def build_input_node_longitudinal(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.utils.input_files as input_files
        from clinica.iotools.utils.data_handling import (
            check_relative_volume_location_in_world_coordinate_system,
        )
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.inputs import clinica_file_reader

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
            synchronize=True,
        )

        all_errors = []
        try:
            read_parameters_node.inputs.pet, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.bids_directory,
                input_files.bids_pet_nii(self.parameters["acq_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.orig_nu, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_ORIG_NU,
            )

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_SURF_R,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_SURF_L,
            )

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_DESTRIEUX_PARC_L,
            )

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_DESTRIEUX_PARC_R,
            )

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_DESIKAN_PARC_L,
            )

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_LONG_DESIKAN_PARC_R,
            )

        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system(
            "T1w-MRI (orig_nu.mgz)",
            read_parameters_node.inputs.orig_nu,
            self.parameters["acq_label"] + " PET",
            read_parameters_node.inputs.pet,
            self.bids_directory,
            self.parameters["acq_label"],
            skip_question=self.parameters["skip_question"],
        )

        # fmt: off
        self.connect(
            [
                (read_parameters_node, self.input_node, [("pet", "pet"),
                                                         ("orig_nu", "orig_nu"),
                                                         ("white_surface_left", "white_surface_left"),
                                                         ("white_surface_right", "white_surface_right"),
                                                         ("destrieux_left", "destrieux_left"),
                                                         ("destrieux_right", "destrieux_right"),
                                                         ("desikan_left", "desikan_left"),
                                                         ("desikan_right", "desikan_right")])
            ]
        )
        # fmt: on

    def build_input_node_cross_sectional(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.utils.input_files as input_files
        from clinica.iotools.utils.data_handling import (
            check_relative_volume_location_in_world_coordinate_system,
        )
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.inputs import clinica_file_reader

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
            synchronize=True,
        )

        all_errors = []
        try:

            read_parameters_node.inputs.pet, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.bids_directory,
                input_files.bids_pet_nii(self.parameters["acq_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.orig_nu, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_ORIG_NU,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_WM_SURF_R,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_WM_SURF_L,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_DESTRIEUX_PARC_L,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_DESTRIEUX_PARC_R,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_left, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_DESIKAN_PARC_L,
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_right, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                input_files.T1_FS_DESIKAN_PARC_R,
            )
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system(
            "T1w-MRI (orig_nu.mgz)",
            read_parameters_node.inputs.orig_nu,
            self.parameters["acq_label"] + " PET",
            read_parameters_node.inputs.pet,
            self.bids_directory,
            self.parameters["acq_label"],
            skip_question=self.parameters["skip_question"],
        )

        # fmt: off
        self.connect(
            [
                (read_parameters_node, self.input_node, [("pet", "pet"),
                                                         ("orig_nu", "orig_nu"),
                                                         ("white_surface_left", "white_surface_left"),
                                                         ("white_surface_right", "white_surface_right"),
                                                         ("destrieux_left", "destrieux_left"),
                                                         ("destrieux_right", "destrieux_right"),
                                                         ("desikan_left", "desikan_left"),
                                                         ("desikan_right", "desikan_right")])
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.

        The function get_wf constructs a pipeline for one subject (in pet_surface_utils.py) and runs it.
        We use iterables to give to the node all the files and information needed.
        """
        # TODO(@arnaud.marcoux): Convert it to a Node with iterables + MapNodes.
        #   I'm experimenting something to avoid the "MapNode of MapNode" case
        #   with iterables. I'll try to apply it on the tractography pipeline.
        #   Check it out to get inspiration from it when it's ready.

        import os

        import nipype.interfaces.utility as niu
        import nipype.pipeline.engine as npe

        import clinica.pipelines.pet_surface.pet_surface_utils as utils
        from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

        full_pipe = npe.MapNode(
            niu.Function(
                input_names=[
                    "subject_id",
                    "session_id",
                    "caps_dir",
                    "pvc_psf_tsv",
                    "pet",
                    "orig_nu",
                    "white_surface_left",
                    "white_surface_right",
                    "working_directory_subjects",
                    "acq_label",
                    "csv_segmentation",
                    "suvr_reference_region",
                    "matscript_folder_inverse_deformation",
                    "desikan_left",
                    "desikan_right",
                    "destrieux_left",
                    "destrieux_right",
                    "spm_standalone_is_available",
                    "is_longitudinal",
                ],
                output_names=[],
                function=utils.get_wf,
            ),
            name="full_pipeline_mapnode",
            iterfield=[
                "subject_id",
                "session_id",
                "pet",
                "orig_nu",
                "white_surface_left",
                "white_surface_right",
                "desikan_left",
                "desikan_right",
                "destrieux_left",
                "destrieux_right",
            ],
        )

        full_pipe.inputs.subject_id = self.subjects
        full_pipe.inputs.session_id = self.sessions
        full_pipe.inputs.caps_dir = self.caps_directory
        full_pipe.inputs.pvc_psf_tsv = self.parameters["pvc_psf_tsv"]
        full_pipe.inputs.working_directory_subjects = os.path.join(
            self.base_dir, self.name
        )
        full_pipe.inputs.acq_label = self.parameters["acq_label"]
        full_pipe.inputs.suvr_reference_region = self.parameters[
            "suvr_reference_region"
        ]
        full_pipe.inputs.csv_segmentation = os.path.abspath(
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "..",
                "..",
                "resources",
                "label_conversion_gtmsegmentation.csv",
            )
        )
        full_pipe.inputs.matscript_folder_inverse_deformation = os.path.abspath(
            os.path.dirname(os.path.realpath(__file__))
        )
        full_pipe.inputs.is_longitudinal = self.parameters["longitudinal"]

        # This section of code determines whether to use SPM standalone or not
        if spm_standalone_is_available():
            use_spm_standalone()
            full_pipe.inputs.spm_standalone_is_available = True
        else:
            full_pipe.inputs.spm_standalone_is_available = False

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, full_pipe, [("pet", "pet"),
                                              ("white_surface_left", "white_surface_left"),
                                              ("white_surface_right", "white_surface_right"),
                                              ("orig_nu", "orig_nu"),
                                              ("destrieux_left", "destrieux_left"),
                                              ("destrieux_right", "destrieux_right"),
                                              ("desikan_left", "desikan_left"),
                                              ("desikan_right", "desikan_right")])
            ]
        )
        # fmt: off
