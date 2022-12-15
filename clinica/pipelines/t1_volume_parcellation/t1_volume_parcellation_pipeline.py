import clinica.pipelines.engine as cpe


class T1VolumeParcellation(cpe.Pipeline):
    """T1VolumeParcellation - Computation of mean GM concentration for a set of regions.

    Returns:
        A clinica pipeline object containing the T1VolumeParcellation pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.atlas import T1_VOLUME_ATLASES
        from clinica.utils.group import check_group_label

        self.parameters.setdefault("group_label", None)
        check_group_label(self.parameters["group_label"])

        self.parameters.setdefault("atlases", T1_VOLUME_ATLASES)
        self.parameters.setdefault("modulate", True)

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["file_list", "atlas_list"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from clinica.utils.input_files import t1_volume_template_tpm_in_mni
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (
            print_groups_in_caps_directory,
            print_images_to_process,
        )

        # Check that group already exists
        if not os.path.exists(
            os.path.join(
                self.caps_directory, "groups", f"group-{self.parameters['group_label']}"
            )
        ):
            print_groups_in_caps_directory(self.caps_directory)
            raise ClinicaException(
                f"Group {self.parameters['group_label']} does not exist. "
                "Did you run t1-volume or t1-volume-create-dartel pipeline?"
            )

        try:
            gm_mni, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                t1_volume_template_tpm_in_mni(
                    group_label=self.parameters["group_label"],
                    tissue_number=1,
                    modulation=self.parameters["modulate"],
                ),
            )
        except ClinicaException as e:
            final_error_str = "Clinica faced error(s) while trying to read files in your CAPS directory.\n"
            final_error_str += str(e)
            raise ClinicaCAPSError(final_error_str)

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
        )
        read_parameters_node.inputs.file_list = gm_mni
        read_parameters_node.inputs.atlas_list = self.parameters["atlases"]

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint("The pipeline will last a few seconds per image.")

        self.connect(
            [
                (read_parameters_node, self.input_node, [("file_list", "file_list")]),
                (read_parameters_node, self.input_node, [("atlas_list", "atlas_list")]),
            ]
        )

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from ..t1_volume_parcellation import (
            t1_volume_parcellation_utils as parcellation_utils,
        )

        atlas_stats_node = npe.MapNode(
            nutil.Function(
                input_names=["in_image", "atlas_list"],
                output_names=["atlas_statistics"],
                function=parcellation_utils.atlas_statistics,
            ),
            name="atlas_stats_node",
            iterfield=["in_image"],
        )
        outputnode = npe.Node(
            nutil.IdentityInterface(fields=["atlas_statistics"]),
            name="outputnode",
            mandatory_inputs=True,
        )

        datasink = npe.Node(nio.DataSink(), name="datasink")

        datasink.inputs.base_directory = self.caps_directory
        datasink.inputs.parameterization = True
        datasink.inputs.regexp_substitutions = [
            (
                r"(.*)(atlas_statistics)/.*/(sub-(.*)_ses-(.*)_T1.*)$",
                r"\1/subjects/sub-\4/ses-\5/t1/spm/dartel/group-"
                + self.parameters["group_label"]
                + r"/\2/\3",
            )
        ]

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, atlas_stats_node, [("file_list", "in_image")]),
                (self.input_node, atlas_stats_node, [("atlas_list", "atlas_list")]),
                (atlas_stats_node, outputnode, [("atlas_statistics", "atlas_statistics")]),
                (outputnode, datasink, [("atlas_statistics", "atlas_statistics")]),
            ]
        )
        # fmt: on
