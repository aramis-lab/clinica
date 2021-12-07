import clinica.pipelines.engine as cpe


class T1VolumeRegisterDartel(cpe.Pipeline):
    """T1VolumeExistingDartel - Reuse existing Dartel template.

    Returns:
        A clinica pipeline object containing the T1VolumeExistingDartel pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.group import check_group_label

        if "group_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory group_label key in pipeline parameter.")
        self.parameters.setdefault("tissues", [1, 2, 3])

        check_group_label(self.parameters["group_label"])

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["dartel_input_images", "dartel_iteration_templates"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["dartel_flow_fields"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from clinica.utils.input_files import (
            t1_volume_dartel_input_tissue,
            t1_volume_i_th_iteration_group_template,
        )
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        from clinica.utils.ux import print_images_to_process

        read_input_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
        )

        all_errors = []

        # Dartel Input Tissues
        # ====================
        d_input = []
        for tissue_number in self.parameters["tissues"]:
            try:
                current_file, _ = clinica_file_reader(
                    self.subjects,
                    self.sessions,
                    self.caps_directory,
                    t1_volume_dartel_input_tissue(tissue_number),
                )
                d_input.append(current_file)
            except ClinicaException as e:
                all_errors.append(e)

        # Dartel Templates
        # ================
        dartel_iter_templates = []
        for i in range(1, 7):
            try:
                current_iter = clinica_group_reader(
                    self.caps_directory,
                    t1_volume_i_th_iteration_group_template(
                        self.parameters["group_label"], i
                    ),
                )

                dartel_iter_templates.append(current_iter)
            except ClinicaException as e:
                all_errors.append(e)

        if len(all_errors) > 0:
            error_message = "Clinica faced error(s) while trying to read files in your CAPS/BIDS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaCAPSError(error_message)

        read_input_node.inputs.dartel_input_images = d_input
        read_input_node.inputs.dartel_iteration_templates = dartel_iter_templates

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)

        # fmt: off
        self.connect(
            [
                (read_input_node, self.input_node, [("dartel_input_images", "dartel_input_images")]),
                (read_input_node, self.input_node, [("dartel_iteration_templates", "dartel_iteration_templates")]),
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import re

        import nipype.interfaces.io as nio
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import zip_nii

        # Writing flowfields into CAPS
        # ============================
        write_flowfields_node = npe.MapNode(
            name="write_flowfields_node",
            iterfield=["container", "flow_fields"],
            interface=nio.DataSink(infields=["flow_fields"]),
        )
        write_flowfields_node.inputs.base_directory = self.caps_directory
        write_flowfields_node.inputs.parameterization = False
        write_flowfields_node.inputs.container = [
            "subjects/"
            + self.subjects[i]
            + "/"
            + self.sessions[i]
            + "/t1/spm/dartel/group-"
            + self.parameters["group_label"]
            for i in range(len(self.subjects))
        ]
        write_flowfields_node.inputs.regexp_substitutions = [
            (r"(.*)_Template(\.nii(\.gz)?)$", r"\1\2"),
            (r"(.*)c1(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-graymatter\3"),
            (r"(.*)c2(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-whitematter\3"),
            (r"(.*)c3(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-csf\3"),
            (r"(.*)c4(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-bone\3"),
            (r"(.*)c5(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-softtissue\3"),
            (r"(.*)c6(sub-.*)(\.nii(\.gz)?)$", r"\1\2_segm-background\3"),
            (r"(.*)r(sub-.*)(\.nii(\.gz)?)$", r"\1\2\3"),
            (r"(.*)_dartelinput(\.nii(\.gz)?)$", r"\1\2"),
            (
                r"(.*)flow_fields/u_(sub-.*)_segm-.*(\.nii(\.gz)?)$",
                r"\1\2_target-"
                + re.escape(self.parameters["group_label"])
                + r"_transformation-forward_deformation\3",
            ),
            (r"trait_added", r""),
        ]

        # fmt: off
        self.connect(
            [
                (self.output_node, write_flowfields_node, [(("dartel_flow_fields", zip_nii, True), "flow_fields")])
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_utils as utils
        from clinica.utils.filemanip import unzip_nii

        # Unzipping
        # =========
        unzip_dartel_input_node = npe.MapNode(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_dartel_input_node",
            iterfield=["in_file"],
        )
        unzip_templates_node = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_templates_node",
        )
        # DARTEL with existing template
        # =============================
        dartel_existing_template = npe.MapNode(
            utils.DARTELExistingTemplate(),
            name="dartel_existing_template",
            iterfield=["image_files"],
        )

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, unzip_dartel_input_node, [("dartel_input_images", "in_file")]),
                (self.input_node, unzip_templates_node, [("dartel_iteration_templates", "in_file")]),
                (unzip_dartel_input_node, dartel_existing_template, [(("out_file", utils.prepare_dartel_input_images), "image_files")]),
                (unzip_templates_node, dartel_existing_template, [(("out_file", utils.create_iteration_parameters, None), "iteration_parameters")]),
                (dartel_existing_template, self.output_node, [("dartel_flow_fields", "dartel_flow_fields")]),
            ]
        )
        # fmt: on
