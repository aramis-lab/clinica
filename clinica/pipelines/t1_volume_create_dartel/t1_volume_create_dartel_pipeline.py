import clinica.pipelines.engine as cpe


class T1VolumeCreateDartel(cpe.Pipeline):
    """T1VolumeCreateDartel - Create new Dartel template.

    Returns:
        A clinica pipeline object containing the T1VolumeCreateDartel pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.group import check_group_label

        if "group_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory group_label key in pipeline parameter.")
        self.parameters.setdefault("dartel_tissues", [1, 2, 3])

        check_group_label(self.parameters["group_label"])

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["dartel_inputs"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["final_template_file", "template_files", "dartel_flow_fields"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os
        import sys

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.input_files import t1_volume_dartel_input_tissue
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (
            print_begin_image,
            print_groups_in_caps_directory,
            print_images_to_process,
        )

        representative_output = os.path.join(
            self.caps_directory,
            "groups",
            f"group-{self.parameters['group_label']}",
            "t1",
            f"group-{self.parameters['group_label']}_template.nii.gz",
        )
        if os.path.exists(representative_output):
            cprint(
                msg=(
                    f"DARTEL template for {self.parameters['group_label']} already exists. "
                    "Currently, Clinica does not propose to overwrite outputs for this pipeline."
                ),
                lvl="warning",
            )
            print_groups_in_caps_directory(self.caps_directory)
            sys.exit(0)

        # Check that there is at least 2 subjects
        if len(self.subjects) <= 1:
            raise ClinicaException(
                "This pipeline needs at least 2 images to create DARTEL "
                f"template but Clinica only found {len(self.subjects)}."
            )

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
        )
        all_errors = []
        d_input = []
        for tissue_number in self.parameters["dartel_tissues"]:
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

        # Raise all errors if some happened
        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise RuntimeError(error_message)

        # d_input is a list of size len(self.parameters['dartel_tissues'])
        #     Each element of this list is a list of size len(self.subjects)
        read_parameters_node.inputs.dartel_inputs = d_input

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                "Computational time for DARTEL creation will depend on the number of images."
            )
            print_begin_image(f"group-{self.parameters['group_label']}")

        # fmt: off
        self.connect(
            [
                (read_parameters_node, self.input_node, [("dartel_inputs", "dartel_input_images")])
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import os.path as op
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

        # Writing templates into CAPS
        # ===========================
        write_template_node = npe.Node(nio.DataSink(), name="write_template_node")
        write_template_node.inputs.parameterization = False
        write_template_node.inputs.base_directory = self.caps_directory
        write_template_node.inputs.container = op.join(
            "groups", f"group-{self.parameters['group_label']}", "t1"
        )
        write_template_node.inputs.regexp_substitutions = [
            (
                r"(.*)final_template_file/.*(\.nii(\.gz)?)$",
                r"\1group-"
                + re.escape(self.parameters["group_label"])
                + r"_template\2",
            ),
            (
                r"(.*)template_files/.*([0-9])(\.nii(\.gz)?)$",
                r"\1group-"
                + re.escape(self.parameters["group_label"])
                + r"_iteration-\2_template\3",
            ),
        ]

        # fmt: off
        self.connect(
            [
                (self.output_node, write_flowfields_node, [(("dartel_flow_fields", zip_nii, True), "flow_fields")]),
                (self.output_node, write_template_node, [(("final_template_file", zip_nii, True), "final_template_file"),
                                                         (("template_files", zip_nii, True), "template_files")],
                ),
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.spm as spm
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import unzip_nii
        from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

        if spm_standalone_is_available():
            use_spm_standalone()

        # Unzipping
        # =========
        unzip_node = npe.MapNode(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_node",
            iterfield=["in_file"],
        )

        # DARTEL template
        # ===============
        dartel_template = npe.Node(spm.DARTEL(), name="dartel_template")

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, unzip_node, [("dartel_input_images", "in_file")]),
                (unzip_node, dartel_template, [("out_file", "image_files")]),
                (dartel_template, self.output_node, [("dartel_flow_fields", "dartel_flow_fields"),
                                                     ("final_template_file", "final_template_file"),
                                                     ("template_files", "template_files")]),
            ]
        )
        # fmt: on
