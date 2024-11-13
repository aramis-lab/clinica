from typing import List

from clinica.pipelines.engine import GroupPipeline


class T1VolumeCreateDartel(GroupPipeline):
    """T1VolumeCreateDartel - Create new Dartel template.

    Returns:
        A clinica pipeline object containing the T1VolumeCreateDartel pipeline.
    """

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        self.parameters.setdefault("dartel_tissues", [1, 2, 3])

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["dartel_inputs"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return ["final_template_file", "template_files", "dartel_flow_fields"]

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import sys

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.input_files import get_t1_volume_dartel_input_tissue
        from clinica.utils.inputs import clinica_list_of_files_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (
            print_begin_image,
            print_groups_in_caps_directory,
            print_images_to_process,
        )

        representative_output = (
            self.group_directory / "t1" / f"{self.group_id}_template.nii.gz"
        )
        if representative_output.exists():
            cprint(
                msg=(
                    f"DARTEL template for {self.group_label} already exists. "
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
        patterns = [
            get_t1_volume_dartel_input_tissue(tissue_number)
            for tissue_number in self.parameters["dartel_tissues"]
        ]
        try:
            d_input = clinica_list_of_files_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                patterns,
            )
            # d_input is a list of size len(self.parameters['dartel_tissues'])
            #     Each element of this list is a list of size len(self.subjects)
            read_parameters_node.inputs.dartel_inputs = d_input
        except ClinicaException as e:
            raise RuntimeError(e)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                "Computational time for DARTEL creation will depend on the number of images."
            )
            print_begin_image(str(self.group_id))

        self.connect(
            [
                (
                    read_parameters_node,
                    self.input_node,
                    [("dartel_inputs", "dartel_input_images")],
                )
            ]
        )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import os.path as op
        import re

        import nipype.interfaces.io as nio
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import zip_nii

        write_flowfields_node = npe.MapNode(
            name="write_flowfields_node",
            iterfield=["container", "flow_fields"],
            interface=nio.DataSink(infields=["flow_fields"]),
        )
        write_flowfields_node.inputs.base_directory = str(self.caps_directory)
        write_flowfields_node.inputs.parameterization = False
        write_flowfields_node.inputs.container = [
            "subjects/"
            + self.subjects[i]
            + "/"
            + self.sessions[i]
            + f"/t1/spm/dartel/{self.group_id}"
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
                + re.escape(str(self.group_label))
                + r"_transformation-forward_deformation\3",
            ),
            (r"trait_added", r""),
        ]

        write_template_node = npe.Node(nio.DataSink(), name="write_template_node")
        write_template_node.inputs.parameterization = False
        write_template_node.inputs.base_directory = str(self.caps_directory)
        write_template_node.inputs.container = op.join(
            "groups", str(self.group_id), "t1"
        )
        write_template_node.inputs.regexp_substitutions = [
            (
                r"(.*)final_template_file/.*(\.nii(\.gz)?)$",
                r"\1group-" + re.escape(str(self.group_label)) + r"_template\2",
            ),
            (
                r"(.*)template_files/.*([0-9])(\.nii(\.gz)?)$",
                r"\1group-"
                + re.escape(str(self.group_label))
                + r"_iteration-\2_template\3",
            ),
        ]

        self.connect(
            [
                (
                    self.output_node,
                    write_flowfields_node,
                    [(("dartel_flow_fields", zip_nii, True), "flow_fields")],
                ),
                (
                    self.output_node,
                    write_template_node,
                    [
                        (("final_template_file", zip_nii, True), "final_template_file"),
                        (("template_files", zip_nii, True), "template_files"),
                    ],
                ),
            ]
        )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.spm as spm
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import unzip_nii
        from clinica.utils.spm import use_spm_standalone_if_available

        use_spm_standalone_if_available()

        unzip_node = npe.MapNode(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_node",
            iterfield=["in_file"],
        )

        dartel_template = npe.Node(spm.DARTEL(), name="dartel_template")

        self.connect(
            [
                (self.input_node, unzip_node, [("dartel_input_images", "in_file")]),
                (unzip_node, dartel_template, [("out_file", "image_files")]),
                (
                    dartel_template,
                    self.output_node,
                    [
                        ("dartel_flow_fields", "dartel_flow_fields"),
                        ("final_template_file", "final_template_file"),
                        ("template_files", "template_files"),
                    ],
                ),
            ]
        )
