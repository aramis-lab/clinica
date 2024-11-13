from typing import List

from clinica.pipelines.engine import GroupPipeline


class T1VolumeParcellation(GroupPipeline):
    """T1VolumeParcellation - Computation of mean GM concentration for a set of regions.

    Returns:
        A clinica pipeline object containing the T1VolumeParcellation pipeline.
    """

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        from clinica.utils.atlas import T1AndPetVolumeAtlasName

        self.parameters.setdefault("atlases", T1AndPetVolumeAtlasName)
        self.parameters.setdefault("modulate", True)

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["file_list", "atlas_list"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return []

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from clinica.utils.input_files import t1_volume_template_tpm_in_mni
        from clinica.utils.inputs import clinica_file_filter, clinica_file_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (
            print_groups_in_caps_directory,
            print_images_to_process,
        )

        if not self.group_directory.exists():
            print_groups_in_caps_directory(self.caps_directory)
            raise ClinicaException(
                f"Group {self.group_label} does not exist. "
                "Did you run t1-volume or t1-volume-create-dartel pipeline?"
            )

        gm_mni, self.subjects, self.sessions = clinica_file_filter(
            self.subjects,
            self.sessions,
            self.caps_directory,
            t1_volume_template_tpm_in_mni(
                group_label=self.group_label,
                tissue_number=1,
                modulation=self.parameters["modulate"],
            ),
        )

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

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        pass

    def _build_core_nodes(self):
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
        datasink.inputs.base_directory = str(self.caps_directory)
        datasink.inputs.parameterization = True
        datasink.inputs.regexp_substitutions = [
            (
                r"(.*)(atlas_statistics)/.*/(sub-(.*)_ses-(.*)_T1.*)$",
                rf"\1/subjects/sub-\4/ses-\5/t1/spm/dartel/{self.group_id}" + r"/\2/\3",
            )
        ]
        self.connect(
            [
                (self.input_node, atlas_stats_node, [("file_list", "in_image")]),
                (self.input_node, atlas_stats_node, [("atlas_list", "atlas_list")]),
                (
                    atlas_stats_node,
                    outputnode,
                    [("atlas_statistics", "atlas_statistics")],
                ),
                (outputnode, datasink, [("atlas_statistics", "atlas_statistics")]),
            ]
        )
