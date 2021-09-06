# -*- coding: utf-8 -*-


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config

import clinica.pipelines.engine as cpe

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class DeepLearningPrepareData(cpe.Pipeline):
    """Deeplearning prepare data - MRI in nifti format are transformed into
    PyTorch tensors. The transformation is applied to: the whole volume, a
    selection of 3D patches, or slices extracted from the 3D volume. By default
    it uses the cropped version of the MRI (see option "--use_uncropper_image")


    Returns:
        A clinica pipeline object containing the Deeplearning prepare data pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["input_nifti"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["image_id"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        from os import path

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import (
            ClinicaBIDSError,
            ClinicaCAPSError,
            ClinicaException,
        )
        from clinica.utils.input_files import (
            T1W_EXTENSIVE,
            T1W_LINEAR,
            T1W_LINEAR_CROPPED,
            pet_linear_nii,
        )
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        from .deeplearning_prepare_data_utils import check_mask_list

        # Select the correct filetype corresponding to modality
        if self.parameters.get("modality") == "t1-linear":
            if self.parameters.get("use_uncropped_image"):
                FILE_TYPE = T1W_LINEAR
            else:
                FILE_TYPE = T1W_LINEAR_CROPPED
        if self.parameters.get("modality") == "t1-extensive":
            FILE_TYPE = T1W_EXTENSIVE

        if self.parameters.get("modality") == "pet-linear":
            FILE_TYPE = pet_linear_nii(
                self.parameters.get("acq_label"),
                self.parameters.get("suvr_reference_region"),
                self.parameters.get("use_uncropped_image"),
            )
        if self.parameters.get("modality") == "custom":
            FILE_TYPE = {
                "pattern": f"*{self.parameters.get('custom_suffix')}",
                "description": "Custom suffix",
            }

        # Input file:
        try:
            input_files = clinica_file_reader(
                self.subjects, self.sessions, self.caps_directory, FILE_TYPE
            )
        except ClinicaException as e:
            err = (
                "Clinica faced error(s) while trying to read files in your CAPS directory.\n"
                + str(e)
            )
            raise ClinicaBIDSError(err)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint("The pipeline will last approximately 30 seconds per image.")

        if self.parameters.get("extract_method") == "slice":
            self.slice_direction = self.parameters.get("slice_direction")
            self.slice_mode = self.parameters.get("slice_mode")
        else:
            self.slice_direction = "axial"
            self.slice_mode = "rgb"

        if self.parameters.get("extract_method") == "patch":
            self.patch_size = self.parameters.get("patch_size")
            self.stride_size = self.parameters.get("stride_size")
        else:
            self.patch_size = 50
            self.stride_size = 50

        # Load the corresponding masks
        if self.parameters.get("extract_method") == "roi":
            self.roi_list = self.parameters.get("roi_list")

            if self.parameters.get("modality") == "custom":
                self.mask_pattern = self.parameters.get("custom_mask_pattern")
                self.template = self.parameters.get("custom_template")
                if not self.template:
                    raise ValueError(
                        "A custom template must be defined when the modality is set to custom."
                    )
            else:
                self.mask_pattern = None
                from .deeplearning_prepare_data_utils import TEMPLATE_DICT

                self.template = TEMPLATE_DICT[self.parameters.get("modality")]

            self.masks_location = path.join(
                self.caps_directory, "masks", f"tpl-{self.template}"
            )

            if not self.roi_list:
                raise ValueError("A list of regions must be given.")
            else:
                check_mask_list(
                    self.masks_location,
                    self.roi_list,
                    self.mask_pattern,
                    not self.parameters.get("use_uncropped_image"),
                )
        else:
            self.masks_location = ""
            self.mask_pattern = None
            self.roi_list = []

        # The reading node
        # -------------------------
        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("input_nifti", input_files),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )

        self.connect(
            [
                (read_node, self.input_node, [("input_nifti", "input_nifti")]),
            ]
        )

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import DataSink

        from clinica.utils.filemanip import get_subject_id
        from clinica.utils.nipype import container_from_filename, fix_join

        # Write node
        # ----------------------
        write_node = npe.Node(name="WriteCaps", interface=DataSink())
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False

        # Get subject ID node
        # ----------------------
        image_id_node = npe.Node(
            interface=nutil.Function(
                input_names=["bids_or_caps_file"],
                output_names=["image_id"],
                function=get_subject_id,
            ),
            name="ImageID",
        )

        # Find container path from input filename
        # ----------------------
        container_path = npe.Node(
            nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="ContainerPath",
        )

        # fmt: off
        self.connect(
            [
                (self.input_node, image_id_node, [("input_nifti", "bids_or_caps_file")]),
                (self.input_node, container_path, [("input_nifti", "bids_or_caps_filename")]),
                # (image_id_node, write_node, [('image_id', '@image_id')]),
                (image_id_node, write_node, [("image_id", "@image_id")]),
            ]
        )
        # fmt: on

        subfolder = "image_based"
        if self.parameters.get("extract_method") == "slice":
            subfolder = "slice_based"
            # fmt: off
            self.connect(
                [
                    (self.output_node, write_node, [("slices_rgb_output", "@slices_rgb_output")]),
                    (self.output_node, write_node, [("slices_original_output", "@slices_original_output")]),
                ]
            )
            # fmt: on

        elif self.parameters.get("extract_method") == "patch":
            subfolder = "patch_based"
            # fmt: off
            self.connect(
                [
                    (self.output_node, write_node, [("patches_output", "@patches_output")])
                ]
            )
            # fmt: on
        elif self.parameters.get("extract_method") == "roi":
            subfolder = "roi_based"
            # fmt: off
            self.connect(
                [
                    (self.output_node, write_node, [("roi_output", "@roi_output")])
                ]
            )
            # fmt: on
        else:
            # fmt: off
            self.connect(
                [
                    (self.output_node, write_node, [("output_pt_file", "@output_pt_file")])
                ]
            )
            # fmt: on

        mod_subfolder = ""
        if self.parameters.get("modality") == "t1-linear":
            mod_subfolder = "t1_linear"
        if self.parameters.get("modality") == "t1-extensive":
            mod_subfolder = "t1_extensive"
        if self.parameters.get("modality") == "pet-linear":
            mod_subfolder = "pet_linear"
        if self.parameters.get("modality") == "custom":
            mod_subfolder = "custom"

        # fmt: off
        self.connect(
            [
                (container_path, write_node, [
                    (("container", fix_join, "deeplearning_prepare_data", subfolder, mod_subfolder), "container")]),
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .deeplearning_prepare_data_utils import (
            extract_patches,
            extract_roi,
            extract_slices,
            save_as_pt,
        )

        # The processing nodes
        # Node to save input in nii.gz format into pytorch .pt format
        # ----------------------
        save_as_pt = npe.MapNode(
            name="save_as_pt",
            iterfield=["input_img"],
            interface=nutil.Function(
                function=save_as_pt,
                input_names=["input_img"],
                output_names=["output_file"],
            ),
        )

        # Extract slices node (options: 3 directions, mode)
        # ----------------------
        extract_slices_node = npe.MapNode(
            name="extract_slices",
            iterfield=["input_tensor"],
            interface=nutil.Function(
                function=extract_slices,
                input_names=["input_tensor", "slice_direction", "slice_mode"],
                output_names=["output_file_rgb", "output_file_original"],
            ),
        )

        extract_slices_node.inputs.slice_direction = self.slice_direction
        extract_slices_node.inputs.slice_mode = self.slice_mode

        # Extract patches node (options, patch size and stride size)
        # ----------------------
        extract_patches_node = npe.MapNode(
            name="extract_patches",
            iterfield=["input_tensor"],
            interface=nutil.Function(
                function=extract_patches,
                input_names=["input_tensor", "patch_size", "stride_size"],
                output_names=["output_patch"],
            ),
        )

        extract_patches_node.inputs.patch_size = self.patch_size
        extract_patches_node.inputs.stride_size = self.stride_size

        # Extract ROi node
        extract_roi_node = npe.MapNode(
            name="extract_ROI",
            iterfield=["input_tensor"],
            interface=nutil.Function(
                function=extract_roi,
                input_names=[
                    "input_tensor",
                    "masks_location",
                    "mask_pattern",
                    "cropped_input",
                    "roi_list",
                    "uncrop_output",
                ],
                output_names=["output_roi"],
            ),
        )
        extract_roi_node.inputs.masks_location = self.masks_location
        extract_roi_node.inputs.mask_pattern = self.mask_pattern
        extract_roi_node.inputs.cropped_input = not self.parameters.get(
            "use_uncropped_image"
        )
        extract_roi_node.inputs.roi_list = self.roi_list
        extract_roi_node.inputs.uncrop_output = self.parameters.get("roi_uncrop_output")

        # Connections
        # ----------------------
        # fmt: off
        self.connect(
            [
                (self.input_node, save_as_pt, [("input_nifti", "input_img")]),
            ]
        )

        if self.parameters.get("extract_method") == "slice":
            self.connect(
                [
                    (save_as_pt, extract_slices_node, [("output_file", "input_tensor")]),
                    (extract_slices_node, self.output_node, [("output_file_rgb", "slices_rgb_output")]),
                    (extract_slices_node, self.output_node, [("output_file_original", "slices_original_output")]),
                ]
            )
        elif self.parameters.get("extract_method") == "patch":
            self.connect(
                [
                    (save_as_pt, extract_patches_node, [("output_file", "input_tensor")]),
                    (extract_patches_node, self.output_node, [("output_patch", "patches_output")]),
                ]
            )
        elif self.parameters.get("extract_method") == "roi":
            self.connect(
                [
                    (save_as_pt, extract_roi_node, [("output_file", "input_tensor")]),
                    (extract_roi_node, self.output_node, [("output_roi", "roi_output")]),
                ]
            )
        else:
            self.connect(
                [
                    (save_as_pt, self.output_node, [("output_file", "output_pt_file")]),
                ]
            )
        # fmt: on
