from nipype import config

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class PETVolume(cpe.Pipeline):
    """PETVolume - Volume-based processing of PET images using SPM.

    Returns:
        A clinica pipeline object containing the PETVolume pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.atlas import PET_VOLUME_ATLASES
        from clinica.utils.group import check_group_label

        self.parameters.setdefault("group_label", None)
        check_group_label(self.parameters["group_label"])
        if "acq_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory acq_label key in pipeline parameter.")

        self.parameters.setdefault("pvc_psf_tsv", None)
        self.parameters.setdefault("mask_tissues", [1, 2, 3])
        self.parameters.setdefault("mask_threshold", 0.3)
        self.parameters.setdefault("pvc_mask_tissues", [1, 2, 3])
        self.parameters.setdefault("smooth", [8])
        self.parameters.setdefault("atlases", PET_VOLUME_ATLASES)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return [
            "pet_image",
            "t1_image_native",
            "mask_tissues",
            "pvc_mask_tissues",
            "psf",
            "flow_fields",
            "dartel_template",
            "reference_mask",
        ]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return [
            "pet_t1_native",
            "pet_mni",
            "pet_suvr",
            "binary_mask",
            "pet_suvr_masked",
            "pet_suvr_masked_smoothed",
            "pet_pvc",
            "pet_pvc_mni",
            "pet_pvc_suvr",
            "pet_pvc_suvr_masked",
            "pet_pvc_suvr_masked_smoothed",
            "atlas_statistics",
            "pvc_atlas_statistics",
        ]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os
        from os.path import exists, join

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.iotools.utils.data_handling import (
            check_relative_volume_location_in_world_coordinate_system,
        )
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.input_files import (
            T1W_NII,
            bids_pet_nii,
            t1_volume_deformation_to_template,
            t1_volume_final_group_template,
            t1_volume_native_tpm,
            t1_volume_native_tpm_in_mni,
        )
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        from clinica.utils.pet import get_suvr_mask, read_psf_information
        from clinica.utils.stream import cprint
        from clinica.utils.ux import (
            print_groups_in_caps_directory,
            print_images_to_process,
        )

        # Check that group already exists
        if not exists(
            join(
                self.caps_directory, "groups", f"group-{self.parameters['group_label']}"
            )
        ):
            print_groups_in_caps_directory(self.caps_directory)
            raise ClinicaException(
                f"Group {self.parameters['group_label']} does not exist. "
                "Did you run t1-volume or t1-volume-create-dartel pipeline?"
            )

        # Tissues DataGrabber
        # ====================
        all_errors = []

        # Grab reference mask
        reference_mask_file = str(
            get_suvr_mask(self.parameters["suvr_reference_region"])
        )

        # PET from BIDS directory
        try:
            pet_bids, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.bids_directory,
                bids_pet_nii(self.parameters["acq_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        # Native T1w-MRI
        try:
            t1w_bids, _ = clinica_file_reader(
                self.subjects, self.sessions, self.bids_directory, T1W_NII
            )
        except ClinicaException as e:
            all_errors.append(e)

        # mask_tissues
        tissues_input = []
        for tissue_number in self.parameters["mask_tissues"]:
            try:
                current_file, _ = clinica_file_reader(
                    self.subjects,
                    self.sessions,
                    self.caps_directory,
                    t1_volume_native_tpm_in_mni(tissue_number, False),
                )
                tissues_input.append(current_file)
            except ClinicaException as e:
                all_errors.append(e)
        # Tissues_input has a length of len(self.parameters['mask_tissues']). Each of these elements has a size of
        # len(self.subjects). We want the opposite: a list of size len(self.subjects) whose elements have a size of
        # len(self.parameters['mask_tissues']. The trick is to iter on elements with zip(*my_list)
        tissues_input_final = []
        for subject_tissue_list in zip(*tissues_input):
            tissues_input_final.append(subject_tissue_list)
        tissues_input = tissues_input_final

        # Flowfields
        try:
            flowfields_caps, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                t1_volume_deformation_to_template(self.parameters["group_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        # Dartel Template
        try:
            final_template = clinica_group_reader(
                self.caps_directory,
                t1_volume_final_group_template(self.parameters["group_label"]),
            )
        except ClinicaException as e:
            all_errors.append(e)

        if self.parameters["pvc_psf_tsv"] is not None:
            iterables_psf = read_psf_information(
                self.parameters["pvc_psf_tsv"],
                self.subjects,
                self.sessions,
                self.parameters["acq_label"],
            )
            self.parameters["apply_pvc"] = True
        else:
            iterables_psf = [[]] * len(self.subjects)
            self.parameters["apply_pvc"] = False

        if self.parameters["apply_pvc"]:
            # pvc tissues input
            pvc_tissues_input = []
            for tissue_number in self.parameters["pvc_mask_tissues"]:
                try:
                    current_file, _ = clinica_file_reader(
                        self.subjects,
                        self.sessions,
                        self.caps_directory,
                        t1_volume_native_tpm(tissue_number),
                    )
                    pvc_tissues_input.append(current_file)
                except ClinicaException as e:
                    all_errors.append(e)

            if len(all_errors) == 0:
                pvc_tissues_input_final = []
                for subject_tissue_list in zip(*pvc_tissues_input):
                    pvc_tissues_input_final.append(subject_tissue_list)
                pvc_tissues_input = pvc_tissues_input_final
        else:
            pvc_tissues_input = []

        if len(all_errors) > 0:
            error_message = "Clinica faced error(s) while trying to read files in your CAPS/BIDS directories.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system(
            "T1w-MRI",
            t1w_bids,
            self.parameters["acq_label"] + " PET",
            pet_bids,
            self.bids_directory,
            self.parameters["acq_label"],
            skip_question=self.parameters["skip_question"],
        )

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
            cprint("The pipeline will last approximately 10 minutes per image.")

        read_input_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
            iterables=[
                ("pet_image", pet_bids),
                ("t1_image_native", t1w_bids),
                ("mask_tissues", tissues_input),
                ("psf", iterables_psf),
                ("flow_fields", flowfields_caps),
                ("pvc_mask_tissues", pvc_tissues_input),
            ],
            synchronize=True,
        )

        read_input_node.inputs.reference_mask = reference_mask_file
        read_input_node.inputs.dartel_template = final_template

        # fmt: off
        self.connect(
            [
                (read_input_node, self.input_node, [("pet_image", "pet_image"),
                                                    ("t1_image_native", "t1_image_native"),
                                                    ("mask_tissues", "mask_tissues"),
                                                    ("flow_fields", "flow_fields"),
                                                    ("dartel_template", "dartel_template"),
                                                    ("reference_mask", "reference_mask"),
                                                    ("psf", "psf"),
                                                    ("pvc_mask_tissues", "pvc_mask_tissues")])
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import re

        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import zip_nii
        from clinica.utils.nipype import container_from_filename, fix_join

        # Find container path from pet filename
        # =====================================
        container_path = npe.Node(
            nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="container_path",
        )
        container_path.inputs.threshold = self.parameters["mask_threshold"]

        # Writing all images into CAPS
        # ============================
        write_images_node = npe.Node(name="write_caps_node", interface=nio.DataSink())
        write_images_node.inputs.base_directory = str(self.caps_directory)
        write_images_node.inputs.parameterization = False
        write_images_node.inputs.regexp_substitutions = [
            (r"(.*/)pet_t1_native/r(sub-.*)(\.nii(\.gz)?)$", r"\1\2_space-T1w_pet\3"),
            (
                r"(.*/)pet_pvc/pvc-rbv_r(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-T1w_pvc-rbv_pet\3",
            ),
            (
                r"(.*/)pet_mni/wr(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_pet\3",
            ),
            (
                r"(.*/)pet_pvc_mni/wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_pvc-rbv_pet\3",
            ),
            (
                r"(.*/)pet_suvr/suvr_wr(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_pet\3",
            ),
            (
                r"(.*/)pet_pvc_suvr/suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_pvc-rbv_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_pet\3",
            ),
            (
                r"(.*/)pet_suvr_masked/masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_mask-brain_pet\3",
            ),
            (
                r"(.*/)pet_pvc_suvr_masked/masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$",
                r"\1\2_space-Ixi549Space_pvc-rbv_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_mask-brain_pet\3",
            ),
            (
                r"(.*/)pet_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$",
                r"\1\3_space-Ixi549Space_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_mask-brain_\2_pet\4",
            ),
            (
                r"(.*/)pet_pvc_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$",
                r"\1\3_space-Ixi549Space_pvc-rbv_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"_mask-brain_\2_pet\4",
            ),
            (
                r"(.*/)binary_mask/(sub-.*_T1w_).*(space-[a-zA-Z0-9]+).*(_brainmask\.nii(\.gz)?)$",
                r"\1\2\3\4",
            ),
        ]

        # Writing atlas statistics into CAPS
        # ==================================
        write_atlas_node = npe.Node(name="write_atlas_node", interface=nio.DataSink())
        write_atlas_node.inputs.base_directory = str(self.caps_directory)
        write_atlas_node.inputs.parameterization = False
        write_atlas_node.inputs.regexp_substitutions = [
            (
                r"(.*/atlas_statistics/)suvr_wr(sub-.*)(_statistics\.tsv)$",
                r"\1\2"
                + r"_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"\3",
            ),
            (
                r"(.*/)pvc_(atlas_statistics/)suvr_wpvc-rbv_r(sub-.*)(_statistics\.tsv)$",
                r"\1\2\3"
                + r"_pvc-rbv_suvr-"
                + re.escape(self.parameters["suvr_reference_region"])
                + r"\4",
            ),
        ]

        # fmt: off
        self.connect(
            [
                (self.input_node, container_path, [("pet_image", "bids_or_caps_filename")]),
                (container_path, write_images_node, [(("container", fix_join, "pet", "preprocessing", f"group-{self.parameters['group_label']}"), "container")]),
                (self.output_node, write_images_node, [(("pet_t1_native", zip_nii, True), "pet_t1_native"),
                                                       (("pet_mni", zip_nii, True), "pet_mni"),
                                                       (("pet_suvr", zip_nii, True), "pet_suvr"),
                                                       (("binary_mask", zip_nii, True), "binary_mask"),
                                                       (("pet_suvr_masked", zip_nii, True), "pet_suvr_masked"),
                                                       (("pet_suvr_masked_smoothed", zip_nii, True), "pet_suvr_masked_smoothed"),
                                                       (("pet_pvc", zip_nii, True), "pet_pvc"),
                                                       (("pet_pvc_mni", zip_nii, True), "pet_pvc_mni"),
                                                       (("pet_pvc_suvr", zip_nii, True), "pet_pvc_suvr"),
                                                       (("pet_pvc_suvr_masked", zip_nii, True), "pet_pvc_suvr_masked"),
                                                       (("pet_pvc_suvr_masked_smoothed", zip_nii, True), "pet_pvc_suvr_masked_smoothed")]),
                (container_path, write_atlas_node, [(("container", fix_join, "pet", "preprocessing", f"group-{self.parameters['group_label']}"), "container")]),
                (self.output_node, write_atlas_node, [("atlas_statistics", "atlas_statistics"),
                                                      ("pvc_atlas_statistics", "pvc_atlas_statistics")]),
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.spm as spm
        import nipype.interfaces.spm.utils as spmutils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.petpvc import PETPVC

        from clinica.utils.filemanip import unzip_nii
        from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

        from .pet_volume_utils import (
            apply_binary_mask,
            atlas_statistics,
            create_binary_mask,
            create_pvc_mask,
            get_from_list,
            init_input_node,
            normalize_to_reference,
            pet_pvc_name,
        )

        if spm_standalone_is_available():
            use_spm_standalone()

        # Initialize pipeline
        # ===================
        init_node = npe.Node(
            interface=nutil.Function(
                input_names=["pet_nii"],
                output_names=["pet_nii"],
                function=init_input_node,
            ),
            name="init_pipeline",
        )

        # Unzipping
        # =========
        unzip_pet_image = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_pet_image",
        )

        unzip_t1_image_native = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_t1_image_native",
        )

        unzip_flow_fields = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_flow_fields",
        )

        unzip_dartel_template = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_dartel_template",
        )

        unzip_reference_mask = npe.Node(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_reference_mask",
        )

        unzip_mask_tissues = npe.MapNode(
            nutil.Function(
                input_names=["in_file"], output_names=["out_file"], function=unzip_nii
            ),
            name="unzip_mask_tissues",
            iterfield=["in_file"],
        )

        # Coregister PET into T1 native space
        # ===================================
        coreg_pet_t1 = npe.Node(spm.Coregister(), name="coreg_pet_t1")

        # Spatially normalize PET into MNI
        # ================================
        dartel_mni_reg = npe.Node(spm.DARTELNorm2MNI(), name="dartel_mni_reg")
        dartel_mni_reg.inputs.modulate = False
        dartel_mni_reg.inputs.fwhm = 0

        # Reslice reference region mask into PET
        # ======================================
        reslice = npe.Node(spmutils.Reslice(), name="reslice")

        # Normalize PET values according to reference region
        # ==================================================
        norm_to_ref = npe.Node(
            nutil.Function(
                input_names=["pet_image", "region_mask"],
                output_names=["suvr_pet_path"],
                function=normalize_to_reference,
            ),
            name="norm_to_ref",
        )

        # Create binary mask from segmented tissues
        # =========================================
        binary_mask = npe.Node(
            nutil.Function(
                input_names=["tissues", "threshold"],
                output_names=["out_mask"],
                function=create_binary_mask,
            ),
            name="binary_mask",
        )
        binary_mask.inputs.threshold = self.parameters["mask_threshold"]

        # Mask PET image
        # ==============
        apply_mask = npe.Node(
            nutil.Function(
                input_names=["image", "binary_mask"],
                output_names=["masked_image_path"],
                function=apply_binary_mask,
            ),
            name="apply_mask",
        )

        # Smoothing
        # =========
        if self.parameters["smooth"] is not None and len(self.parameters["smooth"]) > 0:
            smoothing_node = npe.MapNode(
                spm.Smooth(), name="smoothing_node", iterfield=["fwhm", "out_prefix"]
            )
            smoothing_node.inputs.fwhm = [[x, x, x] for x in self.parameters["smooth"]]
            smoothing_node.inputs.out_prefix = [
                "fwhm-" + str(x) + "mm_" for x in self.parameters["smooth"]
            ]
            # fmt: off
            self.connect(
                [
                    (apply_mask, smoothing_node, [("masked_image_path", "in_files")]),
                    (smoothing_node, self.output_node, [("smoothed_files", "pet_suvr_masked_smoothed")]),
                ]
            )
            # fmt: on
        else:
            self.output_node.inputs.pet_suvr_masked_smoothed = [[]]

        # Atlas Statistics
        # ================
        atlas_stats_node = npe.MapNode(
            nutil.Function(
                input_names=["in_image", "in_atlas_list"],
                output_names=["atlas_statistics"],
                function=atlas_statistics,
            ),
            name="atlas_stats_node",
            iterfield=["in_image"],
        )
        atlas_stats_node.inputs.in_atlas_list = self.parameters["atlases"]

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, init_node, [("pet_image", "pet_nii")]),
                (init_node, unzip_pet_image, [("pet_nii", "in_file")]),
                (self.input_node, unzip_t1_image_native, [("t1_image_native", "in_file")]),
                (self.input_node, unzip_flow_fields, [("flow_fields", "in_file")]),
                (self.input_node, unzip_dartel_template, [("dartel_template", "in_file")]),
                (self.input_node, unzip_reference_mask, [("reference_mask", "in_file")]),
                (self.input_node, unzip_mask_tissues, [("mask_tissues", "in_file")]),
                (unzip_pet_image, coreg_pet_t1, [("out_file", "source")]),
                (unzip_t1_image_native, coreg_pet_t1, [("out_file", "target")]),
                (unzip_flow_fields, dartel_mni_reg, [("out_file", "flowfield_files")]),
                (unzip_dartel_template, dartel_mni_reg, [("out_file", "template_file")]),
                (unzip_reference_mask, reslice, [("out_file", "in_file")]),
                (unzip_mask_tissues, binary_mask, [("out_file", "tissues")]),
                (coreg_pet_t1, dartel_mni_reg, [("coregistered_source", "apply_to_files")]),
                (dartel_mni_reg, reslice, [("normalized_files", "space_defining")]),
                (dartel_mni_reg, norm_to_ref, [("normalized_files", "pet_image")]),
                (reslice, norm_to_ref, [("out_file", "region_mask")]),
                (norm_to_ref, apply_mask, [("suvr_pet_path", "image")]),
                (binary_mask, apply_mask, [("out_mask", "binary_mask")]),
                (norm_to_ref, atlas_stats_node, [("suvr_pet_path", "in_image")]),
                (coreg_pet_t1, self.output_node, [("coregistered_source", "pet_t1_native")]),
                (dartel_mni_reg, self.output_node, [("normalized_files", "pet_mni")]),
                (norm_to_ref, self.output_node, [("suvr_pet_path", "pet_suvr")]),
                (binary_mask, self.output_node, [("out_mask", "binary_mask")]),
                (apply_mask, self.output_node, [("masked_image_path", "pet_suvr_masked")]),
                (atlas_stats_node, self.output_node, [("atlas_statistics", "atlas_statistics")]),
            ]
        )
        # fmt: on

        # PVC
        # ==========
        if self.parameters["apply_pvc"]:
            # Unzipping
            # =========
            unzip_pvc_mask_tissues = npe.MapNode(
                nutil.Function(
                    input_names=["in_file"],
                    output_names=["out_file"],
                    function=unzip_nii,
                ),
                name="unzip_pvc_mask_tissues",
                iterfield=["in_file"],
            )

            # Creating Mask to use in PVC
            # ===========================
            pvc_mask = npe.Node(
                nutil.Function(
                    input_names=["tissues"],
                    output_names=["out_mask"],
                    function=create_pvc_mask,
                ),
                name="pvc_mask",
            )
            # PET PVC
            # =======
            petpvc = npe.Node(PETPVC(), name="pvc")
            petpvc.inputs.pvc = "RBV"
            petpvc.inputs.out_file = "pvc.nii"

            # Spatially normalize PET into MNI
            # ================================
            dartel_mni_reg_pvc = npe.Node(
                spm.DARTELNorm2MNI(), name="dartel_mni_reg_pvc"
            )
            dartel_mni_reg_pvc.inputs.modulate = False
            dartel_mni_reg_pvc.inputs.fwhm = 0

            # Reslice reference region mask into PET
            # ======================================
            reslice_pvc = npe.Node(spmutils.Reslice(), name="reslice_pvc")

            # Normalize PET values according to reference region
            # ==================================================
            norm_to_ref_pvc = npe.Node(
                nutil.Function(
                    input_names=["pet_image", "region_mask"],
                    output_names=["suvr_pet_path"],
                    function=normalize_to_reference,
                ),
                name="norm_to_ref_pvc",
            )

            # Mask PET image
            # ==============
            apply_mask_pvc = npe.Node(
                nutil.Function(
                    input_names=["image", "binary_mask"],
                    output_names=["masked_image_path"],
                    function=apply_binary_mask,
                ),
                name="apply_mask_pvc",
            )
            # Smoothing
            # =========
            if (
                self.parameters["smooth"] is not None
                and len(self.parameters["smooth"]) > 0
            ):
                smoothing_pvc = npe.MapNode(
                    spm.Smooth(), name="smoothing_pvc", iterfield=["fwhm", "out_prefix"]
                )
                smoothing_pvc.inputs.fwhm = [
                    [x, x, x] for x in self.parameters["smooth"]
                ]
                smoothing_pvc.inputs.out_prefix = [
                    "fwhm-" + str(x) + "mm_" for x in self.parameters["smooth"]
                ]
                # fmt: off
                self.connect(
                    [
                        (apply_mask_pvc, smoothing_pvc, [("masked_image_path", "in_files")]),
                        (smoothing_pvc, self.output_node, [("smoothed_files", "pet_pvc_suvr_masked_smoothed")]),
                    ]
                )
                # fmt: on
            else:
                self.output_node.inputs.pet_pvc_suvr_masked_smoothed = [[]]
            # Atlas Statistics
            # ================
            atlas_stats_pvc = npe.MapNode(
                nutil.Function(
                    input_names=["in_image", "in_atlas_list"],
                    output_names=["atlas_statistics"],
                    function=atlas_statistics,
                ),
                name="atlas_stats_pvc",
                iterfield=["in_image"],
            )
            atlas_stats_pvc.inputs.in_atlas_list = self.parameters["atlases"]

            # Connection
            # ==========
            # fmt: off
            self.connect(
                [
                    (self.input_node, unzip_pvc_mask_tissues, [("pvc_mask_tissues", "in_file")]),
                    (unzip_pvc_mask_tissues, pvc_mask, [("out_file", "tissues")]),
                    (unzip_flow_fields, dartel_mni_reg_pvc, [("out_file", "flowfield_files")]),
                    (unzip_dartel_template, dartel_mni_reg_pvc, [("out_file", "template_file")]),
                    (unzip_reference_mask, reslice_pvc, [("out_file", "in_file")]),
                    (coreg_pet_t1, petpvc, [("coregistered_source", "in_file"),
                                            (("coregistered_source", pet_pvc_name, "RBV"), "out_file")]),
                    (pvc_mask, petpvc, [("out_mask", "mask_file")]),
                    (self.input_node, petpvc, [(("psf", get_from_list, 0), "fwhm_x"),
                                               (("psf", get_from_list, 1), "fwhm_y"),
                                               (("psf", get_from_list, 2), "fwhm_z")]),
                    (petpvc, dartel_mni_reg_pvc, [("out_file", "apply_to_files")]),
                    (dartel_mni_reg_pvc, reslice_pvc, [("normalized_files", "space_defining")]),
                    (dartel_mni_reg_pvc, norm_to_ref_pvc, [("normalized_files", "pet_image")]),
                    (reslice_pvc, norm_to_ref_pvc, [("out_file", "region_mask")]),
                    (norm_to_ref_pvc, apply_mask_pvc, [("suvr_pet_path", "image")]),
                    (binary_mask, apply_mask_pvc, [("out_mask", "binary_mask")]),
                    (norm_to_ref_pvc, atlas_stats_pvc, [("suvr_pet_path", "in_image")]),
                    (petpvc, self.output_node, [("out_file", "pet_pvc")]),
                    (dartel_mni_reg_pvc, self.output_node, [("normalized_files", "pet_pvc_mni")]),
                    (norm_to_ref_pvc, self.output_node, [("suvr_pet_path", "pet_pvc_suvr")]),
                    (apply_mask_pvc, self.output_node, [("masked_image_path", "pet_pvc_suvr_masked")]),
                    (atlas_stats_pvc, self.output_node, [("atlas_statistics", "pvc_atlas_statistics")]),
                ]
            )
            # fmt: on
        else:
            self.output_node.inputs.pet_pvc = [[]]
            self.output_node.inputs.pet_pvc_mni = [[]]
            self.output_node.inputs.pet_pvc_suvr = [[]]
            self.output_node.inputs.pet_pvc_suvr_masked = [[]]
            self.output_node.inputs.pvc_atlas_statistics = [[]]
            self.output_node.inputs.pet_pvc_suvr_masked_smoothed = [[]]
