from nipype import config

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class DwiDti(cpe.Pipeline):
    """DTI-based processing of DWI datasets.

    Returns:
        A clinica pipeline object containing the DwiDti pipeline.

    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""

    def check_custom_dependencies(self):
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        return ["preproc_dwi", "preproc_bvec", "preproc_bval", "b0_mask"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        output_list_dti = [
            "dti",
            "fa",
            "md",
            "ad",
            "rd",
            "decfa",
            "registered_fa",
            "registered_md",
            "registered_ad",
            "registered_rd",
            "statistics_fa",
            "statistics_md",
            "statistics_ad",
            "statistics_rd",
            "b_spline_transform",
            "affine_matrix",
        ]

        return output_list_dti

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.utils.input_files as input_files
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.inputs import clinica_list_of_files_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        list_caps_files = clinica_list_of_files_reader(
            self.subjects,
            self.sessions,
            self.caps_directory,
            [
                input_files.DWI_PREPROC_NII,
                input_files.DWI_PREPROC_BVEC,
                input_files.DWI_PREPROC_BVAL,
                input_files.DWI_PREPROC_BRAINMASK,
            ],
            raise_exception=True,
        )

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_participants_sessions(
            self.subjects, self.sessions, folder_participants_tsv
        )

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                f"List available in {os.path.join(folder_participants_tsv, 'participants.tsv')}"
            )
            cprint("The pipeline will last approximately 20 minutes per image.")

        read_input_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
            iterables=[
                ("preproc_dwi", list_caps_files[0]),
                ("preproc_bvec", list_caps_files[1]),
                ("preproc_bval", list_caps_files[2]),
                ("b0_mask", list_caps_files[3]),
            ],
            synchronize=True,
        )

        self.connect(
            [
                (read_input_node, self.input_node, [("b0_mask", "b0_mask")]),
                (read_input_node, self.input_node, [("preproc_dwi", "preproc_dwi")]),
                (read_input_node, self.input_node, [("preproc_bval", "preproc_bval")]),
                (read_input_node, self.input_node, [("preproc_bvec", "preproc_bvec")]),
            ]
        )

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.nipype import container_from_filename, fix_join

        from .dwi_dti_utils import rename_into_caps

        # Find container path from filename
        container_path = npe.Node(
            nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="container_path",
        )

        rename_into_caps = npe.Node(
            nutil.Function(
                input_names=[
                    "in_caps_dwi",
                    "in_norm_fa",
                    "in_norm_md",
                    "in_norm_ad",
                    "in_norm_rd",
                    "in_b_spline_transform",
                    "in_affine_matrix",
                ],
                output_names=[
                    "out_caps_fa",
                    "out_caps_md",
                    "out_caps_ad",
                    "out_caps_rd",
                    "out_caps_b_spline_transform",
                    "out_caps_affine_matrix",
                ],
                function=rename_into_caps,
            ),
            name="rename_into_caps",
        )

        # Writing results into CAPS
        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        # fmt: off
        self.connect(
            [
                (self.input_node, container_path, [("preproc_dwi", "bids_or_caps_filename")]),

                (container_path, write_results, [(("container", fix_join, "dwi", "dti_based_processing"), "container")]),
                (self.output_node, write_results, [("dti", "native_space.@dti")]),
                (self.output_node, write_results, [("fa", "native_space.@fa"),
                                                   ("md", "native_space.@md"),
                                                   ("ad", "native_space.@ad"),
                                                   ("rd", "native_space.@rd"),
                                                   ("decfa", "native_space.@decfa")]),

                (self.input_node, rename_into_caps, [("preproc_dwi", "in_caps_dwi")]),
                (self.output_node, rename_into_caps, [("registered_fa", "in_norm_fa"),
                                                      ("registered_md", "in_norm_md"),
                                                      ("registered_ad", "in_norm_ad"),
                                                      ("registered_rd", "in_norm_rd"),
                                                      ("affine_matrix", "in_affine_matrix"),
                                                      ("b_spline_transform", "in_b_spline_transform")]),

                (rename_into_caps, write_results, [("out_caps_fa", "normalized_space.@registered_fa"),
                                                   ("out_caps_md", "normalized_space.@registered_md"),
                                                   ("out_caps_ad", "normalized_space.@registered_ad"),
                                                   ("out_caps_rd", "normalized_space.@registered_rd"),
                                                   ("out_caps_affine_matrix", "normalized_space.@affine_matrix"),
                                                   ("out_caps_b_spline_transform", "normalized_space.@b_spline_transform")]),

                (self.output_node, write_results, [("statistics_fa", "atlas_statistics.@statistics_fa"),
                                                   ("statistics_md", "atlas_statistics.@statistics_md"),
                                                   ("statistics_ad", "atlas_statistics.@statistics_ad"),
                                                   ("statistics_rd", "atlas_statistics.@statistics_rd")])
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import os

        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix as mrtrix
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.ants import ApplyTransforms, RegistrationSynQuick
        from nipype.interfaces.mrtrix3 import TensorMetrics
        from nipype.interfaces.mrtrix.preprocess import DWI2Tensor

        from clinica.utils.check_dependency import check_environment_variable
        from clinica.utils.dwi import extract_bids_identifier_from_filename

        from .dwi_dti_utils import (
            get_ants_transforms,
            get_caps_filenames,
            print_begin_pipeline,
            print_end_pipeline,
            statistics_on_atlases,
        )

        # Nodes creation
        # ==============
        get_bids_identifier = npe.Node(
            interface=nutil.Function(
                input_names=["caps_dwi_filename"],
                output_names=["bids_identifier"],
                function=extract_bids_identifier_from_filename,
            ),
            name="0-Get_BIDS_Identifier",
        )

        get_caps_filenames = npe.Node(
            interface=nutil.Function(
                input_names=["caps_dwi_filename"],
                output_names=[
                    "bids_source",
                    "out_dti",
                    "out_fa",
                    "out_md",
                    "out_ad",
                    "out_rd",
                    "out_evec",
                ],
                function=get_caps_filenames,
            ),
            name="0-CAPS_Filenames",
        )

        convert_gradients = npe.Node(
            interface=mrtrix.FSL2MRTrix(), name="0-Convert_FSL_Gradient"
        )

        dwi_to_dti = npe.Node(interface=DWI2Tensor(), name="1-Compute_DTI")

        dti_to_metrics = npe.Node(interface=TensorMetrics(), name="2-DTI-based_Metrics")

        register_fa = npe.Node(interface=RegistrationSynQuick(), name="3a-Register_FA")
        fsl_dir = check_environment_variable("FSLDIR", "FSL")
        fa_map = os.path.join(
            fsl_dir, "data", "atlases", "JHU", "JHU-ICBM-FA-1mm.nii.gz"
        )
        register_fa.inputs.fixed_image = fa_map

        ants_transforms = npe.Node(
            interface=nutil.Function(
                input_names=["in_affine_transformation", "in_bspline_transformation"],
                output_names=["transforms"],
                function=get_ants_transforms,
            ),
            name="combine_ants_transforms",
        )

        apply_ants_registration = npe.Node(
            interface=ApplyTransforms(), name="apply_ants_registration"
        )
        apply_ants_registration.inputs.dimension = 3
        apply_ants_registration.inputs.input_image_type = 0
        apply_ants_registration.inputs.interpolation = "Linear"
        apply_ants_registration.inputs.reference_image = fa_map

        apply_ants_registration_for_md = apply_ants_registration.clone(
            "3b-Apply_ANTs_Registration_MD"
        )
        apply_ants_registration_for_ad = apply_ants_registration.clone(
            "3b-Apply_ANTs_Registration_AD"
        )
        apply_ants_registration_for_rd = apply_ants_registration.clone(
            "3b-Apply_ANTs_Registration_RD"
        )

        thres_map = npe.Node(
            fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="RemoveNegative"
        )
        thres_norm_fa = thres_map.clone("3c-RemoveNegative_FA")
        thres_norm_md = thres_map.clone("3c-RemoveNegative_MD")
        thres_norm_ad = thres_map.clone("3c-RemoveNegative_AD")
        thres_norm_rd = thres_map.clone("3c-RemoveNegative_RD")

        scalar_analysis = npe.Node(
            interface=nutil.Function(
                input_names=["in_registered_map", "name_map", "prefix_file"],
                output_names=["atlas_statistics_list"],
                function=statistics_on_atlases,
            ),
            name="4-Scalar_Analysis",
        )
        scalar_analysis_fa = scalar_analysis.clone("4-Scalar_Analysis_FA")
        scalar_analysis_fa.inputs.name_map = "FA"
        scalar_analysis_md = scalar_analysis.clone("4-Scalar_Analysis_MD")
        scalar_analysis_md.inputs.name_map = "MD"
        scalar_analysis_ad = scalar_analysis.clone("4-Scalar_Analysis_AD")
        scalar_analysis_ad.inputs.name_map = "AD"
        scalar_analysis_rd = scalar_analysis.clone("4-Scalar_Analysis_RD")
        scalar_analysis_rd.inputs.name_map = "RD"

        thres_map = npe.Node(
            fsl.Threshold(thresh=0.0), iterfield=["in_file"], name="5-Remove_Negative"
        )
        thres_fa = thres_map.clone("5-Remove_Negative_FA")
        thres_md = thres_map.clone("5-Remove_Negative_MD")
        thres_ad = thres_map.clone("5-Remove_Negative_AD")
        thres_rd = thres_map.clone("5-Remove_Negative_RD")

        print_begin_message = npe.Node(
            interface=nutil.Function(
                input_names=["in_bids_or_caps_file"], function=print_begin_pipeline
            ),
            name="Write-Begin_Message",
        )

        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["in_bids_or_caps_file", "final_file_1", "final_file_2"],
                function=print_end_pipeline,
            ),
            name="Write-End_Message",
        )

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, get_caps_filenames, [("preproc_dwi", "caps_dwi_filename")]),
                # Print begin message
                (self.input_node, print_begin_message, [("preproc_dwi", "in_bids_or_caps_file")]),
                # Get BIDS/CAPS identifier from filename
                (self.input_node, get_bids_identifier, [("preproc_dwi", "caps_dwi_filename")]),
                # Convert FSL gradient files (bval/bvec) to MRtrix format
                (self.input_node, convert_gradients, [("preproc_bval", "bval_file"),
                                                      ("preproc_bvec", "bvec_file")]),
                # Computation of the DTI model
                (self.input_node, dwi_to_dti, [("b0_mask", "mask"),
                                               ("preproc_dwi", "in_file")]),
                (convert_gradients, dwi_to_dti, [("encoding_file", "encoding_file")]),
                (get_caps_filenames, dwi_to_dti, [("out_dti", "out_filename")]),
                # Computation of the different metrics from the DTI
                (get_caps_filenames, dti_to_metrics, [("out_fa", "out_fa")]),
                (get_caps_filenames, dti_to_metrics, [("out_md", "out_adc")]),
                (get_caps_filenames, dti_to_metrics, [("out_ad", "out_ad")]),
                (get_caps_filenames, dti_to_metrics, [("out_rd", "out_rd")]),
                (get_caps_filenames, dti_to_metrics, [("out_evec", "out_evec")]),
                (self.input_node, dti_to_metrics, [("b0_mask", "in_mask")]),
                (dwi_to_dti, dti_to_metrics, [("tensor", "in_file")]),
                # Registration of FA-map onto the atlas:
                (dti_to_metrics, register_fa, [("out_fa", "moving_image")]),
                # Apply deformation field on MD, AD & RD:
                (register_fa, ants_transforms, [("out_matrix", "in_affine_transformation")]),
                (register_fa, ants_transforms, [("forward_warp_field", "in_bspline_transformation")]),

                (dti_to_metrics, apply_ants_registration_for_md, [("out_adc", "input_image")]),
                (ants_transforms, apply_ants_registration_for_md, [("transforms", "transforms")]),

                (dti_to_metrics, apply_ants_registration_for_ad, [("out_ad", "input_image")]),
                (ants_transforms, apply_ants_registration_for_ad, [("transforms", "transforms")]),

                (dti_to_metrics, apply_ants_registration_for_rd, [("out_rd", "input_image")]),
                (ants_transforms, apply_ants_registration_for_rd, [("transforms", "transforms")]),
                # Remove negative values from the DTI maps:
                (register_fa, thres_norm_fa, [("warped_image", "in_file")]),
                (apply_ants_registration_for_md, thres_norm_md, [("output_image", "in_file")]),
                (apply_ants_registration_for_rd, thres_norm_rd, [("output_image", "in_file")]),
                (apply_ants_registration_for_ad, thres_norm_ad, [("output_image", "in_file")]),
                # Generate regional TSV files
                (get_bids_identifier, scalar_analysis_fa, [("bids_identifier", "prefix_file")]),
                (thres_norm_fa, scalar_analysis_fa, [("out_file", "in_registered_map")]),
                (get_bids_identifier, scalar_analysis_md, [("bids_identifier", "prefix_file")]),
                (thres_norm_md, scalar_analysis_md, [("out_file", "in_registered_map")]),
                (get_bids_identifier, scalar_analysis_ad, [("bids_identifier", "prefix_file")]),
                (thres_norm_ad, scalar_analysis_ad, [("out_file", "in_registered_map")]),
                (get_bids_identifier, scalar_analysis_rd, [("bids_identifier", "prefix_file")]),
                (thres_norm_rd, scalar_analysis_rd, [("out_file", "in_registered_map")]),
                # Remove negative values from the DTI maps:
                (get_caps_filenames, thres_fa, [("out_fa", "out_file")]),
                (dti_to_metrics, thres_fa, [("out_fa", "in_file")]),

                (get_caps_filenames, thres_md, [("out_md", "out_file")]),
                (dti_to_metrics, thres_md, [("out_adc", "in_file")]),

                (get_caps_filenames, thres_ad, [("out_ad", "out_file")]),
                (dti_to_metrics, thres_ad, [("out_ad", "in_file")]),

                (get_caps_filenames, thres_rd, [("out_rd", "out_file")]),
                (dti_to_metrics, thres_rd, [("out_rd", "in_file")]),

                # Outputnode
                (dwi_to_dti, self.output_node, [("tensor", "dti")]),
                (thres_fa, self.output_node, [("out_file", "fa")]),
                (thres_md, self.output_node, [("out_file", "md")]),
                (thres_ad, self.output_node, [("out_file", "ad")]),
                (thres_rd, self.output_node, [("out_file", "rd")]),
                (dti_to_metrics, self.output_node, [("out_evec", "decfa")]),

                (register_fa, self.output_node, [("out_matrix", "affine_matrix")]),
                (register_fa, self.output_node, [("forward_warp_field", "b_spline_transform")]),

                (thres_norm_fa, self.output_node, [("out_file", "registered_fa")]),
                (thres_norm_md, self.output_node, [("out_file", "registered_md")]),
                (thres_norm_ad, self.output_node, [("out_file", "registered_ad")]),
                (thres_norm_rd, self.output_node, [("out_file", "registered_rd")]),

                (scalar_analysis_fa, self.output_node, [("atlas_statistics_list", "statistics_fa")]),
                (scalar_analysis_md, self.output_node, [("atlas_statistics_list", "statistics_md")]),
                (scalar_analysis_ad, self.output_node, [("atlas_statistics_list", "statistics_ad")]),
                (scalar_analysis_rd, self.output_node, [("atlas_statistics_list", "statistics_rd")]),
                # Print end message
                (self.input_node, print_end_message, [("preproc_dwi", "in_bids_or_caps_file")]),
                (thres_rd, print_end_message, [("out_file", "final_file_1")]),
                (scalar_analysis_rd, print_end_message, [("atlas_statistics_list", "final_file_2")]),
            ]
        )
        # fmt: on
