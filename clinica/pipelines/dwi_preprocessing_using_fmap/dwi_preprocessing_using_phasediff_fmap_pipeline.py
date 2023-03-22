from nipype import config

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class DwiPreprocessingUsingPhaseDiffFMap(cpe.Pipeline):
    """DWI Preprocessing using phase difference fieldmap.

    Ideas for improvement:
        - Use promising sdcflows workflows and/or dMRIprep

    Note:
        Some reading regarding the reproducibility of FSL eddy command:
        https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;1ccf038f.1608

    Returns:
        A clinica pipeline object containing the DwiPreprocessingUsingPhaseDiffFMap pipeline.
    """

    @staticmethod
    def get_processed_images(caps_directory, subjects, sessions):
        import os

        from clinica.utils.filemanip import extract_image_ids
        from clinica.utils.input_files import DWI_PREPROC_NII
        from clinica.utils.inputs import clinica_file_reader

        image_ids = []
        if os.path.isdir(caps_directory):
            preproc_files, _ = clinica_file_reader(
                subjects, sessions, caps_directory, DWI_PREPROC_NII, False
            )
            image_ids = extract_image_ids(preproc_files)
        return image_ids

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.stream import cprint

        self.parameters.setdefault("low_bval", 5)
        low_bval = self.parameters["low_bval"]
        if low_bval < 0:
            raise ValueError(
                f"The low_bval is negative ({low_bval}): it should be zero or close to zero."
            )
        if self.parameters["low_bval"] > 100:
            cprint(
                f"The low_bval parameter is {low_bval}: it should be close to zero.",
                lvl="warning",
            )

        self.parameters.setdefault("use_cuda", False)
        self.parameters.setdefault("initrand", False)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        from clinica.utils.check_dependency import is_binary_present
        from clinica.utils.exceptions import ClinicaMissingDependencyError

        if self.parameters["use_cuda"]:
            if not is_binary_present("eddy_cuda"):
                raise ClinicaMissingDependencyError(
                    "[Error] FSL eddy with CUDA was set but Clinica could not find eddy_cuda in your PATH environment. "
                    "Check that  https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#The_eddy_executables is correctly set."
                )

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            List[str]: The list of inputs for the DwiPreprocessingUsingPhaseDiffFMap pipeline namely:
                * dwi: Path of the diffusion weighted image in BIDS format
                * bvec: Path of the bvec file in BIDS format
                * bval: Path of the bval file in BIDS format
                * dwi_json: Path of the DWI JSON file in BIDS format and containing
                    TotalReadoutTime and PhaseEncodingDirection metadata (see BIDS specifications)
                * fmap_magnitude: Path of the (1st) magnitude image in BIDS format
                * fmap_phasediff: Path of the phase difference image in BIDS format
                * fmap_phasediff_json: Path of the phase difference JSON file in BIDS format
                    and containing EchoTime1 & EchoTime2 metadata (see BIDS specifications)
        """
        return [
            "dwi",
            "bvec",
            "bval",
            "dwi_json",
            "fmap_magnitude",
            "fmap_phasediff",
            "fmap_phasediff_json",
        ]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            List[str]: The list of outputs generated by the DwiPreprocessingUsingPhaseDiffFMap pipeline namely:
                * preproc_dwi: Path of the preprocessed DWI
                * preproc_bvec: Path of the preprocessed bvec
                * preproc_bval: Path of the preprocessed bval
                * b0_mask: Path of the b0 brainmask
                * magnitude_on_b0: Path of the smoothed calibrated fmap on b0 space
                * calibrated_fmap_on_b0: Path of the calibrated fmap on b0 space
                * smoothed_fmap_on_b0: Path of the magnitude fmap on b0 space
        """
        return [
            "preproc_dwi",
            "preproc_bvec",
            "preproc_bval",
            "b0_mask",
            "magnitude_on_b0",
            "calibrated_fmap_on_b0",
            "smoothed_fmap_on_b0",
        ]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.input_files import (
            DWI_BVAL,
            DWI_BVEC,
            DWI_JSON,
            DWI_NII,
            FMAP_MAGNITUDE1_NII,
            FMAP_PHASEDIFF_JSON,
            FMAP_PHASEDIFF_NII,
        )
        from clinica.utils.inputs import clinica_list_of_files_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        list_bids_files = clinica_list_of_files_reader(
            self.subjects,
            self.sessions,
            self.bids_directory,
            [
                DWI_NII,
                DWI_BVEC,
                DWI_BVAL,
                DWI_JSON,
                FMAP_MAGNITUDE1_NII,
                FMAP_PHASEDIFF_NII,
                FMAP_PHASEDIFF_JSON,
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
            cprint(
                "Computational time will depend of the number of volumes in your DWI dataset and the use of CUDA."
            )

        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("dwi", list_bids_files[0]),
                ("bvec", list_bids_files[1]),
                ("bval", list_bids_files[2]),
                ("dwi_json", list_bids_files[3]),
                ("fmap_magnitude", list_bids_files[4]),
                ("fmap_phasediff", list_bids_files[5]),
                ("fmap_phasediff_json", list_bids_files[6]),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        # fmt: off
        self.connect(
            [
                (read_node, self.input_node, [("dwi", "dwi"),
                                              ("bvec", "bvec"),
                                              ("bval", "bval"),
                                              ("dwi_json", "dwi_json"),
                                              ("fmap_magnitude", "fmap_magnitude"),
                                              ("fmap_phasediff", "fmap_phasediff"),
                                              ("fmap_phasediff_json", "fmap_phasediff_json")]),
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.nipype import container_from_filename, fix_join

        from .dwi_preprocessing_using_phasediff_fmap_utils import rename_into_caps

        # Find container path from DWI filename
        # =====================================
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
                    "in_bids_dwi",
                    "fname_dwi",
                    "fname_bval",
                    "fname_bvec",
                    "fname_brainmask",
                    "fname_magnitude",
                    "fname_fmap",
                    "fname_smoothed_fmap",
                ],
                output_names=[
                    "out_caps_dwi",
                    "out_caps_bval",
                    "out_caps_bvec",
                    "out_caps_brainmask",
                    "out_caps_magnitude",
                    "out_caps_fmap",
                    "out_caps_smoothed_fmap",
                ],
                function=rename_into_caps,
            ),
            name="rename_into_caps",
        )

        # Writing results into CAPS
        # =========================
        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        # fmt: off
        self.connect(
            [
                (self.input_node, container_path, [("dwi", "bids_or_caps_filename")]),
                (self.input_node, rename_into_caps, [("dwi", "in_bids_dwi")]),
                (self.output_node, rename_into_caps, [("preproc_dwi", "fname_dwi"),
                                                      ("preproc_bval", "fname_bval"),
                                                      ("preproc_bvec", "fname_bvec"),
                                                      ("b0_mask", "fname_brainmask"),
                                                      ("magnitude_on_b0", "fname_magnitude"),
                                                      ("calibrated_fmap_on_b0", "fname_fmap"),
                                                      ("smoothed_fmap_on_b0", "fname_smoothed_fmap")]),
                (container_path, write_results, [(("container", fix_join, "dwi"), "container")]),
                (rename_into_caps, write_results, [("out_caps_dwi", "preprocessing.@preproc_dwi"),
                                                   ("out_caps_bval", "preprocessing.@preproc_bval"),
                                                   ("out_caps_bvec", "preprocessing.@preproc_bvec"),
                                                   ("out_caps_brainmask", "preprocessing.@b0_mask"),
                                                   ("out_caps_magnitude", "preprocessing.@magnitude"),
                                                   ("out_caps_fmap", "preprocessing.@fmap"),
                                                   ("out_caps_smoothed_fmap", "preprocessing.@smoothed_fmap")]),
            ]
        )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.ants as ants
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix3 as mrtrix3
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.utility as niu
        import nipype.pipeline.engine as npe

        from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
            eddy_fsl_pipeline,
        )
        from clinica.utils.dwi import compute_average_b0

        from .dwi_preprocessing_using_phasediff_fmap_utils import (
            init_input_node,
            print_end_pipeline,
        )
        from .dwi_preprocessing_using_phasediff_fmap_workflows import (
            compute_reference_b0,
            prepare_phasediff_fmap,
        )

        # Step 0: Initialization
        # ======================
        # Initialize input parameters and print begin message
        init_node = npe.Node(
            interface=nutil.Function(
                input_names=self.get_input_fields(),
                output_names=[
                    "image_id",
                    "dwi",
                    "bvec",
                    "bval",
                    "total_readout_time",
                    "phase_encoding_direction",
                    "fmap_magnitude",
                    "fmap_phasediff",
                    "delta_echo_time",
                ],
                function=init_input_node,
            ),
            name="0-InitNode",
        )

        reference_b0 = compute_reference_b0(
            low_bval=self.parameters["low_bval"],
            use_cuda=self.parameters["use_cuda"],
            initrand=self.parameters["initrand"],
        )

        # Step 2: Calibrate and register FMap
        # ===================================
        # Bias field correction of the magnitude image
        bias_mag_fmap = npe.Node(
            ants.N4BiasFieldCorrection(dimension=3), name="2a-N4MagnitudeFmap"
        )
        # Brain extraction of the magnitude image
        bet_mag_fmap = npe.Node(
            fsl.BET(frac=0.4, mask=True), name="2b-BetN4MagnitudeFmap"
        )

        # Calibrate FMap
        calibrate_fmap = prepare_phasediff_fmap(name="2c-CalibrateFMap")

        # Register the BET magnitude fmap onto the BET b0
        bet_mag_fmap2b0 = npe.Node(
            interface=fsl.FLIRT(), name="2d-RegistrationBetMagToB0"
        )
        bet_mag_fmap2b0.inputs.dof = 6
        bet_mag_fmap2b0.inputs.output_type = "NIFTI_GZ"

        # Apply the transformation on the calibrated fmap
        fmap2b0 = npe.Node(interface=fsl.ApplyXFM(), name="2e-1-FMapToB0")
        fmap2b0.inputs.output_type = "NIFTI_GZ"

        # Apply the transformation on the magnitude image
        mag_fmap2b0 = fmap2b0.clone("2e-2-MagFMapToB0")

        # Smooth the registered (calibrated) fmap
        smoothing = npe.Node(interface=fsl.maths.IsotropicSmooth(), name="2f-Smoothing")
        smoothing.inputs.sigma = 4.0

        # Step 3: Run FSL eddy
        # ====================
        eddy = eddy_fsl_pipeline(
            low_bval=self.parameters["low_bval"],
            use_cuda=self.parameters["use_cuda"],
            initrand=self.parameters["initrand"],
        )

        # Step 4: Bias correction
        # =======================
        # Use implementation detailed in (Jeurissen et al., 2014)
        bias = npe.Node(mrtrix3.DWIBiasCorrect(use_ants=True), name="4-RemoveBias")

        # Step 5: Final brainmask
        # =======================
        # Compute average b0 on corrected dataset (for brain mask extraction)
        compute_avg_b0 = npe.Node(
            niu.Function(
                input_names=["in_dwi", "in_bval"],
                output_names=["out_b0_average"],
                function=compute_average_b0,
            ),
            name="5a-ComputeB0Average",
        )
        compute_avg_b0.inputs.low_bval = self.parameters["low_bval"]

        # Compute b0 mask on corrected avg b0
        mask_avg_b0 = npe.Node(fsl.BET(mask=True, robust=True), name="5b-MaskB0")

        # Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["image_id", "final_file"], function=print_end_pipeline
            ),
            name="99-WriteEndMessage",
        )

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                # Step 0: Initialization
                # ======================
                # Initialize input parameters and print begin message
                (self.input_node, init_node, [("dwi", "dwi"),
                                              ("bvec", "bvec"),
                                              ("bval", "bval"),
                                              ("dwi_json", "dwi_json"),
                                              ("fmap_magnitude", "fmap_magnitude"),
                                              ("fmap_phasediff", "fmap_phasediff"),
                                              ("fmap_phasediff_json", "fmap_phasediff_json")]),
                # Step 1: Computation of the reference b0 (i.e. average b0 but with EPI distortions)
                # =======================================
                (init_node, reference_b0, [("bval", "inputnode.b_values"),
                                           ("bvec", "inputnode.b_vectors"),
                                           ("dwi", "inputnode.dwi"),
                                           ("total_readout_time", "inputnode.total_readout_time"),
                                           ("phase_encoding_direction", "inputnode.phase_encoding_direction"),
                                           ("image_id", "inputnode.image_id")]),
                # Step 2: Calibrate and register FMap
                # ===================================
                # Bias field correction of the magnitude image
                (init_node, bias_mag_fmap, [("fmap_magnitude", "input_image")]),
                # Brain extraction of the magnitude image
                (bias_mag_fmap, bet_mag_fmap, [("output_image", "in_file")]),
                # Calibration of the FMap
                (bet_mag_fmap, calibrate_fmap, [("mask_file", "input_node.fmap_mask"),
                                                ("out_file", "input_node.fmap_magnitude")]),
                (init_node, calibrate_fmap, [("fmap_phasediff", "input_node.fmap_phasediff"),
                                             ("delta_echo_time", "input_node.delta_echo_time")]),
                # Register the BET magnitude fmap onto the BET b0
                (bet_mag_fmap, bet_mag_fmap2b0, [("out_file", "in_file")]),
                (reference_b0, bet_mag_fmap2b0, [("outputnode.reference_b0", "reference")]),
                # Apply the transformation on the magnitude image
                (bet_mag_fmap2b0, mag_fmap2b0, [("out_matrix_file", "in_matrix_file")]),
                (bias_mag_fmap, mag_fmap2b0, [("output_image", "in_file")]),
                (reference_b0, mag_fmap2b0, [("outputnode.out_file", "reference")]),
                # Apply the transformation on the calibrated fmap
                (bet_mag_fmap2b0, fmap2b0, [("out_matrix_file", "in_matrix_file")]),
                (calibrate_fmap, fmap2b0, [("output_node.calibrated_fmap", "in_file")]),
                (reference_b0, fmap2b0, [("outputnode.out_file", "reference")]),
                # # Smooth the registered (calibrated) fmap
                (fmap2b0, smoothing, [("out_file", "in_file")]),

                # Step 3: Run FSL eddy
                # ====================
                (init_node, eddy, [("dwi", "inputnode.in_file"),
                                   ("bval", "inputnode.in_bval"),
                                   ("bvec", "inputnode.in_bvec"),
                                   ("image_id", "inputnode.image_id")]),
                (smoothing, eddy, [("out_file", "inputnode.field")]),
                (reference_b0, eddy, [("brainmask", "inputnode.in_mask")]),

                # Step 4: Bias correction
                # =======================
                (init_node, bias, [("bval", "in_bval")]),
                (eddy, bias, [("outputnode.out_rotated_bvecs", "in_bvec"),
                              ("outputnode.out_corrected", "in_file")]),
                # Step 5: Final brainmask
                # =======================
                # Compute average b0 on corrected dataset (for brain mask extraction)
                (init_node, compute_avg_b0, [("bval", "in_bval")]),
                (bias, compute_avg_b0, [("out_file", "in_dwi")]),
                # Compute b0 mask on corrected avg b0
                (compute_avg_b0, mask_avg_b0, [("reference_b0", "in_file")]),

                # Print end message
                (init_node, print_end_message, [("image_id", "image_id")]),
                (mask_avg_b0, print_end_message, [("mask_file", "final_file")]),

                # Output node
                (init_node, self.output_node, [("bval", "preproc_bval")]),
                (eddy, self.output_node, [("outputnode.out_rotated_bvecs", "preproc_bvec")]),
                (bias, self.output_node, [("out_file", "preproc_dwi")]),
                (mask_avg_b0, self.output_node, [("mask_file", "b0_mask")]),
                (bet_mag_fmap2b0, self.output_node, [("out_file", "magnitude_on_b0")]),
                (fmap2b0, self.output_node, [("out_file", "calibrated_fmap_on_b0")]),
                (smoothing, self.output_node, [("out_file", "smoothed_fmap_on_b0")]),
            ]
        )
        # fmt: on
