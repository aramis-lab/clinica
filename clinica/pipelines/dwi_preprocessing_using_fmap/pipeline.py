from pathlib import Path
from typing import List

from nipype import config

from clinica.pipelines.engine import DWIPreprocessingPipeline

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class DwiPreprocessingUsingPhaseDiffFMap(DWIPreprocessingPipeline):
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
    def get_processed_images(
        caps_directory: Path, subjects: List[str], sessions: List[str]
    ) -> List[str]:
        from clinica.utils.filemanip import extract_image_ids
        from clinica.utils.input_files import DWI_PREPROC_NII
        from clinica.utils.inputs import clinica_file_reader

        image_ids: List[str] = []
        if caps_directory.is_dir():
            preproc_files, _ = clinica_file_reader(
                subjects, sessions, caps_directory, DWI_PREPROC_NII, False
            )
            image_ids = extract_image_ids(preproc_files)
        return image_ids

    def get_input_fields(self) -> List[str]:
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

    def get_output_fields(self) -> List[str]:
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

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
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
        save_participants_sessions(
            self.subjects, self.sessions, self.base_dir / self.name
        )
        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                f"List available in {self.base_dir / self.name / 'participants.tsv'}"
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
        self.connect(
            [
                (
                    read_node,
                    self.input_node,
                    [
                        ("dwi", "dwi"),
                        ("bvec", "bvec"),
                        ("bval", "bval"),
                        ("dwi_json", "dwi_json"),
                        ("fmap_magnitude", "fmap_magnitude"),
                        ("fmap_phasediff", "fmap_phasediff"),
                        ("fmap_phasediff_json", "fmap_phasediff_json"),
                    ],
                ),
            ]
        )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.nipype import container_from_filename, fix_join

        from .utils import rename_into_caps

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

        write_results = npe.Node(name="write_results", interface=nio.DataSink())
        write_results.inputs.base_directory = str(self.caps_directory)
        write_results.inputs.parameterization = False

        self.connect(
            [
                (self.input_node, container_path, [("dwi", "bids_or_caps_filename")]),
                (self.input_node, rename_into_caps, [("dwi", "in_bids_dwi")]),
                (
                    self.output_node,
                    rename_into_caps,
                    [
                        ("preproc_dwi", "fname_dwi"),
                        ("preproc_bval", "fname_bval"),
                        ("preproc_bvec", "fname_bvec"),
                        ("b0_mask", "fname_brainmask"),
                        ("magnitude_on_b0", "fname_magnitude"),
                        ("calibrated_fmap_on_b0", "fname_fmap"),
                        ("smoothed_fmap_on_b0", "fname_smoothed_fmap"),
                    ],
                ),
                (
                    container_path,
                    write_results,
                    [(("container", fix_join, "dwi"), "container")],
                ),
                (
                    rename_into_caps,
                    write_results,
                    [
                        ("out_caps_dwi", "preprocessing.@preproc_dwi"),
                        ("out_caps_bval", "preprocessing.@preproc_bval"),
                        ("out_caps_bvec", "preprocessing.@preproc_bvec"),
                        ("out_caps_brainmask", "preprocessing.@b0_mask"),
                        ("out_caps_magnitude", "preprocessing.@magnitude"),
                        ("out_caps_fmap", "preprocessing.@fmap"),
                        ("out_caps_smoothed_fmap", "preprocessing.@smoothed_fmap"),
                    ],
                ),
            ]
        )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix3 as mrtrix3
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.pipelines.dwi_preprocessing_using_t1.workflows import (
            eddy_fsl_pipeline,
        )
        from clinica.utils.dwi import compute_average_b0

        from .utils import init_input_node, print_end_pipeline
        from .workflows import calibrate_and_register_fmap, compute_reference_b0

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
            b_value_threshold=self.parameters["low_bval"],
            use_cuda=self.parameters["use_cuda"],
            initrand=self.parameters["initrand"],
        )
        fmap_calibration_and_registration = calibrate_and_register_fmap(
            base_dir=self.base_dir
        )

        # Step 3: Run FSL eddy
        # ====================
        eddy = eddy_fsl_pipeline(
            use_cuda=self.parameters["use_cuda"],
            initrand=self.parameters["initrand"],
            image_id=True,
            field=True,
            compute_mask=False,
        )

        # Step 4: Bias correction
        # =======================
        # Use implementation detailed in (Jeurissen et al., 2014)
        bias = npe.Node(mrtrix3.DWIBiasCorrect(use_ants=True), name="4-RemoveBias")

        # Step 5: Final brainmask
        # =======================
        # Compute average b0 on corrected dataset (for brain mask extraction)
        compute_avg_b0 = npe.Node(
            nutil.Function(
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

        connections = [
            (
                self.input_node,
                init_node,
                [
                    ("dwi", "dwi"),
                    ("bvec", "bvec"),
                    ("bval", "bval"),
                    ("dwi_json", "dwi_json"),
                    ("fmap_magnitude", "fmap_magnitude"),
                    ("fmap_phasediff", "fmap_phasediff"),
                    ("fmap_phasediff_json", "fmap_phasediff_json"),
                ],
            ),
            (
                init_node,
                reference_b0,
                [
                    ("bval", "inputnode.b_values_filename"),
                    ("bvec", "inputnode.b_vectors_filename"),
                    ("dwi", "inputnode.dwi_filename"),
                    ("total_readout_time", "inputnode.total_readout_time"),
                    ("phase_encoding_direction", "inputnode.phase_encoding_direction"),
                    ("image_id", "inputnode.image_id"),
                ],
            ),
            (
                init_node,
                fmap_calibration_and_registration,
                [
                    ("fmap_magnitude", "inputnode.bias_magnitude_fmap"),
                    ("fmap_phasediff", "inputnode.fmap_phasediff"),
                    ("delta_echo_time", "inputnode.delta_echo_time"),
                ],
            ),
            (
                reference_b0,
                fmap_calibration_and_registration,
                [("outputnode.reference_b0", "inputnode.reference_b0")],
            ),
            (
                init_node,
                eddy,
                [
                    ("dwi", "inputnode.dwi_filename"),
                    ("bval", "inputnode.b_values_filename"),
                    ("bvec", "inputnode.b_vectors_filename"),
                    ("image_id", "inputnode.image_id"),
                ],
            ),
            (
                fmap_calibration_and_registration,
                eddy,
                [("outputnode.smooth_calibrated_fmap", "inputnode.field")],
            ),
            (reference_b0, eddy, [("outputnode.brainmask", "inputnode.in_mask")]),
            # Step 4: Bias correction
            # =======================
            (init_node, bias, [("bval", "in_bval")]),
            (
                eddy,
                bias,
                [
                    ("outputnode.out_rotated_bvecs", "in_bvec"),
                    ("outputnode.out_corrected", "in_file"),
                ],
            ),
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
            (
                eddy,
                self.output_node,
                [("outputnode.out_rotated_bvecs", "preproc_bvec")],
            ),
            (bias, self.output_node, [("out_file", "preproc_dwi")]),
            (mask_avg_b0, self.output_node, [("mask_file", "b0_mask")]),
            (
                fmap_calibration_and_registration,
                self.output_node,
                [
                    (
                        "outputnode.bet_magnitude_fmap_registered_onto_b0",
                        "magnitude_on_b0",
                    )
                ],
            ),
            (
                fmap_calibration_and_registration,
                self.output_node,
                [("outputnode.registered_calibrated_fmap", "calibrated_fmap_on_b0")],
            ),
            (
                fmap_calibration_and_registration,
                self.output_node,
                [("outputnode.smooth_calibrated_fmap", "smoothed_fmap_on_b0")],
            ),
        ]

        self.connect(connections)
