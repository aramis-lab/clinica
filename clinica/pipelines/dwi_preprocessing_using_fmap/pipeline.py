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

        Returns
        -------
        List[str] :
            The list of inputs for the DwiPreprocessingUsingPhaseDiffFMap pipeline namely:
                - "dwi_filename" : The path of the diffusion weighted image in BIDS format.
                - "b_vectors_filename" : The path of the b-vectors file in BIDS format.
                - "b_values_filename": The path of the b-values file in BIDS format.
                - "dwi_json_filename" : The path of the DWI JSON file in BIDS format and containing
                  'TotalReadoutTime' and 'PhaseEncodingDirection' metadata (see BIDS specifications).
                - "fmap_magnitude_filename": The path to the (1st) magnitude image in BIDS format.
                - "fmap_phasediff_filename" : The path of the phase difference image in BIDS format.
                - "fmap_phasediff_json_filename" : The path of the phase difference JSON file in BIDS
                  format and containing 'EchoTime1' & 'EchoTime2' metadata (see BIDS specifications).
        """
        return [
            "dwi_filename",
            "b_vectors_filename",
            "b_values_filename",
            "dwi_json_filename",
            "fmap_magnitude_filename",
            "fmap_phasediff_filename",
            "fmap_phasediff_json_filename",
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
                (name, value)
                for name, value in zip(self.get_input_fields(), list_bids_files)
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (read_node, self.input_node, [(x, x) for x in self.get_input_fields()]),
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
                    "dwi_filename",
                    "dwi_preproc_filename",
                    "b_values_preproc_filename",
                    "b_vectors_preproc_filename",
                    "b0_brain_mask_filename",
                    "calibrated_magnitude_image_filename",
                    "calibrated_field_map_image_filename",
                    "calibrated_smoothed_field_map_image_filename",
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
                (
                    self.input_node,
                    container_path,
                    [("dwi_filename", "bids_or_caps_filename")],
                ),
                (self.input_node, rename_into_caps, [("dwi_filename", "dwi_filename")]),
                (
                    self.output_node,
                    rename_into_caps,
                    [
                        ("preproc_dwi", "dwi_preproc_filename"),
                        ("preproc_bval", "b_values_preproc_filename"),
                        ("preproc_bvec", "b_vectors_preproc_filename"),
                        ("b0_mask", "b0_brain_mask_filename"),
                        ("magnitude_on_b0", "calibrated_magnitude_image_filename"),
                        (
                            "calibrated_fmap_on_b0",
                            "calibrated_field_map_image_filename",
                        ),
                        (
                            "smoothed_fmap_on_b0",
                            "calibrated_smoothed_field_map_image_filename",
                        ),
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
        from clinica.utils.dwi import compute_average_b0_task

        from .utils import init_input_node, print_end_pipeline
        from .workflows import calibrate_and_register_fmap, compute_reference_b0

        init_node = npe.Node(
            interface=nutil.Function(
                input_names=self.get_input_fields(),
                output_names=[
                    "image_id",
                    "dwi_filename",
                    "b_vectors_filename",
                    "b_values_filename",
                    "total_readout_time",
                    "phase_encoding_direction",
                    "fmap_magnitude_filename",
                    "fmap_phasediff_filename",
                    "delta_echo_time",
                ],
                function=init_input_node,
            ),
            name="0-InitNode",
        )

        reference_b0 = compute_reference_b0(
            self.base_dir,
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
            self.base_dir,
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
                input_names=["dwi_filename", "b_value_filename"],
                output_names=["reference_b0"],
                function=compute_average_b0_task,
            ),
            name="5a-ComputeB0Average",
        )
        compute_avg_b0.inputs.b_value_threshold = self.parameters["low_bval"]

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
            (self.input_node, init_node, [(x, x) for x in self.get_input_fields()]),
            (
                init_node,
                reference_b0,
                [
                    (x, f"inputnode.{x}")
                    for x in (
                        "b_values_filename",
                        "b_vectors_filename",
                        "dwi_filename",
                        "total_readout_time",
                        "phase_encoding_direction",
                        "image_id",
                    )
                ],
            ),
            (
                init_node,
                fmap_calibration_and_registration,
                [
                    ("fmap_magnitude_filename", "inputnode.bias_magnitude_fmap"),
                    ("fmap_phasediff_filename", "inputnode.fmap_phasediff"),
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
                    (x, f"inputnode.{x}")
                    for x in (
                        "dwi_filename",
                        "b_values_filename",
                        "b_vectors_filename",
                        "image_id",
                    )
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
            (init_node, bias, [("b_values_filename", "in_bval")]),
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
            (init_node, compute_avg_b0, [("b_values_filename", "b_value_filename")]),
            (bias, compute_avg_b0, [("out_file", "dwi_filename")]),
            # Compute b0 mask on corrected avg b0
            (compute_avg_b0, mask_avg_b0, [("reference_b0", "in_file")]),
            # Print end message
            (init_node, print_end_message, [("image_id", "image_id")]),
            (mask_avg_b0, print_end_message, [("mask_file", "final_file")]),
            # Output node
            (init_node, self.output_node, [("b_values_filename", "preproc_bval")]),
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
