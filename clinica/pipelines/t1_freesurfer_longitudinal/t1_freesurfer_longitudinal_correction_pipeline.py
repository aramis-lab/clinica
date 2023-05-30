import clinica.pipelines.engine as cpe


class T1FreeSurferLongitudinalCorrection(cpe.Pipeline):
    """FreeSurfer Longitudinal correction class.

    Returns:
        A clinica pipeline object containing the T1FreeSurferLongitudinalCorrection pipeline
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""

    def check_custom_dependencies(self):
        """Check dependencies that cannot be listed in `info.json`."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the T1FreeSurferLongitudinalCorrection pipeline is:
                * participant_id (str): Participant ID (e.g. "sub-CLNC01")
                * session_id (str): Session ID associated to `participant_id` (e.g. "ses-M000")
                * long_id (str): Longitudinal ID associated to `participant_id` (e.g. "long-M000" or "long-M000M018")

        Returns:
            A list of (string) input fields name.
        """
        return ["participant_id", "session_id", "long_id"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Note:
            The list of inputs of the T1FreeSurferLongitudinalCorrection pipeline is:
                * subject_id (str): FreeSurfer ID
                    (e.g. "sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018")

        Returns:
            A list of (string) output fields name.
        """
        return ["subject_id"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.input_files import T1_FS_DESTRIEUX, T1_FS_T_DESTRIEUX
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint

        from .longitudinal_utils import (
            extract_subject_session_longitudinal_ids_from_filename,
            grab_image_ids_from_caps_directory,
            read_part_sess_long_ids_from_tsv,
            save_part_sess_long_ids_to_tsv,
        )
        from .t1_freesurfer_longitudinal_correction_utils import get_processed_images

        #  We re-initialize subjects, sessions and add long IDS since the latter is not handled by Clinica
        if self.tsv_file:
            (
                self.subjects,
                self.sessions,
                list_long_id,
            ) = read_part_sess_long_ids_from_tsv(self.tsv_file)
        else:
            (
                self.subjects,
                self.sessions,
                list_long_id,
            ) = grab_image_ids_from_caps_directory(self.caps_directory)
        # TODO: Find a way to distinguish cross-sectional from longitudinal pipelines in Clinica Core

        # Display image(s) already present in CAPS folder
        # ===============================================
        processed_ids = get_processed_images(
            self.caps_directory, self.subjects, self.sessions, list_long_id
        )
        if len(processed_ids) > 0:
            cprint(
                msg=f"Clinica found {len(processed_ids)} images(s) already processed in CAPS directory:",
                lvl="warning",
            )
            for image_id in processed_ids:
                cprint(f"{image_id.replace('_', ' | ')}", lvl="warning")
            if self.overwrite_caps:
                output_folder = "<CAPS>/subjects/<participant_id>/<session_id>/t1/<long_id>/freesurfer_longitudinal/"
                cprint(
                    msg=f"Output folders in {output_folder} will be recreated.",
                    lvl="warning",
                )
            else:
                cprint("Image(s) will be ignored by Clinica.", lvl="warning")
                input_ids = [
                    p_id + "_" + s_id + "_" + l_id
                    for p_id, s_id, l_id in zip(
                        self.subjects, self.sessions, list_long_id
                    )
                ]
                to_process_ids = list(set(input_ids) - set(processed_ids))
                (
                    self.subjects,
                    self.sessions,
                    list_long_id,
                ) = extract_subject_session_longitudinal_ids_from_filename(
                    to_process_ids
                )

        all_errors = []
        try:
            # Check that t1-freesurfer has run on the CAPS directory
            clinica_file_reader(
                self.subjects, self.sessions, self.caps_directory, T1_FS_DESTRIEUX
            )
        except ClinicaException as e:
            all_errors.append(e)

        try:
            # Check that t1-freesurfer-template has run on the CAPS directory
            clinica_file_reader(
                self.subjects, list_long_id, self.caps_directory, T1_FS_T_DESTRIEUX
            )
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = "Clinica faced errors while trying to read files in your CAPS directory.\n"
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_part_sess_long_ids_to_tsv(
            self.subjects, self.sessions, list_long_id, folder_participants_tsv
        )

        def print_images_to_process(
            list_participant_id, list_session_id, list_longitudinal_id
        ):
            cprint(
                f"The pipeline will be run on the following {len(list_participant_id)} image(s):"
            )
            for (p_id, s_id, l_id) in zip(
                list_participant_id, list_session_id, list_longitudinal_id
            ):
                cprint(f"\t{p_id} | {s_id} | {l_id}")

        if len(self.subjects):
            # TODO: Generalize long IDs to the message display
            print_images_to_process(self.subjects, self.sessions, list_long_id)
            cprint(
                "List available in %s"
                % os.path.join(folder_participants_tsv, "participants.tsv")
            )
            cprint("The pipeline will last approximately 10 hours per participant.")

        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("participant_id", self.subjects),
                ("session_id", self.sessions),
                ("long_id", list_long_id),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (read_node, self.input_node, [("participant_id", "participant_id")]),
                (read_node, self.input_node, [("session_id", "session_id")]),
                (read_node, self.input_node, [("long_id", "long_id")]),
            ]
        )

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .t1_freesurfer_longitudinal_correction_utils import save_to_caps

        save_to_caps = npe.Node(
            interface=nutil.Function(
                input_names=["source_dir", "subject_id", "caps_dir", "overwrite_caps"],
                output_names=["image_id"],
                function=save_to_caps,
            ),
            name="SaveToCaps",
        )
        save_to_caps.inputs.source_dir = os.path.join(
            self.base_dir, self.name, "ReconAll"
        )
        save_to_caps.inputs.caps_dir = self.caps_directory
        save_to_caps.inputs.overwrite_caps = self.overwrite_caps

        self.connect(
            [
                (self.output_node, save_to_caps, [("subject_id", "subject_id")]),
            ]
        )

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .t1_freesurfer_longitudinal_correction_utils import (
            init_input_node,
            move_subjects_dir_to_source_dir,
            run_recon_all_long,
            write_tsv_files,
        )

        # Nodes declaration
        # =================
        # Initialize the pipeline
        init_input = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "caps_dir",
                    "participant_id",
                    "session_id",
                    "long_id",
                    "output_dir",
                ],
                output_names=["subjects_dir"],
                function=init_input_node,
            ),
            name="0-InitPipeline",
        )
        init_input.inputs.caps_dir = self.caps_directory
        init_input.inputs.output_dir = os.path.join(
            self.base_dir, self.name, "ReconAll"
        )

        # Run recon-all command
        recon_all = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "subjects_dir",
                    "participant_id",
                    "session_id",
                    "long_id",
                    "directive",
                ],
                output_names=["subject_id"],
                function=run_recon_all_long,
            ),
            name="1-SegmentationReconAll",
        )
        recon_all.inputs.directive = "-all"

        # Generate TSV files containing a summary of the regional statistics
        create_tsv = npe.Node(
            interface=nutil.Function(
                input_names=["subjects_dir", "subject_id"],
                output_names=["subject_id"],
                function=write_tsv_files,
            ),
            name="2-CreateTsvFiles",
        )

        # Move $SUBJECT_DIR to source_dir (1 time point case)
        move_subjects_dir = npe.Node(
            interface=nutil.Function(
                input_names=["subjects_dir", "source_dir", "subject_id"],
                output_names=["subject_id"],
                function=move_subjects_dir_to_source_dir,
            ),
            name="3-MoveSubjectsDir",
        )
        move_subjects_dir.inputs.source_dir = os.path.join(
            self.base_dir, self.name, "ReconAll"
        )

        # Connections
        # ===========
        # fmt: off
        self.connect(
            [
                # Initialize the pipeline
                (self.input_node, init_input, [("participant_id", "participant_id"),
                                               ("session_id", "session_id"),
                                               ("long_id", "long_id")]),
                # Run recon-all command
                (init_input, recon_all, [("subjects_dir", "subjects_dir")]),
                (self.input_node, recon_all, [("participant_id", "participant_id"),
                                              ("session_id", "session_id"),
                                              ("long_id", "long_id")]),
                # Generate TSV files
                (init_input, create_tsv, [("subjects_dir", "subjects_dir")]),
                (recon_all, create_tsv, [("subject_id", "subject_id")]),
                # Move $SUBJECT_DIR to source_dir (macOS case)
                (init_input, move_subjects_dir, [("subjects_dir", "subjects_dir")]),
                (create_tsv, move_subjects_dir, [("subject_id", "subject_id")]),
                # Output node
                (move_subjects_dir, self.output_node, [("subject_id", "subject_id")]),
            ]
        )
