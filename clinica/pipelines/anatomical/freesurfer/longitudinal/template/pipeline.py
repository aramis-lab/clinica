from pathlib import Path
from typing import List

from clinica.pipelines.engine import Pipeline
from clinica.utils.input_files import (
    Parcellation,
    get_t1_freesurfer_segmentation,
    get_t1_freesurfer_template,
)


class T1FreeSurferTemplate(Pipeline):
    """FreeSurfer Longitudinal template class.

    Returns:
        A clinica pipeline object containing the T1FreeSurferTemplate pipeline.
    """

    @staticmethod
    def get_processed_images(
        caps_directory: Path, subjects: List[str], sessions: List[str]
    ) -> List[str]:
        import re

        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.longitudinal import get_long_id
        from clinica.utils.participant import get_unique_subjects

        [list_participant_id, list_list_session_ids] = get_unique_subjects(
            subjects, sessions
        )
        list_long_id = [
            get_long_id(list_session_ids) for list_session_ids in list_list_session_ids
        ]
        image_ids: List[str] = []
        if caps_directory.is_dir():
            t1_freesurfer_files, _ = clinica_file_reader(
                list_participant_id,
                list_long_id,
                caps_directory,
                get_t1_freesurfer_template(Parcellation.DESTRIEUX),
            )
            image_ids = [
                re.search(r"(sub-[a-zA-Z0-9]+)_(long-[a-zA-Z0-9]+)", file).group()
                for file in t1_freesurfer_files
            ]
        return image_ids

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        pass

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that cannot be listed in `info.json`."""
        pass

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Notes
        -----
        The list of inputs of the T1FreeSurferTemplate pipeline is:
            * participant_id (str): Participant ID (e.g. "sub-CLNC01")
            * list_session_ids (List[str]: List of session IDs associated to `participant_id`
                (e.g. ["ses-M000"] or ["ses-M000", "ses-M018", "ses-M036"])

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["participant_id", "list_session_ids"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Notes
        -----
        The list of inputs of the T1FreeSurferTemplate pipeline is:
            * image_id (str): Image ID (e.g. "sub-CLNC01_long-M000M018")

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return ["image_id"]

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.pipelines.anatomical.freesurfer.longitudinal.utils import (
            save_part_sess_long_ids_to_tsv,
        )
        from clinica.utils.filemanip import extract_subjects_sessions_from_filename
        from clinica.utils.inputs import clinica_file_filter
        from clinica.utils.longitudinal import (
            get_long_id,
            get_participants_long_id,
            read_sessions,
        )
        from clinica.utils.participant import (
            get_unique_subjects,
            unique_subjects_sessions_to_subjects_sessions,
        )
        from clinica.utils.stream import cprint

        from .utils import extract_participant_long_ids_from_filename

        output_ids = self.get_processed_images(
            self.caps_directory, self.subjects, self.sessions
        )
        (
            processed_participants,
            processed_long_sessions,
        ) = extract_participant_long_ids_from_filename(output_ids)
        if len(processed_participants) > 0:
            cprint(
                msg=(
                    f"Clinica found {len(processed_participants)} participant(s) "
                    "already processed in CAPS directory:"
                ),
                lvl="warning",
            )
            for p_id, l_id in zip(processed_participants, processed_long_sessions):
                cprint(f"{p_id} | {l_id}", lvl="warning")
            if self.overwrite_caps:
                output_folder = "<CAPS>/subjects/<participant_id>/<long_id>/freesurfer_unbiased_template/"
                cprint(
                    f"Output folders in {output_folder} will be recreated.",
                    lvl="warning",
                )
            else:
                cprint("Participant(s) will be ignored by Clinica.", lvl="warning")
                input_ids = [
                    f"{p_id}_{s_id}" for p_id, s_id in zip(self.subjects, self.sessions)
                ]
                processed_sessions_per_participant = [
                    read_sessions(self.caps_directory, p_id, l_id)
                    for (p_id, l_id) in zip(
                        processed_participants, processed_long_sessions
                    )
                ]
                participants, sessions = unique_subjects_sessions_to_subjects_sessions(
                    processed_participants, processed_sessions_per_participant
                )
                processed_ids = [
                    f"{p_id}_{s_id}" for p_id, s_id in zip(participants, sessions)
                ]
                to_process_ids = list(set(input_ids) - set(processed_ids))
                self.subjects, self.sessions = extract_subjects_sessions_from_filename(
                    to_process_ids
                )
        _, self.subjects, self.sessions = clinica_file_filter(
            self.subjects,
            self.sessions,
            self.caps_directory,
            get_t1_freesurfer_segmentation(Parcellation.DESTRIEUX),
        )
        long_ids = get_participants_long_id(self.subjects, self.sessions)
        save_part_sess_long_ids_to_tsv(
            self.subjects, self.sessions, long_ids, self.base_dir / self.name
        )
        [list_participant_id, list_list_session_ids] = get_unique_subjects(
            self.subjects, self.sessions
        )
        list_long_id = [
            get_long_id(list_session_ids) for list_session_ids in list_list_session_ids
        ]

        def print_images_to_process(
            unique_part_list: list[str],
            per_part_session_list: list[str],
            list_part_long_id: list[str],
        ):
            cprint(
                f"The pipeline will be run on the following {len(unique_part_list)} participant(s):"
            )
            for part_id, list_sess_id, list_id in zip(
                unique_part_list, per_part_session_list, list_part_long_id
            ):
                sessions_participant = ", ".join(s_id for s_id in list_sess_id)
                cprint(f"\t{part_id} | {sessions_participant} | {list_id}")

        if len(self.subjects):
            # TODO: Generalize long IDs to the message display
            print_images_to_process(
                list_participant_id, list_list_session_ids, list_long_id
            )
            cprint(
                f"List available in {self.base_dir / self.name / 'participants.tsv'}"
            )
            cprint("The pipeline will last approximately 10 hours per participant.")

        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("participant_id", list_participant_id),
                ("list_session_ids", list_list_session_ids),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (read_node, self.input_node, [("participant_id", "participant_id")]),
                (
                    read_node,
                    self.input_node,
                    [("list_session_ids", "list_session_ids")],
                ),
            ]
        )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .tasks import save_to_caps_task

        save_to_caps = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "source_dir",
                    "image_id",
                    "list_session_ids",
                    "caps_dir",
                    "overwrite_caps",
                ],
                output_names=["image_id"],
                function=save_to_caps_task,
            ),
            name="SaveToCaps",
        )
        save_to_caps.inputs.source_dir = str(self.base_dir / self.name / "ReconAll")
        save_to_caps.inputs.caps_dir = str(self.caps_directory)
        save_to_caps.inputs.overwrite_caps = False

        self.connect(
            [
                (
                    self.input_node,
                    save_to_caps,
                    [("list_session_ids", "list_session_ids")],
                ),
                (self.output_node, save_to_caps, [("image_id", "image_id")]),
            ]
        )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .tasks import init_input_node_task, move_subjects_dir_to_source_dir_task
        from .utils import run_recon_all_base

        init_input = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "caps_dir",
                    "participant_id",
                    "list_session_ids",
                    "output_dir",
                ],
                output_names=["image_id", "subjects_dir", "flags"],
                function=init_input_node_task,
            ),
            name="0-InitPipeline",
        )
        init_input.inputs.caps_dir = str(self.caps_directory)
        init_input.inputs.output_dir = str(self.base_dir / self.name / "ReconAll")

        recon_all = npe.Node(
            interface=nutil.Function(
                input_names=["subjects_dir", "subject_id", "flags", "directive"],
                output_names=["subject_id"],
                function=run_recon_all_base,
            ),
            name="1-SegmentationReconAll",
        )
        recon_all.inputs.directive = "-all"

        move_subjects_dir = npe.Node(
            interface=nutil.Function(
                input_names=["subjects_dir", "source_dir", "subject_id"],
                output_names=["subject_id"],
                function=move_subjects_dir_to_source_dir_task,
            ),
            name="2-MoveSubjectsDir",
        )
        move_subjects_dir.inputs.source_dir = str(
            self.base_dir / self.name / "ReconAll"
        )

        connections = [
            (self.input_node, init_input, [("participant_id", "participant_id")]),
            (
                self.input_node,
                init_input,
                [("list_session_ids", "list_session_ids")],
            ),
            (
                init_input,
                recon_all,
                [
                    ("subjects_dir", "subjects_dir"),
                    ("image_id", "subject_id"),
                    ("flags", "flags"),
                ],
            ),
            (init_input, move_subjects_dir, [("subjects_dir", "subjects_dir")]),
            (recon_all, move_subjects_dir, [("subject_id", "subject_id")]),
            (move_subjects_dir, self.output_node, [("subject_id", "image_id")]),
        ]
        self.connect(connections)
