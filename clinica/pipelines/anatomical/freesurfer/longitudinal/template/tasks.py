__all__ = [
    "init_input_node_task",
    "move_subjects_dir_to_source_dir_task",
    "save_to_caps_task",
]


def init_input_node_task(
    caps_dir: str, participant_id: str, list_session_ids: list, output_dir: str
) -> tuple:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.template.utils import (
        init_input_node,
    )

    image_id, subjects_dir, flags = init_input_node(
        Path(caps_dir),
        participant_id,
        list_session_ids,
        Path(output_dir),
    )
    return image_id, str(subjects_dir), flags


def move_subjects_dir_to_source_dir_task(
    subjects_dir: str, source_dir: str, subject_id: str
):
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.utils import (
        move_subjects_dir_to_source_dir,
    )

    return move_subjects_dir_to_source_dir(
        Path(subjects_dir),
        Path(source_dir),
        subject_id,
    )


def save_to_caps_task(
    source_dir: str,
    image_id: str,
    list_session_ids: list,
    caps_dir: str,
    overwrite_caps: bool = False,
):
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.template.utils import (
        save_to_caps,
    )

    return save_to_caps(
        Path(source_dir),
        image_id,
        list_session_ids,
        Path(caps_dir),
        overwrite_caps=overwrite_caps,
    )
