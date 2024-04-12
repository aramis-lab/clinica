__all__ = [
    "init_input_node_task",
    "move_subjects_dir_to_source_dir_task",
    "save_to_caps_task",
    "write_tsv_files_task",
]


def save_to_caps_task(
    source_dir: str, subject_id: str, caps_dir: str, overwrite_caps: bool = False
) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.utils import (
        save_to_caps,
    )

    return save_to_caps(
        Path(source_dir),
        subject_id,
        Path(caps_dir),
        overwrite_caps=overwrite_caps,
    )


def move_subjects_dir_to_source_dir_task(
    subjects_dir: str, source_dir: str, subject_id: str
) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.utils import (
        move_subjects_dir_to_source_dir,
    )
    from clinica.pipelines.anatomical.freesurfer.utils import (
        extract_image_id_from_freesurfer_id,
    )

    return move_subjects_dir_to_source_dir(
        Path(subjects_dir),
        Path(source_dir),
        subject_id,
        image_id="_".join(extract_image_id_from_freesurfer_id(subject_id)),
    )


def init_input_node_task(
    caps_dir: str, participant_id: str, session_id: str, long_id: str, output_dir: str
) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.utils import (
        init_input_node,
    )

    return str(
        init_input_node(
            Path(caps_dir),
            participant_id,
            session_id,
            long_id,
            Path(output_dir),
        )
    )


def write_tsv_files_task(subjects_dir: str, subject_id: str) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.utils import (
        write_tsv_files,
    )

    return write_tsv_files(Path(subjects_dir), subject_id)
