def remove_nan_from_image_task(image_path: str) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import remove_nan_from_image

    return str(remove_nan_from_image(Path(image_path)))


def perform_gtmseg_task(
    caps_dir: str, subject_id: str, session_id: str, is_longitudinal: bool
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import perform_gtmseg

    return str(perform_gtmseg(Path(caps_dir), subject_id, session_id, is_longitudinal))
