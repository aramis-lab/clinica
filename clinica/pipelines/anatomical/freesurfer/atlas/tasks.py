__all__ = ["compute_atlas_task", "write_tsv_files_task"]


def compute_atlas_task(
    caps_directory: str, to_process_with_atlases: tuple, path_to_atlas: str
) -> tuple:
    """Adapter for Nipype."""
    from pathlib import Path

    from .utils import compute_atlas

    subject_dir, image_id, atlas = compute_atlas(
        Path(caps_directory),
        to_process_with_atlases,
        Path(path_to_atlas),
    )
    return str(subject_dir), image_id, atlas


def write_tsv_files_task(subject_dir: str, image_id: str, atlas: str) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from .utils import write_tsv_files

    return write_tsv_files(Path(subject_dir), image_id, atlas)
