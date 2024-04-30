__all__ = [
    "init_input_node_task",
    "save_to_caps_task",
    "write_tsv_files_task",
]


def init_input_node_task(t1w: str, recon_all_args: str, output_dir: str) -> tuple:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.t1.utils import init_input_node

    image_id, t1w, flags, subjects_dir = init_input_node(
        Path(t1w), recon_all_args, Path(output_dir)
    )

    return image_id, str(t1w), flags, str(subjects_dir)


def save_to_caps_task(
    source_dir: str, image_id: str, caps_dir: str, overwrite_caps: bool = False
) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.t1.utils import save_to_caps

    return save_to_caps(
        Path(source_dir), image_id, Path(caps_dir), overwrite_caps=overwrite_caps
    )


def write_tsv_files_task(subjects_dir: str, image_id: str) -> str:
    """Adapter for Nipype."""
    from pathlib import Path

    from clinica.pipelines.anatomical.freesurfer.t1.utils import write_tsv_files

    return write_tsv_files(Path(subjects_dir), image_id)
