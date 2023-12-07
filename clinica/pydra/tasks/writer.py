from os import PathLike
from pathlib import Path
from typing import Optional

from pydra.mark import task


@task
def write_bids_file(
    input_file: PathLike,
    dataset_path: PathLike,
    participant_id: str,
    session_id: Optional[str],
    datatype: str,
    suffix: str,
    entities: dict,
) -> Path:
    source_file = Path(input_file)
    source_ext = source_file.name.split(sep=".", maxsplit=1)

    target_dir = (
        Path(dataset_path)
        / f"sub-{participant_id}"
        / (f"ses-{session_id}" if session_id else "")
        / datatype
    )

    bids_stem = "_".join(
        [f"sub-{participant_id}"]
        + ([f"ses-{session_id}"] if session_id else [])
        + ([f"{ent}-{val}" for ent, val in entities.items()] if entities else [])
        + [suffix]
    )

    target_file = target_dir / ".".join([bids_stem, source_ext])

    # Ensure target directory exists.
    target_dir.mkdir(parents=True, exist_ok=True)

    # Write source file to target directory.
    target_file.write_bytes(source_file.read_bytes())

    return target_file


@task
def write_caps_file(
    input_file: PathLike,
    dataset_path: PathLike,
    participant_id: str,
    session_id: str,
    datatype: str,
    caps_pattern: str,
    entities: dict,
) -> Path:
    source_file = Path(input_file)
    source_ext = source_file.name.split(sep=".", maxsplit=1)

    target_dir = (
        Path(dataset_path)
        / f"sub-{participant_id}"
        / (f"ses-{session_id}" if session_id else "")
        / datatype.format(**entities)
    )

    caps_stem = "_".join(
        [f"sub-{participant_id}"]
        + ([f"ses-{session_id}"] if session_id else [])
        + [caps_pattern.format(**entities)]
    )

    target_file = target_dir / ".".join([caps_stem, source_ext])

    # Ensure target directory exists.
    target_dir.mkdir(parents=True, exist_ok=True)

    # Write source file to target directory.
    target_file.write_bytes(source_file.read_bytes())

    return target_file


@task
def write_freesurfer_dir(
    input_file: PathLike,
    dataset_path: PathLike,
    participant_id: str,
    session_id: str,
    patterns: dict,
    entities: dict,
) -> Path:
    ...
