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
    """Task to write a file to a BIDS dataset.

    Examples
    --------
    >>> from shutil import rmtree
    >>> from tempfile import mkdtemp, mkstemp
    >>> test_input_file = Path(mkstemp(suffix=".nii.gz")[1])
    >>> test_output_dir = Path(mkdtemp())
    >>> task = write_bids_file(
    ...     input_file=test_input_file,
    ...     dataset_path=test_output_dir,
    ...     participant_id="P01",
    ...     session_id="M00",
    ...     datatype="anat",
    ...     suffix="T1w",
    ...     entities={"desc": "copy"},
    ... )
    >>> result = task()
    >>> str(next(test_output_dir.rglob("*.nii.gz")).relative_to(test_output_dir))
    'sub-P01/ses-M00/anat/sub-P01_ses-M00_desc-copy_T1w.nii.gz'
    >>> test_input_file.unlink()
    >>> rmtree(test_output_dir)
    """
    source_file = Path(input_file)
    source_ext = source_file.name.split(sep=".", maxsplit=1)[-1]

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
    session_id: Optional[str],
    datatype: str,
    suffix: str,
    entities: dict,
) -> Path:
    """Task to write a file a CAPS dataset.

    Examples
    --------
    >>> from shutil import rmtree
    >>> from tempfile import mkdtemp, mkstemp
    >>> test_input_file = Path(mkstemp(suffix=".mat")[1])
    >>> test_output_dir = Path(mkdtemp())
    >>> task = write_caps_file(
    ...     input_file=test_input_file,
    ...     dataset_path=test_output_dir,
    ...     participant_id="P01",
    ...     session_id="M00",
    ...     datatype="anat",
    ...     suffix="space-{space}_res-{res}_affine",
    ...     entities={"space": "MNI152", "res": "1x1x1"},
    ... )
    >>> result = task()
    >>> str(next(test_output_dir.rglob("*.mat")).relative_to(test_output_dir))
    'sub-P01/ses-M00/anat/sub-P01_ses-M00_space-MNI152_res-1x1x1_affine.mat'
    >>> test_input_file.unlink()
    >>> rmtree(test_output_dir)
    """
    source_file = Path(input_file)
    source_ext = source_file.name.split(sep=".", maxsplit=1)[-1]

    target_dir = (
        Path(dataset_path)
        / f"sub-{participant_id}"
        / (f"ses-{session_id}" if session_id else "")
        / datatype.format(**entities)
    )

    caps_stem = "_".join(
        [f"sub-{participant_id}"]
        + ([f"ses-{session_id}"] if session_id else [])
        + [suffix.format(**entities)]
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
    session_id: Optional[str],
    patterns: dict,
    entities: dict,
) -> Path:
    ...
