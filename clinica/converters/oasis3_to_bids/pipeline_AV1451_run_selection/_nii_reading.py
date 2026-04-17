"""NIfTI and BIDS JSON I/O helpers for the AV1451 pipeline."""

from pathlib import Path

import nibabel as nib
import numpy as np


def get_nii_frame_count(json_path: Path) -> int | None:
    """Load the NIfTI co-located with a BIDS JSON sidecar and return its frame count.

    In BIDS the NIfTI sits in the same directory as the JSON with the same stem.
    Returns None when the NIfTI is missing or cannot be read.
    """
    for ext in (".nii.gz", ".nii"):
        nii_path = json_path.with_suffix(ext)
        if nii_path.exists():
            try:
                img = nib.load(str(nii_path))
                shape = img.shape
                return shape[3] if len(shape) > 3 else 1
            except Exception:
                return None
    return None


def find_bids_nifti(json_path: Path) -> Path | None:
    """Return the NIfTI file co-located with *json_path*, or None if not found."""
    for ext in (".nii.gz", ".nii"):
        candidate = json_path.with_suffix(ext)
        if candidate.exists():
            return candidate
    return None


def find_bids_json(
    bids_dir: Path, subject_id: str, session_id: str, file_name: str
) -> Path | None:
    """Return the BIDS JSON sidecar path for a given scan record, or None if not found."""
    candidate = bids_dir / subject_id / session_id / "pet" / file_name
    return candidate if candidate.exists() else None


def sanitize_nifti_header(nib_img: nib.Nifti1Image) -> nib.Nifti1Image:
    """Ensure header zooms match affine-derived voxel spacing.

    Prevents ANTs errors caused by tiny float precision differences between
    the header pixdim and the spacing computed from the affine.
    """
    affine = nib_img.affine
    spacing = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    new_header = nib_img.header.copy()
    zooms = list(new_header.get_zooms())
    zooms[:3] = spacing
    new_header.set_zooms(zooms)
    return nib.Nifti1Image(np.asarray(nib_img.dataobj), affine, new_header)


def load_frames(
    nii_path: Path, frame_indices: list[int] | None = None
) -> tuple[list[nib.Nifti1Image], list[int]]:
    """Load a 4-D NIfTI and return the requested frames as a list of 3-D images.

    Parameters
    ----------
    nii_path:
        Path to the NIfTI file.
    frame_indices:
        0-based indices of frames to load.  None means all frames.

    Returns
    -------
    (frames, used_indices)
        frames       — list of sanitized 3-D Nifti1Image objects
        used_indices — the actual indices loaded
    """
    img = nib.load(str(nii_path))
    data = np.asarray(img.dataobj, dtype=np.float32)

    if data.ndim == 3:
        return [sanitize_nifti_header(img)], [0]

    if frame_indices is None:
        frame_indices = list(range(data.shape[3]))

    affine = img.affine
    header = img.header.copy()
    frames = [
        sanitize_nifti_header(nib.Nifti1Image(data[..., idx], affine, header))
        for idx in frame_indices
    ]
    return frames, frame_indices
