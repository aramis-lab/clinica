"""PET image processing: late-frame extraction, coregistration, and averaging."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import nibabel as nib
import numpy as np
import pandas as pd

try:
    import ants
except ImportError as exc:
    raise ImportError(
        "antspyx is required for coregistration.\n"
        "Install it with:  pip install antspyx"
    ) from exc

from ._nii_reading import find_bids_json, find_bids_nifti, load_frames

if TYPE_CHECKING:
    from ._tracer_config import TracerConfig


def get_late_frame_indices(
    json_path: Path, late_phase_start_s: float
) -> list[int] | None:
    """Return 0-based indices of frames starting at or after the late-phase threshold.

    Returns None when the FrameTimes structure is absent or no late frames exist.
    """
    with open(json_path, "r", errors="ignore") as f:
        metadata = json.load(f)

    ft = metadata.get("FrameTimes")
    if not ft or not isinstance(ft, dict):
        return None

    inner = ft.get("FrameTimes")
    if not inner or not isinstance(inner, dict):
        return None

    values = inner.get("Values")
    if not values:
        return None

    indices = [i for i, v in enumerate(values) if v[0] >= late_phase_start_s]
    return indices if indices else None


def coregister_and_average(frames: list[nib.Nifti1Image]) -> nib.Nifti1Image:
    """Rigid-body coregister all frames to the first frame, then average them.

    When only one frame is provided it is returned unchanged.
    """
    if len(frames) == 1:
        return frames[0]

    fixed_nib = frames[0]
    fixed_ants = ants.from_nibabel(fixed_nib)

    aligned_data = [np.asarray(fixed_nib.dataobj, dtype=np.float32)]
    for frame_nib in frames[1:]:
        moving_ants = ants.from_nibabel(frame_nib)
        result = ants.registration(
            fixed=fixed_ants,
            moving=moving_ants,
            type_of_transform="Affine",
        )
        warped_nib = ants.to_nibabel(result["warpedmovout"])
        aligned_data.append(np.asarray(warped_nib.dataobj, dtype=np.float32))

    mean_data = np.mean(np.stack(aligned_data, axis=-1), axis=-1).astype(np.float32)
    return nib.Nifti1Image(mean_data, fixed_nib.affine, fixed_nib.header)


def process_scan(
    row: pd.Series, bids_dir: Path, output_dir: Path, tracer_cfg: TracerConfig
) -> bool:
    """Process a single scan record: locate, (extract,) coregister, average, save.

    Parameters
    ----------
    row:
        A row from the usable scan DataFrame.  Must contain Subject_ID,
        Session_ID, File_Name, and Action.
    bids_dir:
        BIDS dataset root directory.
    output_dir:
        Directory where ``*_coreg_avg.nii.gz`` files are written.
    tracer_cfg:
        Tracer-specific configuration.

    Returns
    -------
    True on success (including already-existing outputs), False on any failure.
    """
    from clinica.utils.stream import cprint

    subject_id = row["Subject_ID"]
    session_id = row["Session_ID"]
    file_name = row["File_Name"]
    action = row["Action"]
    tag = f"{subject_id}/{session_id}/{file_name}"

    out_name = Path(file_name).stem + "_coreg_avg.nii.gz"
    out_path = output_dir / out_name
    if out_path.exists():
        cprint(f"[EXISTS] Already processed, skipping: {tag}", lvl="debug")
        return True

    json_path = find_bids_json(bids_dir, subject_id, session_id, file_name)
    if json_path is None:
        cprint(f"[SKIP] JSON not found: {tag}", lvl="warning")
        return False

    nii_path = find_bids_nifti(json_path)
    if nii_path is None:
        cprint(f"[SKIP] NIfTI not found: {tag}", lvl="warning")
        return False

    frame_indices = None
    if action == tracer_cfg.extract_action:
        frame_indices = get_late_frame_indices(
            json_path, tracer_cfg.late_phase_start_s
        )
        if not frame_indices:
            cprint(f"[SKIP] No late-phase frames found in JSON: {tag}", lvl="warning")
            return False
        cprint(
            f"[EXTRACT] {tag}: {len(frame_indices)} late frames "
            f"(indices {frame_indices[0]}-{frame_indices[-1]})",
            lvl="info",
        )

    try:
        frames, _ = load_frames(nii_path, frame_indices)
    except Exception as exc:
        cprint(f"[ERROR] Failed to load {nii_path}: {exc}", lvl="error")
        return False

    cprint(f"[COREG] {tag}: {len(frames)} frame(s)", lvl="info")

    try:
        mean_img = coregister_and_average(frames)
    except Exception as exc:
        cprint(f"[ERROR] Coregistration failed for {tag}: {exc}", lvl="error")
        return False

    nib.save(mean_img, str(out_path))
    cprint(f"[DONE] {out_path}", lvl="info")
    return True


def process_all_scans(
    usable_df: pd.DataFrame,
    bids_dir: Path,
    output_dir: Path,
    tracer_cfg: TracerConfig,
) -> tuple[int, int]:
    """Process all scans in *usable_df* whose Action is in the tracer's usable set.

    Returns
    -------
    (n_ok, n_skip)
        n_ok   — number of scans successfully processed (or already done)
        n_skip — number of scans skipped due to missing files or errors
    """
    to_process = usable_df[usable_df["Action"].isin(tracer_cfg.usable_actions)]
    n_ok = n_skip = 0
    for _, row in to_process.iterrows():
        if process_scan(row, bids_dir, output_dir, tracer_cfg):
            n_ok += 1
        else:
            n_skip += 1
    return n_ok, n_skip
