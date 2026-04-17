"""PET image processing: late-frame extraction, coregistration, and averaging."""

import json
from pathlib import Path

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

from ._nii_reading import find_json, find_nifti, load_frames
from ._scan_classification import USABLE_ACTIONS

LATE_PHASE_START_S = 4500.0  # 75 minutes in seconds

# Maps the Action values from step 1 to a human-readable processing label.
ACTION_TO_PROCESSING = {
    "USE FOR SUVR": "Coregister and Average",
    "Extract 75-100min frames for SUVR": "Extract frames then Coregister and Average",
}


def get_late_frame_indices(json_path: Path) -> list[int] | None:
    """Return 0-based indices of frames starting at or after 75 min post-injection.

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

    indices = [i for i, v in enumerate(values) if v[0] >= LATE_PHASE_START_S]
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


def process_scan(row: pd.Series, data_dir: Path, output_dir: Path) -> bool:
    """Process a single scan record: locate, (extract,) coregister, average, save.

    Parameters
    ----------
    row:
        A row from the usable scan DataFrame.  Must contain Subject_ID,
        Session_ID, File_Name, and Action.
    data_dir:
        Root directory of the raw OASIS-3 data.
    output_dir:
        Directory where ``*_AV1451_coreg_avg.nii.gz`` files are written.

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

    out_name = f"{subject_id}_{session_id}_AV1451_coreg_avg.nii.gz"
    out_path = output_dir / out_name
    if out_path.exists():
        cprint(f"[EXISTS] Already processed, skipping: {tag}", lvl="debug")
        return True

    nii_path = find_nifti(subject_id, session_id, file_name, data_dir)
    if nii_path is None:
        cprint(f"[SKIP] NIfTI not found: {tag}", lvl="warning")
        return False

    frame_indices = None
    if action == "Extract 75-100min frames for SUVR":
        json_path = find_json(subject_id, session_id, file_name, data_dir)
        if json_path is None:
            cprint(f"[SKIP] BIDS JSON not found: {tag}", lvl="warning")
            return False
        frame_indices = get_late_frame_indices(json_path)
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
    usable_df: pd.DataFrame, data_dir: Path, output_dir: Path
) -> tuple[int, int]:
    """Process all scans in *usable_df* whose Action is in USABLE_ACTIONS.

    Returns
    -------
    (n_ok, n_skip)
        n_ok   — number of scans successfully processed (or already done)
        n_skip — number of scans skipped due to missing files or errors
    """
    to_process = usable_df[usable_df["Action"].isin(USABLE_ACTIONS)]
    n_ok = n_skip = 0
    for _, row in to_process.iterrows():
        if process_scan(row, data_dir, output_dir):
            n_ok += 1
        else:
            n_skip += 1
    return n_ok, n_skip
