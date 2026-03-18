"""
Step 4: Coregister and average tau-PET frames.

For each included scan in the master inventory, co-registers all
(late-phase) time frames to the first frame using rigid-body registration
(ANTs) and averages them into a single 3D NIfTI.

Usage:  python 04_coregister_average.py

Input:  inventory/oasis3_av1451_master.csv
Output: processed/coreg_avg/{Subject}_{Session}_AV1451_coreg_avg.nii.gz

Requires: nibabel, numpy, antspyx  (pip install antspyx)
"""

import json
from pathlib import Path

import nibabel as nib
import numpy as np

try:
    import ants
except ImportError as exc:
    raise ImportError(
        "antspyx is required for coregistration.\n"
        "Install it with:  pip install antspyx"
    ) from exc

import pandas as pd

# --- DEFAULT CONFIGURATION ---
DEFAULT_DATA_DIR   = Path(r"D:\Oasis3\raw")
DEFAULT_MASTER_CSV = Path(r"D:\Oasis3\inventory\oasis3_av1451_master.csv")
DEFAULT_OUTPUT_DIR = Path(r"D:\Oasis3\processed\coreg_avg")
LATE_PHASE_START_S = 4500.0   # 75 minutes in seconds

ACTIONS_TO_PROCESS = {
    "Coregister and Average",
    "Extract frames then Coregister and Average",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_nifti(subject_id: str, session_id: str, file_name: str, data_dir: Path):
    """Find the NIfTI file for a scan entry."""
    session_dir = data_dir / f"{subject_id}_AV1451_{session_id}"
    if not session_dir.is_dir():
        return None

    stem = Path(file_name).stem
    for pet_dir in sorted(session_dir.glob("pet*")):
        nii_dir = pet_dir / "NIFTI"
        for ext in (".nii.gz", ".nii"):
            candidate = nii_dir / (stem + ext)
            if candidate.exists():
                return candidate
    return None


def find_json(subject_id: str, session_id: str, file_name: str, data_dir: Path):
    """Find the BIDS JSON sidecar for a scan entry."""
    session_dir = data_dir / f"{subject_id}_AV1451_{session_id}"
    if not session_dir.is_dir():
        return None

    for pet_dir in sorted(session_dir.glob("pet*")):
        candidate = pet_dir / "BIDS" / file_name
        if candidate.exists():
            return candidate
    return None


def get_late_frame_indices(json_path: Path):
    """Return 0-based indices of frames starting >= 75 min post-injection."""
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


def sanitize_nifti_header(nib_img):
    """Ensure header zooms match affine-derived spacing.

    Prevents ANTs errors caused by tiny float precision differences
    between header pixdim and the spacing computed from the affine.
    """
    affine = nib_img.affine
    spacing = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    new_header = nib_img.header.copy()
    zooms = list(new_header.get_zooms())
    zooms[:3] = spacing
    new_header.set_zooms(zooms)
    return nib.Nifti1Image(np.asarray(nib_img.dataobj), affine, new_header)


def load_frames(nii_path: Path, frame_indices=None):
    """Load a 4-D NIfTI and return a list of 3-D nibabel images."""
    img = nib.load(str(nii_path))
    data = np.asarray(img.dataobj, dtype=np.float32)

    if data.ndim == 3:
        return [sanitize_nifti_header(img)], [0]

    if frame_indices is None:
        frame_indices = list(range(data.shape[3]))

    affine = img.affine
    header = img.header.copy()

    frames = []
    for idx in frame_indices:
        frame_img = nib.Nifti1Image(data[..., idx], affine, header)
        frames.append(sanitize_nifti_header(frame_img))

    return frames, frame_indices


def coregister_and_average(frames):
    """Rigid-body coregister all frames to the first, then average."""
    if len(frames) == 1:
        return frames[0]

    fixed_nib  = frames[0]
    fixed_ants = ants.from_nibabel_nifti(fixed_nib)

    aligned_data = [np.asarray(fixed_nib.dataobj, dtype=np.float32)]

    for frame_nib in frames[1:]:
        moving_ants = ants.from_nibabel_nifti(frame_nib)
        result = ants.registration(
            fixed=fixed_ants,
            moving=moving_ants,
            type_of_transform="Affine",
        )
        warped_nib = ants.to_nibabel_nifti(result["warpedmovout"])
        aligned_data.append(np.asarray(warped_nib.dataobj, dtype=np.float32))

    mean_data = np.mean(np.stack(aligned_data, axis=-1), axis=-1).astype(np.float32)
    return nib.Nifti1Image(mean_data, fixed_nib.affine, fixed_nib.header)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(data_dir=None, inventory_dir=None, output_dir=None):
    data_dir = Path(data_dir) if data_dir else DEFAULT_DATA_DIR
    master_csv = Path(inventory_dir) / "oasis3_av1451_master.csv" if inventory_dir else DEFAULT_MASTER_CSV
    output_dir = Path(output_dir) if output_dir else DEFAULT_OUTPUT_DIR

    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(master_csv)
    to_process = df[df["Processing_Action"].isin(ACTIONS_TO_PROCESS)].copy()

    print(f"Master inventory: {len(df)} total rows")
    print(f"Scans to process: {len(to_process)}\n")

    n_ok   = 0
    n_skip = 0

    for _, row in to_process.iterrows():
        subject_id  = row["Subject_ID"]
        session_id  = row["Session_ID"]
        file_name   = row["File_Name"]
        proc_action = row["Processing_Action"]
        tag = f"{subject_id}/{session_id}/{file_name}"

        # --- Skip if already processed ---
        out_name = f"{subject_id}_{session_id}_AV1451_coreg_avg.nii.gz"
        out_path = output_dir / out_name
        if out_path.exists():
            print(f"  [EXISTS] Already done, skipping: {tag}")
            n_ok += 1
            continue

        # --- Locate NIfTI ---
        nii_path = find_nifti(subject_id, session_id, file_name, data_dir)
        if nii_path is None:
            print(f"  [SKIP] NIfTI not found: {tag}")
            n_skip += 1
            continue

        # --- Determine frame indices ---
        frame_indices = None
        if proc_action == "Extract frames then Coregister and Average":
            json_path = find_json(subject_id, session_id, file_name, data_dir)
            if json_path is None:
                print(f"  [SKIP] BIDS JSON not found: {tag}")
                n_skip += 1
                continue
            frame_indices = get_late_frame_indices(json_path)
            if not frame_indices:
                print(f"  [SKIP] No late-phase frames found in JSON: {tag}")
                n_skip += 1
                continue
            print(f"  [EXTRACT] {tag}: {len(frame_indices)} late frames "
                  f"(indices {frame_indices[0]}-{frame_indices[-1]})")

        # --- Load frames ---
        try:
            frames, used_indices = load_frames(nii_path, frame_indices)
        except Exception as exc:
            print(f"  [ERROR] Failed to load {nii_path}: {exc}")
            n_skip += 1
            continue

        print(f"  [COREG]  {tag}: {len(frames)} frame(s)")

        # --- Coregister and average ---
        try:
            mean_img = coregister_and_average(frames)
        except Exception as exc:
            print(f"  [ERROR] Coregistration failed: {tag}: {exc}")
            n_skip += 1
            continue

        # --- Save ---
        nib.save(mean_img, str(out_path))
        print(f"  [DONE]   {out_path}")
        n_ok += 1

    print(f"\nFinished: {n_ok} processed, {n_skip} skipped.")


if __name__ == "__main__":
    main()
