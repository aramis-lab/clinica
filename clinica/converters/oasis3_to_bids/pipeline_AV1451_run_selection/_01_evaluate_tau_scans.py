"""
Step 1: Evaluate all AV1451 tau-PET scans in raw data.

Scans the raw data directory, parses BIDS JSON metadata and NIfTI headers,
classifies each run by acquisition protocol, and writes the inventory CSV.

Usage:  python 01_evaluate_tau_scans.py

Output: inventory/oasis3_av1451_inventory.csv
"""

import csv
import json
from collections import Counter
from pathlib import Path

try:
    import nibabel as nib

    HAS_NIBABEL = True
except ImportError:
    HAS_NIBABEL = False

# --- DEFAULT CONFIGURATION ---
DEFAULT_DATA_DIR = Path(r"E:\Oasis3\raw")
DEFAULT_OUTPUT_DIR = Path(r"E:\Oasis3\inventory")

CSV_FIELDS = [
    "Subject_ID",
    "Session_ID",
    "File_Name",
    "Run_Number",
    "Modality",
    "Frame_Count",
    "First_Frame_Start_Min",
    "Last_Frame_End_Min",
    "Total_Scan_Duration_Min",
    "Frame_Durations_Summary",
    "Injection_Time",
    "Likely_Procedure",
    "Action",
]


def parse_frame_times(metadata):
    """Extract timing info from FrameTimes JSON structure."""
    ft = metadata.get("FrameTimes")
    if not ft or not isinstance(ft, dict):
        return None

    inner = ft.get("FrameTimes")
    if not inner or not isinstance(inner, dict):
        return None

    values = inner.get("Values")
    if not values or not isinstance(values, list) or len(values) == 0:
        return None

    frame_count = len(values)
    first_start_s = values[0][0]
    last_start_s = values[-1][0]
    last_duration_s = values[-1][2]
    last_end_s = last_start_s + last_duration_s

    durations = [v[2] for v in values]
    dur_counts = []
    current_dur = durations[0]
    current_count = 1
    for d in durations[1:]:
        if d == current_dur:
            current_count += 1
        else:
            dur_counts.append(f"{current_count}x{current_dur:.0f}s")
            current_dur = d
            current_count = 1
    dur_counts.append(f"{current_count}x{current_dur:.0f}s")

    return {
        "frame_count": frame_count,
        "first_start_min": round(first_start_s / 60.0, 1),
        "last_end_min": round(last_end_s / 60.0, 1),
        "duration_min": round((last_end_s - first_start_s) / 60.0, 1),
        "durations_summary": ", ".join(dur_counts),
    }


def classify_pet_run(frame_info, nii_frames, modality):
    """
    Classify the run based on OASIS3 AV1451 protocols.

    Protocols extracted from the OASIS3 Imaging Data Dictionary, accessed 17/04/2026
    https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf
    """
    if modality == "CT":
        return "CT Scan", "Ignore"

    if frame_info is None:
        if nii_frames is None:
            return "Unknown (no timing data)", "Check Manually"
        if nii_frames == 1:
            return "Single Frame (pre-summed?)", "Check Manually"
        if nii_frames < 6:
            return "Static (no timing)", "Check Manually"
        if nii_frames == 6:
            # Several cases of late recordings with 6 frames of 5min duration
            return "Proc 2 / Proc 3 Late Static (75-100min)", "USE FOR SUVR"

        return f"Multi-frame ({nii_frames}f, no timing)", "Check Manually"

    n = frame_info["frame_count"]
    start = frame_info["first_start_min"]
    end = frame_info["last_end_min"]

    if start < 2 and n >= 30 and end > 70:
        return "Proc 1 Full (0-100min)", "Extract 75-100min frames for SUVR"

    if start < 2 and end <= 62:
        # Some cases, the run contains 6x5mins frames but start at 0 and end in 30mins
        if frame_info["durations_summary"] == "6x300s" and n == 6:
            return "Proc 2 / Proc Late Static (75-100min)", "USE FOR SUVR"
        return "Proc 3 Early (0-60min)", "Ignore for Standard SUVR"

    if start >= 70:
        return "Proc 2 / Proc 3 Late (75-100min)", "USE FOR SUVR"

    return f"Unclassified ({start:.0f}-{end:.0f}min, {n}f)", "Check Manually"


def get_nii_frame_count(json_path):
    """Load NIfTI to get frame count as fallback."""
    if not HAS_NIBABEL:
        return None

    pet_dir = json_path.parent.parent
    nii_dir = pet_dir / "NIFTI"
    for ext in (".nii.gz", ".nii"):
        nii_path = nii_dir / (json_path.stem + ext)
        if nii_path.exists():
            try:
                img = nib.load(str(nii_path))
                shape = img.shape
                return shape[3] if len(shape) > 3 else 1
            except Exception:
                return None
    return None


def main(data_dir=None, output_dir=None):
    data_dir = Path(data_dir) if data_dir else DEFAULT_DATA_DIR
    output_dir = Path(output_dir) if output_dir else DEFAULT_OUTPUT_DIR
    output_csv = output_dir / "oasis3_av1451_inventory.csv"

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning {data_dir} for AV1451 sessions...")

    av1451_dirs = sorted(
        d for d in data_dir.iterdir() if d.is_dir() and "_AV1451_" in d.name
    )
    print(f"Found {len(av1451_dirs)} AV1451 session directories.")

    results = []
    for session_dir in av1451_dirs:
        parts = session_dir.name.split("_")
        subject_id = parts[0] if parts else "Unknown"
        session_day = parts[2] if len(parts) >= 3 else "Unknown"

        pet_dirs = sorted(
            d for d in session_dir.iterdir() if d.is_dir() and d.name.startswith("pet")
        )

        for pet_dir in pet_dirs:
            bids_dir = pet_dir / "BIDS"
            if not bids_dir.exists():
                continue

            for json_path in bids_dir.glob("*.json"):
                if "dataset_description" in json_path.name:
                    continue

                with open(json_path, "r", errors="ignore") as f:
                    try:
                        metadata = json.load(f)
                    except json.JSONDecodeError:
                        print(f"  Warning: Could not parse {json_path}")
                        continue

                mod_str = metadata.get("Modality", "")
                modality = "CT" if "CT" in mod_str.upper() else "PET"

                frame_info = parse_frame_times(metadata)

                nii_frames = None
                if frame_info is None:
                    nii_frames = get_nii_frame_count(json_path)

                filename = json_path.name
                run_parts = [p for p in filename.split("_") if p.startswith("run-")]
                run_id = run_parts[0] if run_parts else "Unknown"

                likely_proc, action = classify_pet_run(frame_info, nii_frames, modality)
                inj_time = metadata.get("InjectionStart") or metadata.get(
                    "TimeZero", "Unknown"
                )

                results.append(
                    {
                        "Subject_ID": subject_id,
                        "Session_ID": session_day,
                        "File_Name": filename,
                        "Run_Number": run_id,
                        "Modality": modality,
                        "Frame_Count": frame_info["frame_count"]
                        if frame_info
                        else (nii_frames or "Unknown"),
                        "First_Frame_Start_Min": frame_info["first_start_min"]
                        if frame_info
                        else "N/A",
                        "Last_Frame_End_Min": frame_info["last_end_min"]
                        if frame_info
                        else "N/A",
                        "Total_Scan_Duration_Min": frame_info["duration_min"]
                        if frame_info
                        else "N/A",
                        "Frame_Durations_Summary": frame_info["durations_summary"]
                        if frame_info
                        else "N/A",
                        "Injection_Time": inj_time,
                        "Likely_Procedure": likely_proc,
                        "Action": action,
                    }
                )

    if not results:
        print("No AV1451 files found. Check your directory path.")
        return
    import pandas as pd

    df = pd.DataFrame(results)
    print(df["Likely_Procedure"].unique())
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(results)

    action_counts = Counter(r["Action"] for r in results)
    print(f"\nDone! Evaluated {len(results)} AV1451 files.")
    print(f"Output: {output_csv}")
    print("Summary:")
    for action, count in action_counts.most_common():
        print(f"  {action}: {count}")


if __name__ == "__main__":
    main()
