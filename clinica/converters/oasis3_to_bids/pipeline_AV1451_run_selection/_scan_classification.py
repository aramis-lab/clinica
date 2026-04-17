"""AV1451 scan classification: frame-timing parsing and protocol classification."""

import json
from pathlib import Path

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

USABLE_ACTIONS = {"USE FOR SUVR", "Extract 75-100min frames for SUVR"}


def parse_frame_times(metadata: dict) -> dict | None:
    """Extract timing info from the FrameTimes JSON structure.

    Returns a dict with keys:
        frame_count, first_start_min, last_end_min, duration_min, durations_summary
    or None when the structure is absent or malformed.
    """
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


def classify_pet_run(
    frame_info: dict | None, nii_frames: int | None, modality: str
) -> tuple[str, str]:
    """Classify a PET run against the known OASIS-3 AV1451 acquisition protocols.

    Protocols sourced from the OASIS-3 Imaging Data Dictionary v2.3 (p.22).
    https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/
    OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf

    Returns
    -------
    (likely_procedure, action)
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
            return "Proc 2 / Proc 3 Late Static (75-100min)", "USE FOR SUVR"
        return f"Multi-frame ({nii_frames}f, no timing)", "Check Manually"

    n = frame_info["frame_count"]
    start = frame_info["first_start_min"]
    end = frame_info["last_end_min"]

    if start < 2 and n >= 30 and end > 70:
        return "Proc 1 Full (0-100min)", "Extract 75-100min frames for SUVR"

    if start < 2 and end <= 62:
        if frame_info["durations_summary"] == "6x300s" and n == 6:
            return "Proc 2 / Proc Late Static (75-100min)", "USE FOR SUVR"
        return "Proc 3 Early (0-60min)", "Ignore for Standard SUVR"

    if start >= 70 and n == 6:
        return "Proc 2 / Proc 3 Late (75-100min)", "USE FOR SUVR"

    return f"Unclassified ({start:.0f}-{end:.0f}min, {n}f)", "Check Manually"


def evaluate_session_dir(session_dir: Path) -> list[dict]:
    """Evaluate all PET runs in a single AV1451 session directory.

    Parameters
    ----------
    session_dir:
        A directory named like ``OAS30001_AV1451_d0761``.

    Returns
    -------
    List of scan record dicts, one per BIDS JSON sidecar found.
    """
    from ._nii_reading import get_nii_frame_count

    parts = session_dir.name.split("_")
    subject_id = parts[0] if parts else "Unknown"
    session_day = parts[2] if len(parts) >= 3 else "Unknown"

    records = []
    for pet_dir in sorted(
        d for d in session_dir.iterdir() if d.is_dir() and d.name.startswith("pet")
    ):
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
                    continue

            mod_str = metadata.get("Modality", "")
            modality = "CT" if "CT" in mod_str.upper() else "PET"

            frame_info = parse_frame_times(metadata)
            nii_frames = (
                None if frame_info is not None else get_nii_frame_count(json_path)
            )

            filename = json_path.name
            run_parts = [p for p in filename.split("_") if p.startswith("run-")]
            run_id = run_parts[0] if run_parts else "Unknown"

            likely_proc, action = classify_pet_run(frame_info, nii_frames, modality)
            inj_time = metadata.get("InjectionStart") or metadata.get(
                "TimeZero", "Unknown"
            )

            records.append(
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
    return records


def scan_all_sessions(data_dir: Path) -> list[dict]:
    """Scan *data_dir* for all AV1451 session directories and evaluate each one.

    Returns a flat list of scan record dicts across all sessions.
    """
    av1451_dirs = sorted(
        d for d in data_dir.iterdir() if d.is_dir() and "_AV1451_" in d.name
    )
    records = []
    for session_dir in av1451_dirs:
        records.extend(evaluate_session_dir(session_dir))
    return records
