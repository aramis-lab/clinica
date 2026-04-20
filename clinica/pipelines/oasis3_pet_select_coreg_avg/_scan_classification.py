"""PET scan classification: frame-timing parsing and protocol classification."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._tracer_config import TracerConfig

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


def parse_frame_times(metadata: dict) -> dict | None:
    """Extract timing info from the FrameTimes JSON structure.

    Returns a dict with keys:
        frame_count, first_start_min, last_end_min, duration_min, durations_summary
    or None when the structure is absent or malformed.
    """
    values = None
    for parent_key in ("FrameTimes", "Time"):
        ft = metadata.get(parent_key)
        if ft and isinstance(ft, dict):
            inner = ft.get("FrameTimes")
            if inner and isinstance(inner, dict):
                values = inner.get("Values")
                if values:
                    break
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


def evaluate_bids_pet_dir(
    pet_dir: Path, subject_id: str, session_id: str, tracer_cfg: TracerConfig
) -> list[dict]:
    """Evaluate all PET runs for a given tracer in a single BIDS ``pet/`` directory.

    Parameters
    ----------
    pet_dir:
        A BIDS ``pet/`` directory, e.g. ``bids_dir/sub-OAS30001/ses-M126/pet``.
    subject_id:
        BIDS subject label, e.g. ``sub-OAS30001``.
    session_id:
        BIDS session label, e.g. ``ses-M126``.
    tracer_cfg:
        Tracer-specific configuration.

    Returns
    -------
    List of scan record dicts, one per matching JSON sidecar found.
    """
    from ._nii_reading import get_nii_frame_count

    records = []
    for json_path in sorted(
        pet_dir.glob(f"*{tracer_cfg.bids_tracer_tag}*_pet.json")
    ):
        with open(json_path, "r", errors="ignore") as f:
            try:
                metadata = json.load(f)
            except json.JSONDecodeError:
                continue

        modality = metadata.get("Modality", "PET")
        frame_info = parse_frame_times(metadata)
        nii_frames = None if frame_info is not None else get_nii_frame_count(json_path)

        filename = json_path.name
        run_parts = [p for p in filename.split("_") if p.startswith("run-")]
        run_id = run_parts[0] if run_parts else "Unknown"

        likely_proc, action = tracer_cfg.classify_fn(frame_info, nii_frames, modality)
        inj_time = metadata.get("InjectionStart") or metadata.get("TimeZero", "Unknown")

        records.append(
            {
                "Subject_ID": subject_id,
                "Session_ID": session_id,
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


def scan_all_sessions(bids_dir: Path, tracer_cfg: TracerConfig) -> list[dict]:
    """Scan *bids_dir* for all PET sessions of the given tracer and evaluate each one.

    Expects a standard BIDS layout: ``bids_dir/sub-*/ses-*/pet/``.

    Returns a flat list of scan record dicts across all subjects and sessions.
    """
    records = []
    for pet_dir in sorted(bids_dir.glob("sub-*/ses-*/pet")):
        if not pet_dir.is_dir():
            continue
        subject_id = pet_dir.parts[-3]
        session_id = pet_dir.parts[-2]
        records.extend(
            evaluate_bids_pet_dir(pet_dir, subject_id, session_id, tracer_cfg)
        )
    return records
