"""Per-tracer acquisition protocol definitions for OASIS-3 PET pipelines.

Each tracer has its own set of acquisition protocols (documented in the
`OASIS-3 Imaging Data Dictionary v2.3, p.22
<https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf>`_)
and a corresponding SUVR late-phase time window.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from clinica.utils.pet import Tracer


@dataclass(frozen=True)
class TracerConfig:
    """Configuration for a single PET tracer in the OASIS-3 pipeline."""

    tracer: Tracer
    bids_tracer_tag: str
    late_phase_start_min: float
    late_phase_end_min: float
    classify_fn: Callable[[dict | None, int | None, str], tuple[str, str]]

    @property
    def late_phase_start_s(self) -> float:
        return self.late_phase_start_min * 60.0

    @property
    def extract_action(self) -> str:
        start = int(self.late_phase_start_min)
        end = int(self.late_phase_end_min)
        return f"Extract {start}-{end}min frames for SUVR"

    @property
    def usable_actions(self) -> set[str]:
        return {"USE FOR SUVR", self.extract_action}

    @property
    def inventory_csv_name(self) -> str:
        return f"oasis3_{self.tracer.name.lower()}_inventory.csv"

    @property
    def no_usable_csv_name(self) -> str:
        return f"oasis3_{self.tracer.name.lower()}_no_usable_session.csv"


# ---------------------------------------------------------------------------
# AV1451 (Flortaucipir — Tau)
# ---------------------------------------------------------------------------
# Proc 1: Full dynamic 0-100 min → extract 75-100 min late frames
# Proc 2: Late static 75-100 min → directly usable
# Proc 3 early: 0-60 min → ignore


def _classify_av1451(
    frame_info: dict | None, nii_frames: int | None, modality: str
) -> tuple[str, str]:
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

    if start >= 70:
        return "Proc 2 / Proc 3 Late (75-100min)", "USE FOR SUVR"

    return f"Unclassified ({start:.0f}-{end:.0f}min, {n}f)", "Check Manually"


AV1451_CONFIG = TracerConfig(
    tracer=Tracer.AV1451,
    bids_tracer_tag=f"trc-{Tracer.AV1451.value}",
    late_phase_start_min=75.0,
    late_phase_end_min=100.0,
    classify_fn=_classify_av1451,
)


# ---------------------------------------------------------------------------
# AV45 (Florbetapir — Amyloid)
# ---------------------------------------------------------------------------
# Proc 1: 70-min dynamic from injection → extract 50-70 min late frames
# Proc 2: 20-min scan starting at 50 min → directly usable


def _classify_av45(
    frame_info: dict | None, nii_frames: int | None, modality: str
) -> tuple[str, str]:
    if modality == "CT":
        return "CT Scan", "Ignore"

    if frame_info is None:
        if nii_frames is None:
            return "Unknown (no timing data)", "Check Manually"
        if nii_frames == 1:
            return "Single Frame (pre-summed?)", "Check Manually"
        if nii_frames <= 4:
            return "Static (no timing)", "Check Manually"
        return f"Multi-frame ({nii_frames}f, no timing)", "Check Manually"

    n = frame_info["frame_count"]
    start = frame_info["first_start_min"]
    end = frame_info["last_end_min"]

    # Proc 1: full 70-min dynamic from injection
    if start < 2 and n >= 15 and end > 60:
        return "Proc 1 Full Dynamic (0-70min)", "Extract 50-70min frames for SUVR"

    # Proc 2: late scan starting at ~50 min, lasting ~20 min
    if start >= 45:
        return "Proc 2 Late Static (50-70min)", "USE FOR SUVR"

    # Early-only acquisition
    if start < 2 and end <= 45:
        return "Early Phase Only (0-45min)", "Ignore for Standard SUVR"

    return f"Unclassified ({start:.0f}-{end:.0f}min, {n}f)", "Check Manually"


AV45_CONFIG = TracerConfig(
    tracer=Tracer.AV45,
    bids_tracer_tag=f"trc-{Tracer.AV45.value}",
    late_phase_start_min=50.0,
    late_phase_end_min=70.0,
    classify_fn=_classify_av45,
)


# ---------------------------------------------------------------------------
# PIB (Pittsburgh Compound B — Amyloid)
# ---------------------------------------------------------------------------
# Proc 1: 60-min dynamic scan (24x5s + 9x20s + 10x1min + 9x5min)
# Proc 2: 30-min dynamic starting at ~30 min (6x300s) → extract 40-60 min
# SUVR window: 40-60 min


def _classify_pib(
    frame_info: dict | None, nii_frames: int | None, modality: str
) -> tuple[str, str]:
    if modality == "CT":
        return "CT Scan", "Ignore"

    if frame_info is None:
        if nii_frames is None:
            return "Unknown (no timing data)", "Check Manually"
        if nii_frames == 1:
            return "Single Frame (pre-summed?)", "Check Manually"
        if nii_frames <= 4:
            return "Static (no timing)", "Check Manually"
        return f"Multi-frame ({nii_frames}f, no timing)", "Check Manually"

    n = frame_info["frame_count"]
    start = frame_info["first_start_min"]
    end = frame_info["last_end_min"]

    # Standard 60-min dynamic scan
    if start < 2 and n >= 20 and end > 50:
        return "60min Dynamic (0-60min)", "Extract 30-60min frames for SUVR"

    # 30-min dynamic starting at ~30 min (covers 30-60 min SUVR window)
    # 30-min dynamic starting at ~30 min (covers 30-60 min SUVR window)
    if 25 <= start <= 35 and end >= 55:
        return "30min Dynamic (30-60min)", "Extract 30-60min frames for SUVR"

    # Late-phase only (if a subset was acquired separately)
    if start >= 30:
        return "Late Phase (30-60min)", "USE FOR SUVR"

    # Early-only acquisition
    if start < 2 and end <= 35:
        return "Early Phase Only", "Ignore for Standard SUVR"

    return f"Unclassified ({start:.0f}-{end:.0f}min, {n}f)", "Check Manually"


PIB_CONFIG = TracerConfig(
    tracer=Tracer.PIB,
    bids_tracer_tag=f"trc-{Tracer.PIB.value}",
    late_phase_start_min=30.0,
    late_phase_end_min=60.0,
    classify_fn=_classify_pib,
)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

TRACER_CONFIGS: dict[str, TracerConfig] = {
    "AV1451": AV1451_CONFIG,
    "AV45": AV45_CONFIG,
    "PIB": PIB_CONFIG,
}

SUPPORTED_TRACERS = list(TRACER_CONFIGS.keys())


def get_tracer_config(tracer_name: str) -> TracerConfig:
    """Look up a TracerConfig by name (case-insensitive).

    Raises
    ------
    ValueError
        If the tracer name is not supported.
    """
    key = tracer_name.upper()
    if key not in TRACER_CONFIGS:
        raise ValueError(
            f"Unsupported tracer '{tracer_name}'. "
            f"Supported tracers: {', '.join(SUPPORTED_TRACERS)}"
        )
    return TRACER_CONFIGS[key]
