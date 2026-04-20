# OASIS-3 PET Run Selection Pipeline

## Overview

The OASIS-3 PET Run Selection Pipeline evaluates, classifies, and processes
PET scans from the OASIS-3 dataset in BIDS format. It supports multiple
tracers — **AV1451** (flortaucipir, tau), **AV45** (florbetapir, amyloid),
and **PIB** (Pittsburgh Compound B, amyloid) — each with its own acquisition
protocols and SUVR time windows.

The core challenge it addresses is that OASIS-3 PET sessions frequently
contain multiple runs acquired under different protocols (documented in p.22
of the [OASIS-3 Imaging Data Dictionary](https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf)),
where only the late-phase data is suitable for computing Standardized Uptake
Value Ratios (SUVR). The pipeline ensures that, independently of the recording
protocol, only the appropriate late-phase images are selected for processing.

## Supported Tracers

| Tracer | BIDS Tag | SUVR Window | Description |
|---|---|---|---|
| **AV1451** | `trc-18FAV1451` | 75-100 min | Flortaucipir (tau) |
| **AV45** | `trc-18FAV45` | 50-70 min | Florbetapir (amyloid) |
| **PIB** | `trc-11CPIB` | 40-60 min | Pittsburgh Compound B (amyloid) |

### AV1451 Acquisition Protocols

| Classification | Criteria | Action |
|---|---|---|
| **Proc 1 Full (0-100 min)** | Start near 0, >30 frames, end >70 min | Extract 75-100 min frames for SUVR |
| **Proc 2 / Proc 3 Late (75-100 min)** | Start >= 70 min | USE FOR SUVR |
| **Proc 3 Early (0-60 min)** | Start near 0, end <= 62 min | Ignore for Standard SUVR |

### AV45 Acquisition Protocols

| Classification | Criteria | Action |
|---|---|---|
| **Proc 1 Full Dynamic (0-70 min)** | Start near 0, >=15 frames, end >60 min | Extract 50-70 min frames for SUVR |
| **Proc 2 Late Static (50-70 min)** | Start >= 45 min | USE FOR SUVR |
| **Early Phase Only (0-45 min)** | Start near 0, end <= 45 min | Ignore for Standard SUVR |

### PIB Acquisition Protocols

| Classification | Criteria | Action |
|---|---|---|
| **60 min Dynamic (0-60 min)** | Start near 0, >=20 frames, end >50 min | Extract 40-60 min frames for SUVR |
| **Late Phase (40-60 min)** | Start >= 35 min | USE FOR SUVR |
| **Early Phase Only** | Start near 0, end <= 35 min | Ignore for Standard SUVR |

## Pipeline Steps

### Step 1 — Scan Classification (`_scan_classification.py`)

The pipeline scans the BIDS dataset (`sub-*/ses-*/pet/`) for all JSON sidecars
matching the selected tracer's BIDS tag. For each run it parses the `FrameTimes`
metadata to extract timing information: frame count, first frame start time,
last frame end time, total scan duration, and a human-readable summary of frame
durations. When timing metadata is unavailable, it falls back to reading the
NIfTI header via `nibabel` to determine the number of frames.

Each run is then classified using tracer-specific rules defined in
`_tracer_config.py`.

### Step 2 — Session Filtering (`_session_filtering.py`)

Sessions are split into usable and unusable groups. A session is **usable** when
at least one of its runs has a usable action for the selected tracer. Sessions
where no run qualifies are flagged as unusable and logged as warnings.

When an `inventory_dir` is provided, two CSVs are saved (named after the tracer):

- `oasis3_<tracer>_inventory.csv` — full per-run inventory with classification
- `oasis3_<tracer>_no_usable_session.csv` — unusable sessions only

### Step 3 — Coregister and Average (`_pet_processing.py`)

All usable scans are processed:

1. For full-dynamic scans, only frames starting at or after the tracer's
   late-phase threshold are extracted.
2. The selected frames are split into individual 3-D volumes.
3. Rigid-body (affine) coregistration to the first frame is performed using
   ANTs (`antspyx`).
4. All aligned frames are averaged into a single 3-D volume and saved as a
   compressed NIfTI (`*_coreg_avg.nii.gz`).

Already-processed scans are automatically skipped on re-runs, making the
pipeline idempotent.

## Usage

### Python API

```python
from clinica.pipelines.oasis3_pet_select_coreg_avg._pipeline import run_pipeline

run_pipeline(
    bids_dir="path/to/bids",
    output_dir="path/to/output",
    tracer="AV1451",              # or "AV45" or "PIB"
    inventory_dir="path/to/inventory",  # optional
)
```

### CLI

```bash
python -m clinica.pipelines.oasis3_pet_select_coreg_avg BIDS_DIR OUTPUT_DIR --tracer AV1451 [--inventory-dir DIR]
```

## Module Structure

| Module | Responsibility |
|---|---|
| `_pipeline.py` | Top-level orchestration |
| `_tracer_config.py` | Per-tracer protocol definitions and SUVR windows |
| `_scan_classification.py` | Frame-timing parsing and protocol classification |
| `_session_filtering.py` | Session-level usability filtering and logging |
| `_pet_processing.py` | Late-frame extraction, coregistration, and averaging |
| `_nii_reading.py` | NIfTI and BIDS JSON I/O helpers |
| `cli.py` | Click CLI entry point |

## Dependencies

- `pandas` — tabular data handling
- `nibabel` — NIfTI I/O
- `antspyx` — rigid-body coregistration
