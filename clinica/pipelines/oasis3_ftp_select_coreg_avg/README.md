# AV1451 Tau-PET Run Selection Pipeline

## Overview

The AV1451 Run Selection Pipeline evaluates, classifies, and processes AV1451
(flortaucipir) tau-PET scans from the OASIS-3 dataset in BIDS format.
The core challenge it addresses is that OASIS-3 tau-PET sessions frequently
contain multiple runs acquired under different protocols (documented in p.22 of
the [OASIS-3 Imaging Data Dictionary](https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf)),
where only the late-phase (75-100 min post-injection) data is suitable for
computing Standardized Uptake Value Ratios (SUVR).
The pipeline ensures that, independently of the recording protocol, only the
late-phase images are selected for processing.

## Pipeline Steps

### Step 1 — Scan Classification (`_scan_classification.py`)

The pipeline scans the BIDS dataset (`sub-*/ses-*/pet/`) for all JSON sidecars
matching the `trc-18FAV1451` tracer tag. For each run it parses the `FrameTimes`
metadata to extract timing information: frame count, first frame start time,
last frame end time, total scan duration, and a human-readable summary of frame
durations. When timing metadata is unavailable, it falls back to reading the
NIfTI header via `nibabel` to determine the number of frames.

Each run is classified against the known OASIS-3 AV1451 acquisition protocols:

| Classification | Criteria | Action |
|---|---|---|
| **Proc 1 Full (0-100 min)** | Start near 0, >30 frames, end >70 min | Extract 75-100 min frames for SUVR |
| **Proc 2 / Proc 3 Late (75-100 min)** | Start >= 70 min | USE FOR SUVR |
| **Proc 3 Early (0-60 min)** | Start near 0, end <= 62 min | Ignore for Standard SUVR |
| **CT Scan** | Modality = CT | Ignore |

### Step 2 — Session Filtering (`_session_filtering.py`)

Sessions are split into usable and unusable groups. A session is **usable** when
at least one of its runs has an action in `{"USE FOR SUVR", "Extract 75-100min
frames for SUVR"}`. Sessions where no run qualifies are flagged as unusable and
logged as warnings.

When an `inventory_dir` is provided, two CSVs are saved:

- `oasis3_av1451_inventory.csv` — full per-run inventory with classification
- `oasis3_av1451_no_usable_session.csv` — unusable sessions only

### Step 3 — Coregister and Average (`_pet_processing.py`)

All usable scans are processed:

1. For full-dynamic scans (Proc 1), only frames starting at >= 75 min
   post-injection are extracted.
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
from clinica.pipelines.oasis3_ftp_select_coreg_avg import run_pipeline

run_pipeline(
    bids_dir="path/to/bids",
    output_dir="path/to/output",
    inventory_dir="path/to/inventory",  # optional
)
```

### CLI

```bash
python -m clinica.pipelines.oasis3_ftp_select_coreg_avg BIDS_DIR OUTPUT_DIR [--inventory-dir DIR]
```

## Module Structure

| Module | Responsibility |
|---|---|
| `_pipeline.py` | Top-level orchestration |
| `_scan_classification.py` | Frame-timing parsing and protocol classification |
| `_session_filtering.py` | Session-level usability filtering and logging |
| `_pet_processing.py` | Late-frame extraction, coregistration, and averaging |
| `_nii_reading.py` | NIfTI and BIDS JSON I/O helpers |
| `cli.py` | Click CLI entry point |

## Dependencies

- `pandas` — tabular data handling
- `nibabel` — NIfTI I/O
- `antspyx` — rigid-body coregistration
