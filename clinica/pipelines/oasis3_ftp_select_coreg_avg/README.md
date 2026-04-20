# AV1451 Tau-PET Run Selection Pipeline — Process Report

## Overview

The AV1451 Run Selection Pipeline systematically evaluates, classifies, and (after 
a manual revision step) process AV1451 (flortaucipir) tau-PET scans from the OASIS-3. 
The core challenge it addresses is that OASIS-3 tau-PET sessions frequently contain 
multiple runs acquired under different protocols (documented in p.22 in the [OASIS-3 
Imaging Methods and Data Dictionary](https://bpb-us-e2.wpmucdn.com/sites.wustl.edu/dist/6/4383/files/2024/04/OASIS-3_Imaging_Data_Dictionary_v2.3-a93c947a586e7367.pdf)), 
where only the late-phase (75–100 min post-injection) data is suitable for computing 
Standardized Uptake Value Ratios (SUVR). 
The pipeline ensures that, independently of the recording protocol, only the 
late-phase images are selected for processing.

## Step 1 — Evaluate Tau Scans

The pipeline begins by scanning the raw data directory for all session folders containing `AV1451` in their name. For each session, it locates BIDS JSON sidecar files and parses the `FrameTimes` metadata structure to extract timing information: frame count, first frame start time, last frame end time, total scan duration, and a human-readable summary of frame durations. When timing metadata is unavailable, it falls back to reading the NIfTI header via `nibabel` to determine the number of frames.

Each run is then classified against known OASIS-3 AV1451 acquisition protocols. Runs starting near time zero with more than 30 frames spanning over 70 minutes are identified as **Procedure 3** (full dynamic 0–100 min), from which late-phase frames can be extracted. Runs starting at or after 70 minutes are classified as **Procedure 2** (late static, directly usable for SUVR). Early acquisitions (0–25 min) are labeled **Procedure 1** and marked to be ignored. CT scans are also excluded. The output is a comprehensive inventory CSV with one row per run.

## Step 2 — Filter and Flag for Manual Review

The second step reads the inventory and groups runs by session. Sessions that already have a directly usable Procedure 2 scan are marked as ready. The remaining sessions are flagged for review under two filtering levels: a standard filter (no direct SUVR scan) and a strict filter (no direct or extractable scan). A review decisions template CSV is generated, containing columns for the reviewer to specify a `Reviewed_Action`, notes, and date.

Only ~50 sessions require manual revision. The manual revision has already been 
performed. However, sharing of the manually filled ``*_reviewed.csv`` table is 
currently kept private to comply with the Data Usage Agreement. There should be a 
discussion on how is this managed by the clinica pipeline.

## Step 3 — Build Master Inventory

After manual review, the third step merges the original automated inventory with the reviewed decisions file (semicolon-delimited). It computes session-level flags — whether any run is directly usable, extractable, or was approved during manual review — and derives a `Final_Status` and `Processing_Action` for every run. Scans are categorized as Included, Excluded (with reason), Superseded, or Pending Review. The output is a single master CSV that serves as the definitive record for downstream processing.

## Step 4 — Coregister and Average

The final step processes all included scans. For each scan marked as "Coregister and Average", it loads the 4-D NIfTI, splits it into individual 3-D frames, and performs rigid-body coregistration to the first frame using ANTs (`antspyx`). For full-dynamic scans requiring frame extraction, it first selects only frames acquired at or after 75 minutes post-injection. All aligned frames are then averaged into a single 3-D volume and saved as a compressed NIfTI file. Already-processed scans are automatically skipped on re-runs, making the pipeline idempotent.

## Summary

This pipeline transforms an unstructured collection of multi-protocol tau-PET acquisitions into a curated, analysis-ready dataset by combining automated protocol classification, structured manual review, and robust image processing.
