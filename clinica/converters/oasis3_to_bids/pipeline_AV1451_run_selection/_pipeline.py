"""Top-level AV1451 pipeline: classify → filter → coregister/average."""

from pathlib import Path

import pandas as pd


def run_pipeline(
    data_dir: Path,
    output_dir: Path,
    inventory_dir: Path | None = None,
) -> None:
    """Run the full AV1451 tau-PET processing pipeline.

    Steps performed in order:
    1. Scan *data_dir* for AV1451 sessions and classify every PET run.
    2. Split sessions into usable / unusable; log a summary via ``cprint``.
    3. Optionally write intermediate CSVs to *inventory_dir*.
    4. Coregister and average each usable scan into a single 3-D NIfTI.

    Parameters
    ----------
    data_dir:
        Root directory containing raw OASIS-3 AV1451 session folders
        (named ``<SubjectID>_AV1451_<day>``).
    output_dir:
        Destination for ``*_AV1451_coreg_avg.nii.gz`` output files.
    inventory_dir:
        When provided, saves two CSVs:
        - ``oasis3_av1451_inventory.csv``        — full per-run inventory
        - ``oasis3_av1451_no_usable_session.csv`` — unusable sessions only
    """
    from clinica.utils.stream import cprint

    from ._pet_processing import process_all_scans
    from ._scan_classification import scan_all_sessions
    from ._session_filtering import log_session_summary, split_usable_sessions

    cprint("Scanning AV1451 session directories...", lvl="info")
    records = scan_all_sessions(data_dir)

    if not records:
        cprint("No AV1451 files found. Check --data-dir.", lvl="warning")
        return

    df = pd.DataFrame(records)

    usable_df, unusable_df = split_usable_sessions(df)
    log_session_summary(usable_df, unusable_df)

    if inventory_dir is not None:
        inventory_dir = Path(inventory_dir)
        inventory_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(inventory_dir / "oasis3_av1451_inventory.csv", index=False)
        unusable_df.to_csv(
            inventory_dir / "oasis3_av1451_no_usable_session.csv", index=False
        )
        cprint(f"Inventory CSVs saved to {inventory_dir}", lvl="info")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_ok, n_skip = process_all_scans(usable_df, data_dir, output_dir)
    cprint(f"Pipeline complete: {n_ok} processed, {n_skip} skipped.", lvl="info")
