"""Top-level PET pipeline: classify → filter → coregister/average."""

from pathlib import Path

import pandas as pd


def run_pipeline(
    bids_dir: Path,
    output_dir: Path,
    tracer: str = "AV1451",
    inventory_dir: Path | None = None,
    bids_output: bool = False,
) -> None:
    """Run the full PET processing pipeline on a BIDS dataset.

    Steps performed in order:
    1. Scan *bids_dir* for sessions of the given tracer and classify every PET run.
    2. Split sessions into usable / unusable; log a summary via ``cprint``.
    3. Optionally write intermediate CSVs to *inventory_dir*.
    4. Coregister and average each usable scan into a single 3-D NIfTI.

    Parameters
    ----------
    bids_dir:
        BIDS dataset root directory containing ``sub-*/ses-*/pet/`` folders.
    output_dir:
        Destination for ``*_coreg_avg.nii.gz`` output files.
    tracer:
        Tracer name (case-insensitive). Supported: ``AV1451``, ``AV45``, ``PIB``.
    inventory_dir:
        When provided, saves two CSVs (named after the tracer):
        - ``oasis3_<tracer>_inventory.csv``         — full per-run inventory
        - ``oasis3_<tracer>_no_usable_session.csv`` — unusable sessions only
    bids_output:
        When True, averaged images are saved back into the BIDS session ``pet/``
        directory with a BIDS-compliant name (no ``run`` entity, no ``coreg_avg``
        suffix).  *output_dir* is ignored for NIfTI outputs in this mode.
    """
    from clinica.utils.stream import cprint

    from ._pet_processing import process_all_scans
    from ._scan_classification import scan_all_sessions
    from ._session_filtering import log_session_summary, split_usable_sessions
    from ._tracer_config import get_tracer_config

    tracer_cfg = get_tracer_config(tracer)
    tracer_name = tracer_cfg.tracer.name

    cprint(f"Scanning {tracer_name} sessions in BIDS dataset...", lvl="info")
    records = scan_all_sessions(bids_dir, tracer_cfg)

    if not records:
        cprint(
            f"No {tracer_name} files found. Check the bids_dir argument.",
            lvl="warning",
        )
        return

    df_records = pd.DataFrame(records)

    usable_df, unusable_df = split_usable_sessions(df_records, tracer_cfg)
    log_session_summary(usable_df, unusable_df, tracer_name)

    if inventory_dir is not None:
        inventory_dir = Path(inventory_dir)
        inventory_dir.mkdir(parents=True, exist_ok=True)
        df_records.to_csv(inventory_dir / tracer_cfg.inventory_csv_name, index=False)
        unusable_df.to_csv(inventory_dir / tracer_cfg.no_usable_csv_name, index=False)
        cprint(f"Inventory CSVs saved to {inventory_dir}", lvl="info")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_ok, n_skip = process_all_scans(
        usable_df, bids_dir, output_dir, tracer_cfg, bids_output=bids_output
    )
    cprint(f"Pipeline complete: {n_ok} processed, {n_skip} skipped.", lvl="info")
