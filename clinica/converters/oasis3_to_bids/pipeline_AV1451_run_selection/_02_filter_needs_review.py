"""
Step 2: Identify sessions that need manual review.

Reads the inventory from step 1, finds sessions without a directly usable
or extractable scan, and generates a review template CSV.

Usage:  python 02_filter_needs_review.py

Input:  inventory/oasis3_av1451_inventory.csv
Output: inventory/oasis3_av1451_needs_review.csv
        inventory/oasis3_av1451_strict_needs_review.csv
        inventory/oasis3_av1451_review_decisions.csv  (template for manual review)

*** After this step, copy your reviewed decisions file to:
    inventory/oasis3_av1451_reviewed_decisions.csv
    Then run step 3. ***
"""

import os
from pathlib import Path

import pandas as pd

# --- DEFAULT CONFIGURATION ---
DEFAULT_INVENTORY_DIR = Path(r"D:\Oasis3\inventory")

READY_ACTION       = "USE FOR SUVR"
EXTRACTABLE_ACTION = "Extract 75-100min frames for SUVR"

SESSION_KEY        = ["Subject_ID", "Session_ID"]
REVIEW_KEY_COLS    = ["Subject_ID", "Session_ID", "File_Name"]
REVIEW_EXTRA_COLS  = ["Reviewed_Action", "Reviewer_Notes", "Review_Date"]


def main(inventory_dir=None):
    inventory_dir = Path(inventory_dir) if inventory_dir else DEFAULT_INVENTORY_DIR

    input_csv            = inventory_dir / "oasis3_av1451_inventory.csv"
    output_csv           = inventory_dir / "oasis3_av1451_needs_review.csv"
    strict_output_csv    = inventory_dir / "oasis3_av1451_strict_needs_review.csv"
    review_decisions_csv = inventory_dir / "oasis3_av1451_review_decisions.csv"

    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} runs from inventory.")

    # --- Sessions with a ready-to-use run ---
    has_suvr = df.groupby(SESSION_KEY)["Action"].transform(
        lambda actions: (actions == READY_ACTION).any()
    )
    needs_review = df[~has_suvr].copy()
    needs_review.to_csv(output_csv, index=False)

    total_sessions  = df.groupby(SESSION_KEY).ngroups
    ready_sessions  = df[has_suvr].groupby(SESSION_KEY).ngroups
    review_sessions = needs_review.groupby(SESSION_KEY).ngroups

    print(f"\nTotal sessions: {total_sessions}")
    print(f"Sessions with ready SUVR scan: {ready_sessions}")
    print(f"Sessions needing review: {review_sessions}")
    print(f"Saved: {output_csv}")

    # --- Strict filter: also exclude sessions with extractable dynamic scans ---
    usable_actions = {READY_ACTION, EXTRACTABLE_ACTION}
    has_any_usable = df.groupby(SESSION_KEY)["Action"].transform(
        lambda actions: actions.isin(usable_actions).any()
    )
    strict_needs_review = df[~has_any_usable].copy()
    strict_needs_review.to_csv(strict_output_csv, index=False)

    strict_sessions = strict_needs_review.groupby(SESSION_KEY).ngroups
    print(f"\nStrict filter (also excludes extractable scans):")
    print(f"  Sessions needing review: {strict_sessions}")
    print(f"Saved: {strict_output_csv}")

    # --- Create/update review template ---
    init_review_file(strict_needs_review, review_decisions_csv)

    print(f"\n{'='*60}")
    print(f"Next: copy your reviewed decisions file to:")
    print(f"  {inventory_dir / 'oasis3_av1451_reviewed_decisions.csv'}")
    print(f"Then run:  python 03_build_master_inventory.py")
    print(f"{'='*60}")


def init_review_file(strict_df, review_decisions_csv):
    """Create or update the review decisions template."""
    new_entries = strict_df[REVIEW_KEY_COLS].copy()
    for col in REVIEW_EXTRA_COLS:
        new_entries[col] = ""

    if not review_decisions_csv.exists():
        new_entries.to_csv(review_decisions_csv, index=False)
        print(f"\nCreated review template: {review_decisions_csv}")
        print(f"  {len(new_entries)} entries to review.")
        return

    existing = pd.read_csv(review_decisions_csv, dtype=str).fillna("")
    merged = existing.merge(
        new_entries[REVIEW_KEY_COLS], on=REVIEW_KEY_COLS,
        how="outer", indicator=True,
    )
    new_rows = merged[merged["_merge"] == "right_only"].drop(columns="_merge")
    for col in REVIEW_EXTRA_COLS:
        if col not in new_rows.columns:
            new_rows[col] = ""

    updated = pd.concat([existing, new_rows[existing.columns]], ignore_index=True)
    updated.to_csv(review_decisions_csv, index=False)

    reviewed   = (updated["Reviewed_Action"] != "").sum()
    unreviewed = (updated["Reviewed_Action"] == "").sum()
    print(f"\nUpdated review template: {review_decisions_csv}")
    print(f"  Reviewed: {reviewed}  |  Unreviewed: {unreviewed}")
    if len(new_rows) > 0:
        print(f"  New entries added: {len(new_rows)}")


if __name__ == "__main__":
    main()
