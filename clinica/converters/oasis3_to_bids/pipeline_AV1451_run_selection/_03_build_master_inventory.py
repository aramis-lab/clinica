"""
Step 3: Build the master inventory.

Merges the automated inventory (step 1) with the manually reviewed
decisions, determines Final_Status and Processing_Action for every scan.

Usage:  python 03_build_master_inventory.py

Input:  inventory/oasis3_av1451_inventory.csv
        inventory/oasis3_av1451_reviewed_decisions.csv  (semicolon-delimited)
Output: inventory/oasis3_av1451_master.csv
"""

from pathlib import Path

import pandas as pd

# --- DEFAULT CONFIGURATION ---
DEFAULT_INVENTORY_DIR = Path(r"D:\Oasis3\inventory")

SESSION_KEY = ["Subject_ID", "Session_ID"]

ACTION_SUVR    = "USE FOR SUVR"
ACTION_EXTRACT = "Extract 75-100min frames for SUVR"
ACTION_IGNORE  = "Ignore for Standard SUVR"
ACTION_CT      = "CT Scan"


def derive_status_action(row, has_direct_suvr, has_extractable, has_reviewed_suvr):
    auto     = str(row["Automated_Action"]).strip()
    reviewed_raw = row.get("Reviewed_Action")
    reviewed = "" if pd.isna(reviewed_raw) else str(reviewed_raw).strip()
    frames   = row.get("Frame_Count")

    # --- Manual review overrides automated classification ---
    if reviewed:
        if reviewed == ACTION_SUVR:
            if str(frames) == "1":
                return "Included (manual review)", "No action needed (already averaged)"
            return "Included (manual review)", "Coregister and Average"
        if reviewed == ACTION_IGNORE:
            return "Excluded (manual review)", ""
        return "Excluded (manual review)", ""

    # --- CT: always excluded ---
    if auto == ACTION_CT:
        return "Excluded (CT)", ""

    # --- Directly usable (automated) ---
    if auto == ACTION_SUVR:
        if str(frames) == "1":
            return "Included", "No action needed (already averaged)"
        return "Included", "Coregister and Average"

    # --- Full dynamic needing frame extraction ---
    if auto == ACTION_EXTRACT:
        if has_direct_suvr:
            return "Superseded (session has Proc 2 direct)", ""
        return "Included", "Extract frames then Coregister and Average"

    # --- Not useful for SUVR (Proc 1 early, etc.) ---
    if auto == ACTION_IGNORE:
        return "Excluded (not late-phase)", ""

    # --- Unreviewed scans that needed manual review ---
    session_has_usable = has_direct_suvr or has_extractable or has_reviewed_suvr
    if session_has_usable:
        return "Discarded (session has better scan)", ""

    return "Pending Review", ""


def main(inventory_dir=None):
    inventory_dir = Path(inventory_dir) if inventory_dir else DEFAULT_INVENTORY_DIR

    inventory_csv          = inventory_dir / "oasis3_av1451_inventory.csv"
    reviewed_decisions_csv = inventory_dir / "oasis3_av1451_reviewed_decisions.csv"
    output_csv             = inventory_dir / "oasis3_av1451_master.csv"

    inventory = pd.read_csv(inventory_csv)
    print(f"Loaded inventory: {len(inventory)} rows")

    if not reviewed_decisions_csv.exists():
        print(f"\nWARNING: Reviewed decisions file not found at:")
        print(f"  {reviewed_decisions_csv}")
        print(f"Proceeding without manual reviews (only automated classifications).\n")
        reviewed = pd.DataFrame(columns=["Subject_ID", "Session_ID", "File_Name",
                                         "Reviewed_Action", "Reviewer_Notes", "Review_Date"])
    else:
        reviewed = pd.read_csv(reviewed_decisions_csv, sep=";", dtype=str)
        reviewed.columns = [c.strip() for c in reviewed.columns]
        review_cols = ["Subject_ID", "Session_ID", "File_Name",
                       "Reviewed_Action", "Reviewer_Notes", "Review_Date"]
        reviewed = reviewed[[c for c in review_cols if c in reviewed.columns]]
        print(f"Loaded reviewed decisions: {len(reviewed)} rows")

    # Merge (left join: keep all inventory rows)
    df = inventory.merge(reviewed, on=["Subject_ID", "Session_ID", "File_Name"], how="left")

    # Rename inventory Action to avoid ambiguity
    df = df.rename(columns={"Action": "Automated_Action"})

    # --- Session-level flags ---
    has_direct_suvr_s = df.groupby(SESSION_KEY)["Automated_Action"].transform(
        lambda x: (x == ACTION_SUVR).any()
    )
    has_extractable_s = df.groupby(SESSION_KEY)["Automated_Action"].transform(
        lambda x: (x == ACTION_EXTRACT).any()
    )
    has_reviewed_suvr_s = df.groupby(SESSION_KEY)["Reviewed_Action"].transform(
        lambda x: (x.fillna("") == ACTION_SUVR).any()
    )

    # --- Derive per-row status and action ---
    statuses = []
    actions  = []
    for idx, row in df.iterrows():
        status, action = derive_status_action(
            row,
            has_direct_suvr=has_direct_suvr_s[idx],
            has_extractable=has_extractable_s[idx],
            has_reviewed_suvr=has_reviewed_suvr_s[idx],
        )
        statuses.append(status)
        actions.append(action)

    df["Final_Status"]      = statuses
    df["Processing_Action"] = actions

    # --- Column order ---
    output_cols = [
        "Subject_ID", "Session_ID", "File_Name", "Run_Number", "Modality",
        "Frame_Count", "First_Frame_Start_Min", "Last_Frame_End_Min",
        "Total_Scan_Duration_Min", "Frame_Durations_Summary",
        "Injection_Time", "Likely_Procedure", "Automated_Action",
        "Reviewed_Action", "Reviewer_Notes", "Review_Date",
        "Final_Status", "Processing_Action",
    ]
    extra = [c for c in df.columns if c not in output_cols]
    df[output_cols + extra].to_csv(output_csv, index=False)

    # --- Summary ---
    print(f"\nMaster inventory saved to: {output_csv}  ({len(df)} rows)\n")

    print("Final_Status breakdown:")
    for status, count in df["Final_Status"].value_counts().items():
        print(f"  {status}: {count}")

    print("\nProcessing_Action breakdown (included scans only):")
    included = df[df["Final_Status"].str.startswith("Included")]
    if len(included) == 0:
        print("  (none)")
    else:
        for action, count in included["Processing_Action"].value_counts().items():
            print(f"  {action}: {count}")

    print(f"\nNext: python 04_coregister_average.py")


if __name__ == "__main__":
    main()
