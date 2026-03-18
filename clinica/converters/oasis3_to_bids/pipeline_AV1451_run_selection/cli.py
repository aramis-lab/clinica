"""
CLI entry point for the AV1451 Tau-PET Run Selection Pipeline.

Usage examples:

  # Run a single step with custom directories:
  python -m clinica.converters.oasis3_to_bids.pipeline_AV1451_run_selection.cli \
      evaluate --data-dir /data/raw --inventory-dir /data/inventory

  # Run all steps sequentially:
  python -m clinica.converters.oasis3_to_bids.pipeline_AV1451_run_selection.cli \
      all --data-dir /data/raw --inventory-dir /data/inventory --output-dir /data/processed
"""

import argparse
import sys


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="av1451-pipeline",
        description="AV1451 Tau-PET Run Selection Pipeline for OASIS-3.",
    )

    # --- Shared directory arguments ---
    parser.add_argument(
        "--data-dir",
        help="Path to the raw OASIS-3 data directory (default: D:\\Oasis3\\raw).",
    )
    parser.add_argument(
        "--inventory-dir",
        help="Path to the inventory/output CSV directory (default: D:\\Oasis3\\inventory).",
    )
    parser.add_argument(
        "--output-dir",
        help="Path for processed NIfTI output (step 4 only, default: D:\\Oasis3\\processed\\coreg_avg).",
    )

    subparsers = parser.add_subparsers(dest="step", required=True)

    subparsers.add_parser(
        "evaluate",
        help="Step 1: Evaluate all AV1451 tau-PET scans and generate the inventory CSV.",
    )
    subparsers.add_parser(
        "filter",
        help="Step 2: Identify sessions that need manual review.",
    )
    subparsers.add_parser(
        "master",
        help="Step 3: Build the master inventory by merging automated + reviewed decisions.",
    )
    subparsers.add_parser(
        "coregister",
        help="Step 4: Coregister and average late-phase tau-PET frames.",
    )
    subparsers.add_parser(
        "all",
        help="Run steps 1-3 sequentially (step 4 requires explicit invocation).",
    )

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    step = args.step

    if step in ("evaluate", "all"):
        from . import _01_evaluate_tau_scans as step1
        print("=" * 60)
        print("STEP 1: Evaluate tau scans")
        print("=" * 60)
        step1.main(data_dir=args.data_dir, output_dir=args.inventory_dir)

    if step in ("filter", "all"):
        from . import _02_filter_needs_review as step2
        print("\n" + "=" * 60)
        print("STEP 2: Filter & flag for manual review")
        print("=" * 60)
        step2.main(inventory_dir=args.inventory_dir)

    if step in ("master", "all"):
        from . import _03_build_master_inventory as step3
        print("\n" + "=" * 60)
        print("STEP 3: Build master inventory")
        print("=" * 60)
        step3.main(inventory_dir=args.inventory_dir)

    if step == "coregister":
        from . import _04_coregister_average as step4
        print("=" * 60)
        print("STEP 4: Coregister and average")
        print("=" * 60)
        step4.main(
            data_dir=args.data_dir,
            inventory_dir=args.inventory_dir,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
