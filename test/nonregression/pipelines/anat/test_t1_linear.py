from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import (
    compare_folders,
    compare_niftis,
    configure_paths,
)

import pytest


@pytest.mark.fast
def test_t1_linear(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1Linear")
    run_t1_linear(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
def test_flair_linear(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "FlairLinear")
    run_flair_linear(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_linear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    out_caps = output_dir / "caps"
    ref_caps = ref_dir / "caps"

    pipeline = AnatLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(out_caps),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={"uncropped_image": False},
        name="t1-linear",
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(out_caps, ref_caps, output_dir)
    compare_niftis(out_caps, ref_caps)


def run_flair_linear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    out_caps = output_dir / "caps"
    ref_caps = ref_dir / "caps"

    pipeline = AnatLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(out_caps),
        base_dir=fspath(working_dir),
        parameters={"uncropped_image": False},
        name="flair-linear",
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(out_caps, ref_caps, output_dir)
    compare_niftis(out_caps, ref_caps)
