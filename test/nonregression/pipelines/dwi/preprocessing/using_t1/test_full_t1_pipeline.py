"""Non regression test for the full DWIPreprocessingUsingT1 pipeline."""

from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure_large_nifti,
)

import pytest


@pytest.mark.slow
def test_dwi_preprocessing_using_t1(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreprocessingUsingT1"
    )
    run_dwi_preprocessing_using_t1(input_dir, tmp_dir, ref_dir, working_dir)


def run_dwi_preprocessing_using_t1(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
        "random_seed": 42,
    }
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    output_folder = (
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M000"
        / "dwi"
        / "preprocessing"
    )
    out_preprocessed_image = (
        output_folder / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-T1w_preproc.nii.gz"
    )
    ref_preprocessed_image = (
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-T1w_preproc.nii.gz"
    )

    assert similarity_measure_large_nifti(
        out_preprocessed_image, ref_preprocessed_image, 0.99
    )
