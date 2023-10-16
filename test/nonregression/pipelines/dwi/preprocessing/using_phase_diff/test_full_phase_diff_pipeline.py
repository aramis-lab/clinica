from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths, similarity_measure

import pytest


@pytest.mark.slow
def test_dwi_preprocessing_using_phase_diff_field_map(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir,
        tmp_path,
        "DWIPreprocessingUsingPhaseDiffFieldmap",
    )
    run_dwi_preprocessing_using_phase_diff_field_map(
        input_dir, tmp_dir, ref_dir, working_dir
    )


def run_dwi_preprocessing_using_phase_diff_field_map(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
        "delete_cache": True,
    }
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_file = fspath(
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M000"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"
    )
    ref_file = fspath(
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.95)
