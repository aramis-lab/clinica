import shutil
from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure,
)

import numpy as np
import pandas as pd
import pytest


@pytest.mark.slow
def test_dwi_dti(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIDTI")
    run_dwi_dti(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.slow
def test_dwi_connectome(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIConnectome")
    run_dwi_connectome(input_dir, tmp_dir, ref_dir, working_dir)


def run_dwi_dti(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_dti.pipeline import DwiDti
    from clinica.utils.dwi import DTIBasedMeasure

    caps_dir = output_dir / "caps"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    pipeline = DwiDti(
        caps_directory=str(caps_dir),
        tsv_file=str(input_dir / "subjects.tsv"),
        base_dir=str(working_dir),
        parameters={"random_seed": 42},
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    subject_id = "sub-PREVDEMALS0010025PG"
    entities = "ses-M000_space-JHUDTI81_desc-preproc_res-1x1x1"
    output = (
        caps_dir
        / "subjects"
        / subject_id
        / "ses-M000"
        / "dwi"
        / "dti_based_processing"
        / "atlas_statistics"
    )
    for measure in DTIBasedMeasure:
        out_csv = pd.read_csv(
            output / f"{subject_id}_{entities}_map-{measure.value}_statistics.tsv",
            sep="\t",
        )
        ref_csv = pd.read_csv(
            ref_dir / f"{subject_id}_{entities}_map-{measure.value}_statistics.tsv",
            sep="\t",
        )
        assert np.allclose(
            np.array(out_csv.mean_scalar),
            np.array(ref_csv.mean_scalar),
            rtol=0.025,
            equal_nan=True,
        )


def run_dwi_connectome(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_connectome.pipeline import DwiConnectome

    caps_dir = output_dir / "caps"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=str(caps_dir),
        tsv_file=str(input_dir / "subjects.tsv"),
        base_dir=str(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    session_id = "ses-M000"
    subject_id = "sub-PREVDEMALS0010025PG"
    suffix = "dwi_space-b0_model-CSD_diffmodel.nii.gz"
    output_folder = (
        caps_dir
        / "subjects"
        / subject_id
        / session_id
        / "dwi"
        / "connectome_based_processing"
    )
    out_fod_file = output_folder / f"{subject_id}_{session_id}_{suffix}"
    ref_fod_file = ref_dir / f"{subject_id}_{session_id}_{suffix}"

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for atlas in ("desikan", "destrieux"):
        assert similarity_measure(
            output_folder
            / f"{subject_id}_{session_id}_dwi_space-b0_atlas-{atlas}_parcellation.nii.gz",
            ref_dir
            / f"{subject_id}_{session_id}_dwi_space-b0_atlas-{atlas}_parcellation.nii.gz",
            0.955,
        )
