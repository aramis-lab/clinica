# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import fspath
from test.nonregression.testing_tools import *

import pytest

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        "DWIPreprocessingUsingT1",
        "DWIPreprocessingUsingPhaseDiffFieldmap",
        "DWIDTI",
        "DWIConnectome",
    ]
)
def test_name(request):
    return request.param


def test_run_dwi(cmdopt, tmp_path, test_name):
    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)
    working_dir = Path(cmdopt["wd"])

    if test_name == "DWIPreprocessingUsingT1":
        run_DWIPreprocessingUsingT1(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "DWIPreprocessingUsingPhaseDiffFieldmap":
        run_DWIPreprocessingUsingPhaseDiffFieldmap(
            input_dir, tmp_out_dir, ref_dir, working_dir
        )

    elif test_name == "DWIDTI":
        run_DWIDTI(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "DWIConnectome":
        run_DWIConnectome(input_dir, tmp_out_dir, ref_dir, working_dir)

    else:
        print(f"Test {test_name} not available.")
        assert 0


def run_DWIPreprocessingUsingT1(
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

    # Assert
    out_file = fspath(
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M00"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG-ses-M00_dwi_space-T1w_preproc.nii.gz"
    )
    ref_file = fspath(
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M00_dwi_space-T1w_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.97)


def run_DWIPreprocessingUsingPhaseDiffFieldmap(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import warnings

    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    warnings.filterwarnings("ignore")

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
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

    # Assert
    out_file = fspath(
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M00"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG_ses-M00_dwi_space-b0_preproc.nii.gz"
    )
    ref_file = fspath(
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M00_dwi_space-b0_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.95)


def run_DWIDTI(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import numpy as np
    import pandas as pds

    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    pipeline = DwiDti(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    subject_id = "sub-PREVDEMALS0010025PG"
    maps = ["AD", "FA", "MD", "RD"]
    out_files = [
        fspath(
            caps_dir
            / "subjects"
            / subject_id
            / "ses-M00"
            / "dwi"
            / "dti_based_processing"
            / "atlas_statistics"
            / (
                subject_id
                + "_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-"
                + m
                + "_statistics.tsv"
            )
        )
        for m in maps
    ]
    ref_files = [
        fspath(
            ref_dir
            / (
                subject_id
                + "_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-"
                + m
                + "_statistics.tsv"
            )
        )
        for m in maps
    ]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep="\t")
        out_mean_scalar = np.array(out_csv.mean_scalar)
        ref_csv = pds.read_csv(ref_files[i], sep="\t")
        ref_mean_scalar = np.array(ref_csv.mean_scalar)

        assert np.allclose(out_mean_scalar, ref_mean_scalar, rtol=0.025, equal_nan=True)


def run_DWIConnectome(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    atlases = ["desikan", "destrieux"]
    subject_id = "sub-PREVDEMALS0010025PG"
    session_id = "ses-M00"

    out_fod_file = fspath(
        caps_dir
        / "subjects"
        / subject_id
        / session_id
        / "dwi"
        / "connectome_based_processing"
        / (subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz")
    )
    ref_fod_file = fspath(
        ref_dir
        / (subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz"),
    )

    out_parc_files = [
        fspath(
            caps_dir
            / "subjects"
            / subject_id
            / session_id
            / "dwi"
            / "connectome_based_processing"
            / (
                subject_id
                + "_"
                + session_id
                + "_dwi_space-b0_atlas-"
                + a
                + "_parcellation.nii.gz"
            )
        )
        for a in atlases
    ]
    ref_parc_files = [
        fspath(
            ref_dir
            / (
                subject_id
                + "_"
                + session_id
                + "_dwi_space-b0_atlas-"
                + a
                + "_parcellation.nii.gz"
            )
        )
        for a in atlases
    ]

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for i in range(len(out_parc_files)):
        assert similarity_measure(out_parc_files[i], ref_parc_files[i], 0.955)
