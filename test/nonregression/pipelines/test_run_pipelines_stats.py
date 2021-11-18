# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import *

import pytest

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        "StatisticsSurface",
        "StatisticsVolume",
        # "StatisticsVolumeCorrection",
    ]
)
def test_name(request):
    return request.param


def test_run_stats(cmdopt, tmp_path, test_name):
    import shutil

    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)
    working_dir = Path(cmdopt["wd"])

    if test_name == "StatisticsSurface":
        run_StatisticsSurface(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "StatisticsVolume":
        run_StatisticsVolume(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "StatisticsVolumeCorrection":
        run_StatisticsVolume(input_dir, tmp_out_dir, ref_dir, working_dir)

    else:
        print(f"Test {test_name} not available.")
        assert 0


def run_StatisticsSurface(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import numpy as np
    from scipy.io import loadmat

    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import (
        StatisticsSurface,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data": "t1-freesurfer",
        "glm_type": "group_comparison",
        "contrast": "group",
        # Optional parameters
        "covariates": ["age", "sex"],
    }
    pipeline = StatisticsSurface(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 1}, bypass_check=True)

    # Check files
    filename = "group-UnitTest_AD-lt-CN_measure-ct_fwhm-20_correctedPValue.mat"
    out_file = (
        caps_dir
        / "groups"
        / "group-UnitTest"
        / "statistics"
        / "surfstat_group_comparison"
        / filename
    )
    ref_file = ref_dir / filename

    out_file_mat = loadmat(fspath(out_file))["correctedpvaluesstruct"]
    ref_file_mat = loadmat(fspath(ref_file))["correctedpvaluesstruct"]
    for i in range(4):
        assert np.allclose(
            out_file_mat[0][0][i],
            ref_file_mat[0][0][i],
            rtol=1e-8,
            equal_nan=True,
        )


def run_StatisticsVolume(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "group-UnitTest_covariates.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    # Instantiate pipeline and run()
    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data": "pet-volume",
        "contrast": "group",
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": "FDG",
        "use_pvc_data": False,
        "suvr_reference_region": "pons",
    }

    pipeline = StatisticsVolume(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )

    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 2}, bypass_check=True)

    output_t_stat = (
        caps_dir
        / "groups"
        / "group-UnitTest"
        / "statistics_volume"
        / "group_comparison_measure-fdg"
        / "group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii"
    )
    ref_t_stat = (
        ref_dir
        / "caps"
        / "groups"
        / "group-UnitTest"
        / "statistics_volume"
        / "group_comparison_measure-fdg"
        / "group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii"
    )

    assert np.allclose(
        nib.load(fspath(output_t_stat)).get_fdata(dtype="float32"),
        nib.load(fspath(ref_t_stat)).get_fdata(dtype="float32"),
    )

    # Remove data in out folder


def run_StatisticsVolumeCorrection(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import (
        StatisticsVolumeCorrection,
    )

    caps_dir = output_dir / "caps"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    # Instantiate pipeline and run()
    parameters = {
        "t_map": "group-UnitTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics.nii",
        "height_threshold": 3.2422,
        "FWEp": 4.928,
        "FDRp": 4.693,
        "FWEc": 206987,
        "FDRc": 206987,
        "n_cuts": 15,
    }
    pipeline = StatisticsVolumeCorrection(
        caps_directory=fspath(caps_dir),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)
    compare_folders(output_dir / "caps", ref_dir / "caps", output_dir)
