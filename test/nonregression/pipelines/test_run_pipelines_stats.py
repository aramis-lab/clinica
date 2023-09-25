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


@pytest.mark.fast
def test_statistics_surface(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsSurface"
    )
    run_statistics_surface(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
@pytest.mark.xdist_group(name="test-group-using-matlab")
def test_statistics_volume_pet(cmdopt, tmp_path):
    """Test the StatisticsVolume pipeline with inputs from PETVolume.

    This test should run in the same process as test_statistics_volume_t1
    to avoid MATLAB race conditions.
    """
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsVolume"
    )
    run_statistics_volume_pet(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
@pytest.mark.xdist_group(name="test-group-using-matlab")
def test_statistics_volume_t1(cmdopt, tmp_path):
    """Test the StatisticsVolume pipeline with inputs from T1Volume.

    This test should run in the same process as test_statistics_volume_pet
    to avoid MATLAB race conditions.
    """
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsVolume"
    )
    run_statistics_volume_t1(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
def test_statistics_volume_correction(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsVolumeCorrection"
    )
    run_statistics_volume_correction(input_dir, tmp_dir, ref_dir, working_dir)


def run_statistics_surface(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import numpy as np
    from scipy.io import loadmat

    from clinica.pipelines.statistics_surface.pipeline import StatisticsSurface

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {
        "group_label": "UnitTest",
        "orig_input_data": "t1-freesurfer",
        "glm_type": "group_comparison",
        "contrast": "group",
        "covariates": ["age", "sex"],
    }
    pipeline = StatisticsSurface(
        caps_directory=caps_dir,
        tsv_file=tsv,
        base_dir=working_dir,
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 1}, bypass_check=True)

    # Check files
    for contrast in ("AD-lt-CN", "CN-lt-AD"):
        for suffix, struct in zip(
            ["coefficients", "uncorrectedPValue", "FDR", "correctedPValue"],
            ["coef", "uncorrectedpvaluesstruct", "FDR", "correctedpvaluesstruct"],
        ):
            filename = f"group-UnitTest_{contrast}_measure-ct_fwhm-20_{suffix}.mat"
            out_file = (
                caps_dir
                / "groups"
                / "group-UnitTest"
                / "statistics"
                / "surfstat_group_comparison"
                / filename
            )
            ref_file = ref_dir / filename
            out_file_mat = loadmat(out_file)[struct]
            ref_file_mat = loadmat(ref_file)[struct]
            if suffix in ["coefficients", "FDR"]:
                assert np.allclose(
                    out_file_mat, ref_file_mat, rtol=1e-8, equal_nan=True
                )
            else:
                keys_to_compare = ["P", "mask", "thresh"]
                if suffix == "correctedPValue":
                    keys_to_compare.append("C")
                for i in keys_to_compare:
                    assert np.allclose(
                        out_file_mat[0][0][i],
                        ref_file_mat[0][0][i],
                        rtol=1e-8,
                        equal_nan=True,
                    )


def run_statistics_volume_pet(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )
    from clinica.utils.pet import Tracer

    caps_dir = output_dir / "caps"
    tsv = input_dir / "group-UnitTest_covariates.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    # Instantiate pipeline and run()
    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data_volume": "pet-volume",
        "contrast": "group",
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": Tracer.FDG,
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


def run_statistics_volume_t1(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )
    from clinica.utils.pet import Tracer

    caps_dir = output_dir / "caps"
    tsv = input_dir / "group-UnitTest_covariates.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    # Instantiate pipeline and run()
    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data_volume": "t1-volume",
        "contrast": "group",
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
        / "group_comparison_measure-graymatter"
        / "group-UnitTest_AD-lt-CN_measure-graymatter_fwhm-8_TStatistics.nii"
    )
    ref_t_stat = (
        ref_dir
        / "caps"
        / "groups"
        / "group-UnitTest"
        / "statistics_volume"
        / "group_comparison_measure-graymatter"
        / "group-UnitTest_AD-lt-CN_measure-graymatter_fwhm-8_TStatistics.nii"
    )

    assert np.allclose(
        nib.load(fspath(output_t_stat)).get_fdata(dtype="float32"),
        nib.load(fspath(ref_t_stat)).get_fdata(dtype="float32"),
    )


def run_statistics_volume_correction(
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
    compare_folders(output_dir / "caps/groups", ref_dir / "caps/groups", output_dir)
