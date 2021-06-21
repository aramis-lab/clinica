# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import pardir
from test.nonregression.testing_tools import *

# Determine location for working_directory
warnings.filterwarnings("ignore")


def test_run_StatisticsSurface(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import numpy as np
    from scipy.io import loadmat

    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import (
        StatisticsSurface,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "StatisticsSurface")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "StatisticsSurface"))
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data": "t1-freesurfer",
        "glm_type": "group_comparison",
        "contrast": "group",
        # Optional parameters
        "covariates": "age sex",
    }
    pipeline = StatisticsSurface(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "StatisticsSurface"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 1}, bypass_check=True)

    # Check files
    filename = "group-UnitTest_AD-lt-CN_measure-ct_fwhm-20_correctedPValue.mat"
    out_file = join(
        root,
        "out",
        "caps",
        "groups",
        "group-UnitTest",
        "statistics",
        "surfstat_group_comparison",
        filename,
    )
    ref_file = join(root, "ref", filename)

    out_file_mat = loadmat(out_file)["correctedpvaluesstruct"]
    ref_file_mat = loadmat(ref_file)["correctedpvaluesstruct"]
    for i in range(4):
        assert np.allclose(
            out_file_mat[0][0][i], ref_file_mat[0][0][i], rtol=1e-8, equal_nan=True
        )
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "StatisticsSurface"), recreate=False)


def test_run_StatisticsVolume(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "StatisticsVolume")

    # Remove potential residual of previous UT
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "StatisticsVolume"), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

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
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "group-UnitTest_covariates.tsv"),
        base_dir=join(working_dir, "StatisticsVolume"),
        parameters=parameters,
    )

    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 2}, bypass_check=True)

    output_t_stat = join(
        root,
        "out",
        "caps",
        "groups",
        "group-UnitTest",
        "statistics_volume",
        "group_comparison_measure-fdg",
        "group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii",
    )
    ref_t_stat = join(
        root,
        "ref",
        "caps",
        "groups",
        "group-UnitTest",
        "statistics_volume",
        "group_comparison_measure-fdg",
        "group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii",
    )

    assert np.allclose(
        nib.load(output_t_stat).get_fdata(), nib.load(ref_t_stat).get_fdata()
    )

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=True)
    clean_folder(join(working_dir, "StatisticsVolume"), recreate=False)


def test_run_StatisticsVolumeCorrection(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import (
        StatisticsVolumeCorrection,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "StatisticsVolumeCorrection")

    # Remove potential residual of previous UT
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "StatisticsVolumeCorrection"), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

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
        caps_directory=join(root, "out", "caps"),
        base_dir=join(working_dir, "StatisticsVolumeCorrection"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)
    compare_folders(join(root, "out"), join(root, "ref"), "caps")

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=True)
    clean_folder(join(working_dir, "StatisticsVolumeCorrection"), recreate=False)
