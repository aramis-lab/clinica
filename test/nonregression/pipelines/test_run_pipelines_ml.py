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
        "WorkflowsML",
        "SpatialSVM",
    ]
)
def test_name(request):
    return request.param


def test_run_ml(cmdopt, tmp_path, test_name):
    import shutil

    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)
    working_dir = Path(cmdopt["wd"])

    if test_name == "WorkflowsML":
        run_WorkflowsML(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "SpatialSVM":
        run_SpatialSVM(input_dir, tmp_out_dir, ref_dir, working_dir)

    else:
        print(f"Test {test_name} not available.")
        assert 0


def run_WorkflowsML(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import warnings

    from clinica.pipelines.machine_learning.ml_workflows import (
        RegionBasedRepHoldOutLogisticRegression,
        RegionBasedRepHoldOutRandomForest,
        VertexBasedRepHoldOutDualSVM,
        VoxelBasedKFoldDualSVM,
    )

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)

    caps_dir = input_dir / "caps"
    tsv = input_dir / "subjects.tsv"
    diagnoses_tsv = input_dir / "diagnosis.tsv"
    group_label = "allADNIdartel"

    output_dir1 = output_dir / "VertexBasedRepHoldOutDualSVM"

    wf1 = VertexBasedRepHoldOutDualSVM(
        caps_directory=fspath(caps_dir),
        subjects_visits_tsv=fspath(tsv),
        diagnoses_tsv=fspath(diagnoses_tsv),
        group_label=group_label,
        output_dir=fspath(output_dir1),
        image_type="PET",
        acq_label="fdg",
        suvr_reference_region="pons",
        fwhm=20,
        n_threads=2,
        n_iterations=10,
        grid_search_folds=3,
        test_size=0.3,
    )
    wf1.run()

    output_dir2 = output_dir / "RegionBasedRepHoldOutLogisticRegression"

    wf2 = RegionBasedRepHoldOutLogisticRegression(
        caps_directory=fspath(caps_dir),
        subjects_visits_tsv=fspath(tsv),
        diagnoses_tsv=fspath(diagnoses_tsv),
        group_label=group_label,
        image_type="PET",
        atlas="AICHA",
        output_dir=fspath(output_dir2),
        acq_label="fdg",
        suvr_reference_region="pons",
        use_pvc_data=False,
        n_threads=2,
        n_iterations=10,
        grid_search_folds=3,
        test_size=0.3,
    )
    wf2.run()

    output_dir3 = output_dir / "RegionBasedRepHoldOutRandomForest"

    wf3 = RegionBasedRepHoldOutRandomForest(
        caps_directory=fspath(caps_dir),
        subjects_visits_tsv=fspath(tsv),
        diagnoses_tsv=fspath(diagnoses_tsv),
        group_label=group_label,
        image_type="T1w",
        atlas="AAL2",
        output_dir=fspath(output_dir3),
        n_threads=2,
        n_iterations=10,
        grid_search_folds=3,
        test_size=0.3,
    )
    wf3.run()

    output_dir4 = output_dir / "VoxelBasedKFoldDualSVM"

    wf4 = VoxelBasedKFoldDualSVM(
        caps_directory=fspath(caps_dir),
        subjects_visits_tsv=fspath(tsv),
        diagnoses_tsv=fspath(diagnoses_tsv),
        group_label=group_label,
        image_type="PET",
        output_dir=fspath(output_dir4),
        acq_label="fdg",
        suvr_reference_region="pons",
        fwhm=8,
        n_threads=2,
        n_folds=5,
        grid_search_folds=3,
    )
    wf4.run()


def run_SpatialSVM(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:

    import shutil
    from os import fspath

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {"group_label": "ADNIbl", "orig_input_data": "t1-volume"}
    # Instantiate pipeline and run()
    pipeline = SpatialSVM(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = ["sub-ADNI011S0023", "sub-ADNI013S0325"]
    out_data_REG_NIFTI = [
        nib.load(
            caps_dir
            / "subjects"
            / sub
            / "ses-M00"
            / "machine_learning"
            / "input_spatial_svm"
            / "group-ADNIbl"
            / (
                sub
                + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz"
            )
        ).get_fdata(dtype="float32")
        for sub in subjects
    ]
    ref_data_REG_NIFTI = [
        nib.load(
            ref_dir
            / (
                sub
                + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz"
            )
        ).get_fdata(dtype="float32")
        for sub in subjects
    ]
    for i in range(len(out_data_REG_NIFTI)):
        assert np.allclose(
            out_data_REG_NIFTI[i], ref_data_REG_NIFTI[i], rtol=1e-3, equal_nan=True
        )
