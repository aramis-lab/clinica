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
def test_workflows_ml(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "WorkflowsML")
    run_workflows_ml(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
def test_spatial_svm(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "SpatialSVM")
    run_spatial_svm(input_dir, tmp_dir, ref_dir, working_dir)


def run_workflows_ml(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import warnings

    from clinica.pipelines.machine_learning.ml_workflows import (
        RegionBasedRepHoldOutLogisticRegression,
        RegionBasedRepHoldOutRandomForest,
        VertexBasedRepHoldOutDualSVM,
        VoxelBasedKFoldDualSVM,
    )
    from clinica.utils.pet import SUVRReferenceRegion, Tracer

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)

    caps_dir = input_dir / "caps"
    tsv = input_dir / "subjects.tsv"
    diagnoses_tsv = input_dir / "diagnosis.tsv"
    group_label = "allADNIdartel"

    output_dir1 = output_dir / "VertexBasedRepHoldOutDualSVM"

    tracer = Tracer.FDG
    region = SUVRReferenceRegion.PONS

    wf1 = VertexBasedRepHoldOutDualSVM(
        caps_directory=fspath(caps_dir),
        subjects_visits_tsv=fspath(tsv),
        diagnoses_tsv=fspath(diagnoses_tsv),
        group_label=group_label,
        output_dir=fspath(output_dir1),
        image_type="PET",
        acq_label=tracer,
        suvr_reference_region=region,
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
        acq_label=tracer,
        suvr_reference_region=region,
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
        acq_label=tracer,
        suvr_reference_region=region,
        fwhm=8,
        n_threads=2,
        n_folds=5,
        grid_search_folds=3,
    )
    wf4.run()


def run_spatial_svm(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil
    from os import fspath

    import nibabel as nib
    from numpy.testing import assert_allclose

    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    pipeline = SpatialSVM(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        group_label="ADNIbl",
        parameters={"orig_input_data_ml": "t1-volume"},
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    for subject in ("sub-ADNI011S0023", "sub-ADNI013S0325"):
        out_data = nib.load(
            fspath(
                caps_dir
                / "subjects"
                / subject
                / "ses-M000"
                / "machine_learning"
                / "input_spatial_svm"
                / "group-ADNIbl"
                / (
                    subject
                    + "_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz"
                )
            )
        ).get_fdata(dtype="float32")
        ref_data = nib.load(
            fspath(
                ref_dir
                / (
                    subject
                    + "_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz"
                )
            )
        ).get_fdata(dtype="float32")
        assert_allclose(out_data, ref_data, rtol=1e-03, atol=1e-7)
