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


def test_run_WorkflowsML(cmdopt):
    import shutil
    import warnings
    from os.path import abspath, dirname, join

    from clinica.pipelines.machine_learning.ml_workflows import (
        RegionBasedRepHoldOutLogisticRegression,
        RegionBasedRepHoldOutRandomForest,
        VertexBasedRepHoldOutDualSVM,
        VoxelBasedKFoldDualSVM,
    )

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)

    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "WorkflowsML")
    root_input = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root_input = join(root_input, "data", "InputsML")

    caps_dir = join(root_input, "in", "caps")
    tsv = join(root_input, "in", "subjects.tsv")
    diagnoses_tsv = join(root_input, "in", "diagnosis.tsv")
    group_label = "allADNIdartel"

    output_dir1 = join(root, "out", "VertexBasedRepHoldOutDualSVM")
    clean_folder(output_dir1, recreate=True)
    wf1 = VertexBasedRepHoldOutDualSVM(
        caps_directory=caps_dir,
        subjects_visits_tsv=tsv,
        diagnoses_tsv=diagnoses_tsv,
        group_label=group_label,
        output_dir=output_dir1,
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
    shutil.rmtree(output_dir1)

    output_dir2 = join(root, "out", "RegionBasedRepHoldOutLogisticRegression")
    clean_folder(output_dir2, recreate=True)
    wf2 = RegionBasedRepHoldOutLogisticRegression(
        caps_directory=caps_dir,
        subjects_visits_tsv=tsv,
        diagnoses_tsv=diagnoses_tsv,
        group_label=group_label,
        image_type="PET",
        atlas="AICHA",
        output_dir=output_dir2,
        acq_label="fdg",
        suvr_reference_region="pons",
        use_pvc_data=False,
        n_threads=2,
        n_iterations=10,
        grid_search_folds=3,
        test_size=0.3,
    )
    wf2.run()
    shutil.rmtree(output_dir2)

    output_dir3 = join(root, "out", "RegionBasedRepHoldOutRandomForest")
    clean_folder(output_dir3, recreate=True)
    wf3 = RegionBasedRepHoldOutRandomForest(
        caps_directory=caps_dir,
        subjects_visits_tsv=tsv,
        diagnoses_tsv=diagnoses_tsv,
        group_label=group_label,
        image_type="T1w",
        atlas="AAL2",
        output_dir=output_dir3,
        n_threads=2,
        n_iterations=10,
        grid_search_folds=3,
        test_size=0.3,
    )
    wf3.run()
    shutil.rmtree(output_dir3)

    output_dir4 = join(root, "out", "VoxelBasedKFoldDualSVM")
    clean_folder(output_dir4, recreate=True)
    wf4 = VoxelBasedKFoldDualSVM(
        caps_directory=caps_dir,
        subjects_visits_tsv=tsv,
        diagnoses_tsv=diagnoses_tsv,
        group_label=group_label,
        image_type="PET",
        output_dir=output_dir4,
        acq_label="fdg",
        suvr_reference_region="pons",
        fwhm=8,
        n_threads=2,
        n_folds=5,
        grid_search_folds=3,
    )
    wf4.run()
    shutil.rmtree(output_dir4)


def test_run_SpatialSVM(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "SpatialSVM")

    # Remove potential residual of previous UT
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "SpatialSVM"), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {"group_label": "ADNIbl", "orig_input_data": "t1-volume"}
    # Instantiate pipeline and run()
    pipeline = SpatialSVM(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "SpatialSVM"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = ["sub-ADNI011S0023", "sub-ADNI013S0325"]
    out_data_REG_NIFTI = [
        nib.load(
            join(
                root,
                "out",
                "caps",
                "subjects",
                sub,
                "ses-M00",
                "machine_learning",
                "input_spatial_svm",
                "group-ADNIbl",
                sub
                + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz",
            )
        ).get_fdata()
        for sub in subjects
    ]
    ref_data_REG_NIFTI = [
        nib.load(
            join(
                root,
                "ref",
                sub
                + "_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz",
            )
        ).get_fdata()
        for sub in subjects
    ]
    for i in range(len(out_data_REG_NIFTI)):
        assert np.allclose(
            out_data_REG_NIFTI[i], ref_data_REG_NIFTI[i], rtol=1e-3, equal_nan=True
        )

    # Remove data in out folder
    clean_folder(join(root, "out", "caps"), recreate=True)
    clean_folder(join(working_dir, "SpatialSVM"), recreate=False)
