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


def test_run_DWIPreprocessingUsingT1(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DWIPreprocessingUsingT1")

    # Remove old instance of UT
    clean_folder(join(root, "out", "caps"), recreate=True)
    clean_folder(join(working_dir, "DWIPreprocessingUsingT1"))

    parameters = {"low_bval": 5}
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DWIPreprocessingUsingT1"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Assert
    out_file = join(
        root,
        "out",
        "caps",
        "subjects",
        "sub-CAPP01001TMM",
        "ses-M00",
        "dwi",
        "preprocessing",
        "sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz",
    )
    ref_file = join(
        root, "ref", "sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.97)

    # Delete out/caps folder
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "DWIPreprocessingUsingT1"), recreate=False)


def test_run_DWIPreprocessingUsingPhaseDiffFieldmap(cmdopt):
    import warnings
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFieldmap,
    )

    warnings.filterwarnings("ignore")

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DWIPreprocessingUsingPhaseDiffFieldmap")

    # Remove old instance of UT
    clean_folder(join(root, "out", "caps"))
    clean_folder(join(working_dir, "DWIPreprocessingUsingPhaseDiffFieldmap"))

    parameters = {"low_bval": 5}
    pipeline = DwiPreprocessingUsingPhaseDiffFieldmap(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DWIPreprocessingUsingPhaseDiffFieldmap"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Assert
    out_file = join(
        root,
        "out",
        "caps",
        "subjects",
        "sub-CAPP01001TMM",
        "ses-M00",
        "dwi",
        "preprocessing",
        "sub-CAPP01001TMM_ses-M00_dwi_space-b0_preproc.nii.gz",
    )
    ref_file = join(root, "ref", "sub-CAPP01001TMM_ses-M00_dwi_space-b0_preproc.nii.gz")

    assert similarity_measure(out_file, ref_file, 0.95)

    # Delete out/caps folder
    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(
        join(working_dir, "DWIPreprocessingUsingPhaseDiffFieldmap"), recreate=False
    )


def test_run_DWIDTI(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import numpy as np
    import pandas as pds

    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DWIDTI")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "DWIDTI"))
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    pipeline = DwiDti(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DWIDTI"),
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    subject_id = "sub-CAPP01001TMM"
    maps = ["AD", "FA", "MD", "RD"]
    out_files = [
        join(
            root,
            "out",
            "caps",
            "subjects",
            subject_id,
            "ses-M00",
            "dwi",
            "dti_based_processing",
            "atlas_statistics",
            subject_id
            + "_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-"
            + m
            + "_statistics.tsv",
        )
        for m in maps
    ]
    ref_files = [
        join(
            root,
            "ref",
            subject_id
            + "_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-"
            + m
            + "_statistics.tsv",
        )
        for m in maps
    ]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep="\t")
        out_mean_scalar = np.array(out_csv.mean_scalar)
        ref_csv = pds.read_csv(ref_files[i], sep="\t")
        ref_mean_scalar = np.array(ref_csv.mean_scalar)

        assert np.allclose(out_mean_scalar, ref_mean_scalar, rtol=0.025, equal_nan=True)

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "DWIDTI"), recreate=False)


def test_run_DWIConnectome(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    # Initialization
    working_dir = join(abspath(cmdopt), "DWIConnectome")
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DWIConnectome")
    in_tsv = join(root, "in", "subjects.tsv")
    out_caps_dir = join(root, "out", "caps")

    subject_id = "sub-HMTC20110506MEMEPPAT27"
    session_id = "ses-M00"

    clean_folder(out_caps_dir, recreate=False)
    clean_folder(working_dir)
    shutil.copytree(join(root, "in", "caps"), out_caps_dir)

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=out_caps_dir,
        tsv_file=in_tsv,
        base_dir=working_dir,
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    atlases = ["desikan", "destrieux"]

    out_fod_file = join(
        root,
        "out",
        "caps",
        "subjects",
        subject_id,
        session_id,
        "dwi",
        "connectome_based_processing",
        subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz",
    )
    ref_fod_file = join(
        root,
        "ref",
        subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz",
    )

    out_parc_files = [
        join(
            root,
            "out",
            "caps",
            "subjects",
            subject_id,
            session_id,
            "dwi",
            "connectome_based_processing",
            subject_id
            + "_"
            + session_id
            + "_dwi_space-b0_atlas-"
            + a
            + "_parcellation.nii.gz",
        )
        for a in atlases
    ]
    ref_parc_files = [
        join(
            root,
            "ref",
            subject_id
            + "_"
            + session_id
            + "_dwi_space-b0_atlas-"
            + a
            + "_parcellation.nii.gz",
        )
        for a in atlases
    ]

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for i in range(len(out_parc_files)):
        assert similarity_measure(out_parc_files[i], ref_parc_files[i], 0.955)

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(working_dir, recreate=False)
