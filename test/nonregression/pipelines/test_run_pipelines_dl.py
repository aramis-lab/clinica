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


def test_run_DLPrepareData(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.deeplearning_prepare_data.deeplearning_prepare_data_pipeline import (
        DeepLearningPrepareData,
    )

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DeepLearningPrepareData")

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, "DeepLearningPrepareData"))
    clean_folder(join(root, "out", "caps"), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    # Test the transformation of the complete T1 MRI
    parameters = {"modality": "t1-linear", "extract_method": "image"}
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DeepLearningPrepareData"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Test the patch extraction
    parameters = {
        "modality": "t1-linear",
        "extract_method": "patch",
        "patch_size": 50,
        "stride_size": 50,
    }
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DeepLearningPrepareData"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Test the slice extraction
    parameters = {
        "modality": "t1-linear",
        "extract_method": "slice",
        "slice_mode": "rgb",
        "slice_direction": 0,
    }
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DeepLearningPrepareData"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)
    # Check output vs ref

    out_folder = join(root, "out")
    ref_folder = join(root, "ref")

    compare_folders(out_folder, ref_folder, shared_folder_name="caps")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "DeepLearningPrepareData"), recreate=False)
