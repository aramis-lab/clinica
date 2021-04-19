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

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "DeepLearningPrepareData")

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, "DeepLearningPrepareData"))
    clean_folder(join(root, "out", "caps"), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    # Test the patch extraction
    data= [
        {"modality": "t1-linear", "use_uncropped_image": False, "extract_method": "image"},
        {"modality": "t1-linear", "use_uncropped_image": True, "extract_method": "image"},
        {"modality": "t1-linear", "extract_method": "patch", "use_uncropped_image": False, "patch_size": 50, "stride_size": 50},
        {"modality": "t1-linear", "extract_method": "patch", "use_uncropped_image": True, "patch_size": 50, "stride_size": 50},
        {"modality": "t1-linear", "extract_method": "slice", "use_uncropped_image": False, "slice_mode": "rgb", "slice_direction": 0},
        {"modality": "t1-linear", "extract_method": "slice", "use_uncropped_image": True, "slice_mode": "rgb", "slice_direction": 0},
    ]
    for parameters in data:
        DLPrepareData_Generic(root, working_dir, parameters)

    # Check output vs ref
    out_folder = join(root, "out")
    ref_folder = join(root, "ref")

    compare_folders(out_folder, ref_folder, shared_folder_name="caps")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "DeepLearningPrepareData"), recreate=False)



def DLPrepareData_Generic(root, working_dir, parameters):

    from os.path import join
    from clinica.pipelines.deeplearning_prepare_data.deeplearning_prepare_data_pipeline import (
        DeepLearningPrepareData,
    )

    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "DeepLearningPrepareData"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)
