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

    # Prepare test for different parameters

    modalities = ["t1-linear", "pet-linear", "custom"]
    uncropped_image = [True, False]

    image_params = {"extract_method": "image"}
    slice_params = {"extract_method": "patch", "patch_size": 50, "stride_size": 50}
    patch_params = {
        "extract_method": "slice",
        "slice_mode": "rgb",
        "slice_direction": 0,
    }
    roi_params = {
        "extract_method": "roi",
        "roi_list": ["rightHippocampusBox", "leftHippocampusBox"],
        "roi_uncrop_output": False,
    }

    data = [image_params, slice_params, patch_params, roi_params]

    for parameters in data:
        for modality in modalities:
            if modality == "pet-linear":
                parameters["acq_label"] = "av45"
                parameters["suvr_reference_region"] = "pons2"
                parameters["use_uncropped_image"] = False
                DLPrepareData_Generic(root, working_dir, parameters)
            elif modality == "custom":
                parameters["custom_template"] = "Ixi549Space"
                parameters[
                    "custom_suffix"
                ] = "graymatter_space-Ixi549Space_modulated-off_probability.nii.gz"
            elif modality == "t1-linear":
                for flag in uncropped_image:
                    parameters["modality"] = modality
                    parameters["use_uncropped_image"] = flag
                    DLPrepareData_Generic(root, working_dir, parameters)
            else:
                raise NotImplementedError(
                    f"Test for modality {modality} was not implemented."
                )

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
