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


def test_run_PETVolume(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume

    working_dir = cmdopt["wd"]
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "PETVolume")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETVolume"))
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {
        "group_label": "UnitTest",
        "acq_label": "fdg",
        "suvr_reference_region": "pons",
    }
    pipeline = PETVolume(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "PETVolume"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_files = [
        join(
            root,
            "out/caps/subjects/" + sub + "/ses-M00/pet/preprocessing/group-UnitTest",
            sub
            + "_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz",
        )
        for sub in subjects
    ]
    ref_files = [
        join(
            root,
            "ref",
            sub
            + "_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz",
        )
        for sub in subjects
    ]

    for i in range(len(out_files)):
        assert likeliness_measure(
            out_files[i], ref_files[i], (1e-2, 0.25), (1e-1, 0.001)
        )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETVolume"), recreate=False)


def test_run_PETLinear(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.pet_linear.pet_linear_pipeline import PETLinear

    working_dir = cmdopt["wd"]
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "PETLinear")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETLinear"))
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {"acq_label": "fdg", "suvr_reference_region": "pons"}

    # Instantiate pipeline
    pipeline = PETLinear(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "PETLinear"),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    out_folder = join(root, "out")
    ref_folder = join(root, "ref")

    compare_folders(out_folder, ref_folder, tmp_path="caps")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETLinear"), recreate=False)


def test_run_PETSurfaceCrossSectional(cmdopt):
    import shutil
    from os.path import abspath, dirname, join

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    working_dir = cmdopt["wd"]
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "PETSurface")

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETSurface"))
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))

    parameters = {
        "acq_label": "FDG",
        "suvr_reference_region": "pons",
        "pvc_psf_tsv": join(root, "in", "subjects.tsv"),
        "longitudinal": False,
    }
    pipeline = PetSurface(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "out", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        base_dir=join(working_dir, "PETSurface"),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    # Check files
    out_files = [
        join(
            root,
            "out/caps/subjects/sub-ADNI011S4105/ses-M00/pet/surface",
            "sub-ADNI011S4105_ses-M00_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-"
            + h
            + "_fwhm-"
            + str(f)
            + "_projection.mgh",
        )
        for h in ["lh", "rh"]
        for f in [0, 5, 10, 15, 20, 25]
    ]
    ref_files = [
        join(
            root,
            "ref",
            "sub-ADNI011S4105_ses-M00_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-"
            + h
            + "_fwhm-"
            + str(f)
            + "_projection.mgh",
        )
        for h in ["lh", "rh"]
        for f in [0, 5, 10, 15, 20, 25]
    ]

    for i in range(len(out_files)):
        assert np.allclose(
            np.squeeze(nib.load(out_files[i]).get_fdata(dtype="float32")),
            np.squeeze(nib.load(ref_files[i]).get_fdata(dtype="float32")),
            rtol=3e-2,
            equal_nan=True,
        )

    clean_folder(join(root, "out", "caps"), recreate=False)
    clean_folder(join(working_dir, "PETSurface"), recreate=False)


# def test_run_PETSurfaceLongitudinal(cmdopt):
#     from os.path import dirname, join, abspath
#     import shutil
#     import nibabel as nib
#     import numpy as np
#     from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
#
#     working_dir = cmdopt
#     root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
#     root = join(root, 'data', 'PETSurfaceLongitudinal')
#
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'PETSurfaceLongitudinal'))
#     shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))
#
#     parameters = {
#         'acq_label': 'FDG',
#         'suvr_reference_region': 'pons',
#         'pvc_psf_tsv': join(root, 'in', 'subjects.tsv'),
#         'longitudinal': True
#     }
#     pipeline = PetSurface(
#         bids_directory=join(root, 'in', 'bids'),
#         caps_directory=join(root, 'out', 'caps'),
#         tsv_file=join(root, 'in', 'subjects.tsv'),
#         base_dir=join(working_dir, 'PETSurfaceLongitudinal'),
#         parameters=parameters
#     )
#     pipeline.build()
#     pipeline.run(bypass_check=True)
#
#     # Check files
#     part_id = 'sub-ADNI041S1260'
#     sess_id = 'ses-M24'
#     long_id = 'long-M00M06M12M18M24'
#     image_id = part_id + '_' + sess_id + '_' + long_id
#     out_files = [join(root, 'out', 'caps', 'subjects', part_id, sess_id, 'pet', long_id, 'surface_longitudinal',
#                       image_id + '_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
#                       + h + '_fwhm-' + str(f) + '_projection.mgh')
#                  for h in ['lh', 'rh']
#                  for f in [0, 5, 10, 15, 20, 25]]
#     ref_files = [join(root, 'ref',
#                       image_id + '_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
#                       + h + '_fwhm-' + str(f) + '_projection.mgh')
#                  for h in ['lh', 'rh']
#                  for f in [0, 5, 10, 15, 20, 25]]
#
#     # Tolerance values were taken from PETSurface - Cross-sectional case
#     for i in range(len(out_files)):
#         assert np.allclose(np.squeeze(nib.load(out_files[i]).get_fdata(dtype="float32")),
#                            np.squeeze(nib.load(ref_files[i]).get_fdata(dtype="float32")),
#                            rtol=3e-2, equal_nan=True)
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
