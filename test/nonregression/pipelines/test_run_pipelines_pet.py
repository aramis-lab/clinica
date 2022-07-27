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

from clinica.utils.pet import Tracer

warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        "PETVolume",
        "PETLinear",
        "PETSurfaceCrossSectional",
    ]
)
def test_name(request):
    return request.param


def run_PETVolume(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume

    # Arrange
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    tracer = Tracer.FDG

    parameters = {
        "group_label": "UnitTest",
        "acq_label": tracer,
        "suvr_reference_region": "pons",
        "skip_question": False,
    }
    pipeline = PETVolume(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    # Acte
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Assert
    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_files = [
        fspath(
            output_dir
            / "caps"
            / "subjects"
            / sub
            / "ses-M00/pet/preprocessing/group-UnitTest"
            / (
                sub
                + f"_ses-M00_trc-{tracer}_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz"
            )
        )
        for sub in subjects
    ]
    ref_files = [
        fspath(
            ref_dir
            / (
                sub
                + f"_ses-M00_trc-{tracer}_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz"
            )
        )
        for sub in subjects
    ]

    for i in range(len(out_files)):
        assert likeliness_measure(
            out_files[i], ref_files[i], (1e-2, 0.25), (1e-1, 0.001)
        )


def run_PETLinear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:

    import shutil

    from clinica.pipelines.pet_linear.pet_linear_pipeline import PETLinear

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    parameters = {"acq_label": Tracer.FDG, "suvr_reference_region": "pons"}

    # Instantiate pipeline
    pipeline = PETLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(output_dir / "caps", ref_dir / "caps", output_dir)


def run_PETSurfaceCrossSectional(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    tracer = Tracer.FDG

    parameters = {
        "acq_label": tracer,
        "suvr_reference_region": "pons",
        "pvc_psf_tsv": fspath(input_dir / "subjects.tsv"),
        "longitudinal": False,
        "skip_question": False,
    }
    pipeline = PetSurface(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    # Check files
    out_files = [
        (
            output_dir
            / "caps/subjects/sub-ADNI011S4105/ses-M00/pet/surface"
            / (
                f"sub-ADNI011S4105_ses-M00_trc-{tracer}_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-"
                + h
                + "_fwhm-"
                + str(f)
                + "_projection.mgh"
            )
        )
        for h in ["lh", "rh"]
        for f in [0, 5, 10, 15, 20, 25]
    ]
    ref_files = [
        (
            ref_dir
            / (
                f"sub-ADNI011S4105_ses-M00_trc-{tracer}_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-"
                + h
                + "_fwhm-"
                + str(f)
                + "_projection.mgh"
            )
        )
        for h in ["lh", "rh"]
        for f in [0, 5, 10, 15, 20, 25]
    ]

    for i in range(len(out_files)):
        assert np.allclose(
            np.squeeze(nib.load(fspath(out_files[i])).get_fdata(dtype="float32")),
            np.squeeze(nib.load(fspath(ref_files[i])).get_fdata(dtype="float32")),
            rtol=1e-1,
            equal_nan=True,
        )


def test_run_pet(cmdopt, tmp_path, test_name):
    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)
    working_dir = Path(cmdopt["wd"])

    if test_name == "PETVolume":
        run_PETVolume(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "PETLinear":
        run_PETLinear(input_dir, tmp_out_dir, ref_dir, working_dir)

    elif test_name == "PETSurfaceCrossSectional":
        run_PETSurfaceCrossSectional(
            base_dir / "PETSurface" / "in",
            tmp_out_dir,
            base_dir / "PETSurface" / "ref",
            working_dir,
        )
    else:
        print(f"Test {test_name} not available.")
        assert 0


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
#     tracer = Tracer.FDG
#
#     parameters = {
#         'acq_label': Tracer.FDG,
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
#                       image_id + f'_trc-{tracer}_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
#                       + h + '_fwhm-' + str(f) + '_projection.mgh')
#                  for h in ['lh', 'rh']
#                  for f in [0, 5, 10, 15, 20, 25]]
#     ref_files = [join(root, 'ref',
#                       image_id + '_trc-{tracer}_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
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
