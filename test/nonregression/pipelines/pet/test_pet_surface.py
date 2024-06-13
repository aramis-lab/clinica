import shutil
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths

import nibabel as nib
import numpy as np
import pytest

from clinica.utils.pet import SUVRReferenceRegion, Tracer


@pytest.mark.slow
def test_pet_surface(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "PETSurface")
    run_pet_surface(input_dir, tmp_dir, ref_dir, working_dir)


def run_pet_surface(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    tracer = Tracer.FDG
    region = SUVRReferenceRegion.PONS

    parameters = {
        "acq_label": tracer,
        "suvr_reference_region": region,
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

    for hemisphere in ("lh", "rh"):
        for fwhm in (0, 5, 10, 15, 20, 25):
            output_folder = (
                output_dir
                / "caps"
                / "subjects"
                / "sub-ADNI011S4105"
                / "ses-M000"
                / "pet"
                / "surface"
            )
            out = fspath(
                output_folder
                / f"sub-ADNI011S4105_ses-M000_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-{hemisphere}_fwhm-{fwhm}_projection.mgh"
            )
            ref = fspath(
                ref_dir
                / f"sub-ADNI011S4105_ses-M000_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-{hemisphere}_fwhm-{fwhm}_projection.mgh"
            )
            assert np.allclose(
                np.squeeze(nib.load(out).get_fdata(dtype="float32")),
                np.squeeze(nib.load(ref).get_fdata(dtype="float32")),
                rtol=3e-2,
                equal_nan=True,
            )


@pytest.mark.slow
def test_run_pet_surface_longitudinal(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "PETSurfaceLongitudinal"
    )
    run_pet_surface_longitudinal(input_dir, tmp_dir, ref_dir, working_dir)


def run_pet_surface_longitudinal(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    tracer = Tracer.FDG
    region = SUVRReferenceRegion.PONS

    parameters = {
        "acq_label": tracer,
        "suvr_reference_region": region,
        "pvc_psf_tsv": fspath(input_dir / "subjects.tsv"),
        "longitudinal": True,
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

    part_id = "sub-ADNI041S1260"
    sess_id = "ses-M024"
    long_id = "long-M000M006M012M018M024"
    image_id = "_".join((part_id, sess_id, long_id))
    output_folder = (
        output_dir
        / "caps"
        / "subjects"
        / part_id
        / sess_id
        / "pet"
        / long_id
        / "surface_longitudinal"
    )
    for hemisphere in ("lh", "rh"):
        for fwhm in (0, 5, 10, 15, 20, 25):
            out = fspath(
                output_folder
                / f"{image_id}_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-{hemisphere}_fwhm-{fwhm}_projection.mgh"
            )
            ref = fspath(
                ref_dir
                / f"{image_id}_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-{hemisphere}_fwhm-{fwhm}_projection.mgh"
            )
            assert np.allclose(
                np.squeeze(nib.load(out).get_fdata(dtype="float32")),
                np.squeeze(nib.load(ref).get_fdata(dtype="float32")),
                rtol=3e-2,
                equal_nan=True,
            )
