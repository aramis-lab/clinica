from pathlib import Path
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal


@pytest.mark.parametrize("suffix", ["T1w", "FLAIR", "fooo"])
def test_get_substitutions_datasink(suffix):
    from clinica.pipelines.t1_linear.anat_linear_utils import (
        _get_substitutions_datasink,
    )

    bids_image_id = f"sub-ADNI022S0004_ses-M000_{suffix}"

    substitutions = _get_substitutions_datasink(bids_image_id, suffix)

    assert len(substitutions) == 3
    assert substitutions[0] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped_cropped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_{suffix}.nii.gz",
    )
    assert substitutions[1] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}0GenericAffine.mat",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
    )
    assert substitutions[2] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_res-1x1x1_{suffix}.nii.gz",
    )


def n4biasfieldcorrection_mock(
    input_image: Path,
    bspline_fitting_distance: int,
    save_bias: bool = False,
    verbose: bool = False,
):
    """The mock simply returns the input image without any processing."""
    return nib.load(input_image)


def test_run_n4biasfieldcorrection_no_bias_saving(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_utils import run_n4biasfieldcorrection

    data = np.random.random((10, 10, 10))
    nib.save(nib.Nifti1Image(data, np.eye(4)), tmp_path / "test.nii.gz")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with patch("ants.image_write", wraps=nib.save) as image_write_mock:
        with patch(
            "clinica.pipelines.t1_linear.anat_linear_utils._call_n4_bias_field_correction",
            wraps=n4biasfieldcorrection_mock,
        ) as ants_bias_correction_mock:
            bias_corrected_image = run_n4biasfieldcorrection(
                tmp_path / "test.nii.gz",
                bspline_fitting_distance=300,
                output_prefix="sub-01_ses-M000",
                output_dir=output_dir,
            )
            image_write_mock.assert_called_once()
            ants_bias_correction_mock.assert_called_once_with(
                tmp_path / "test.nii.gz",
                300,
                save_bias=False,
                verbose=False,
            )
    # Verify that the bias corrected image exists
    # If all went well, it will be the same as the input image because of the mocks.
    assert [f.name for f in output_dir.iterdir()] == [
        "sub-01_ses-M000_bias_corrected_image.nii.gz"
    ]
    assert bias_corrected_image.exists()
    bias_corrected_nifti = nib.load(bias_corrected_image)
    assert_array_equal(bias_corrected_nifti.affine, np.eye(4))
    assert_array_equal(bias_corrected_nifti.get_fdata(), data)


def test_run_n4biasfieldcorrection(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_utils import run_n4biasfieldcorrection

    data = np.random.random((10, 10, 10))
    nib.save(nib.Nifti1Image(data, np.eye(4)), tmp_path / "test.nii.gz")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with patch("ants.image_write", wraps=nib.save) as image_write_mock:
        with patch(
            "clinica.pipelines.t1_linear.anat_linear_utils._call_n4_bias_field_correction",
            wraps=n4biasfieldcorrection_mock,
        ) as ants_bias_correction_mock:
            bias_corrected_image = run_n4biasfieldcorrection(
                tmp_path / "test.nii.gz",
                bspline_fitting_distance=300,
                output_prefix="sub-01_ses-M000",
                output_dir=output_dir,
                save_bias=True,
                verbose=True,
            )
            image_write_mock.assert_called()
            ants_bias_correction_mock.assert_called_with(
                tmp_path / "test.nii.gz",
                300,
                save_bias=True,
                verbose=True,
            )
    assert set([f.name for f in output_dir.iterdir()]) == {
        "sub-01_ses-M000_bias_corrected_image.nii.gz",
        "sub-01_ses-M000_bias_image.nii.gz",
    }
    assert bias_corrected_image.exists()
    bias_corrected_nifti = nib.load(bias_corrected_image)
    assert_array_equal(bias_corrected_nifti.affine, np.eye(4))
    assert_array_equal(bias_corrected_nifti.get_fdata(), data)


def generate_fake_fixed_and_moving_images(folder: Path):
    data = np.random.random((10, 10, 10))
    nib.save(nib.Nifti1Image(data, np.eye(4)), folder / "fixed.nii.gz")
    nib.save(nib.Nifti1Image(data, np.eye(4)), folder / "moving.nii.gz")


def test_run_ants_registration_error(tmp_path, mocker):
    import re

    from clinica.pipelines.t1_linear.anat_linear_utils import run_ants_registration

    generate_fake_fixed_and_moving_images(tmp_path)
    mocker.patch(
        "clinica.pipelines.t1_linear.anat_linear_utils._call_ants_registration",
        return_value={},
    )
    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "Something went wrong when calling antsRegistration with the following parameters :\n"
            f"- fixed_image = {tmp_path / 'fixed.nii.gz'}\n"
            f"- moving_image = {tmp_path / 'moving.nii.gz'}\n"
            f"- random_seed = 0\n"
            f"- type_of_transformation='antsRegistrationSyN[a]'\n"
        ),
    ):
        run_ants_registration(
            tmp_path / "fixed.nii.gz",
            tmp_path / "moving.nii.gz",
            random_seed=0,
        )


def ants_registration_mock(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    verbose: bool = False,
) -> dict:
    workdir = fixed_image.parent / "workdir"
    workdir.mkdir()
    mocked_transform = workdir / "transform.mat"
    mocked_transform.touch()
    return {
        "warpedmovout": nib.load(fixed_image),
        "fwdtransforms": ["fooo.txt", mocked_transform],
        "foo": "bar",
    }


def test_run_ants_registration(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_utils import run_ants_registration

    output_dir = tmp_path / "out"
    output_dir.mkdir()
    generate_fake_fixed_and_moving_images(tmp_path)

    with patch(
        "clinica.pipelines.t1_linear.anat_linear_utils._call_ants_registration",
        wraps=ants_registration_mock,
    ) as mock1:
        with patch("ants.image_write", wraps=nib.save) as mock2:
            run_ants_registration(
                tmp_path / "fixed.nii.gz",
                tmp_path / "moving.nii.gz",
                random_seed=12,
                output_dir=output_dir,
            )
            mock1.assert_called_once_with(
                tmp_path / "fixed.nii.gz", tmp_path / "moving.nii.gz", 12, verbose=False
            )
            mock2.assert_called_once()
