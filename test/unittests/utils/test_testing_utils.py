from typing import Optional, Tuple

import nibabel as nib
import numpy as np
import pytest

TEST_IMAGE_DATA_SHAPE = (16, 10, 6, 3)


def _get_image(
    data: np.ndarray, affine_matrix: Optional[np.ndarray] = None
) -> nib.Nifti1Image:
    if affine_matrix is None:
        affine_matrix = np.eye(4)
    return nib.Nifti1Image(data, affine=affine_matrix)


def _get_random_image(
    data_shape: Tuple[int, ...], affine_matrix: Optional[np.ndarray] = None
) -> nib.Nifti1Image:
    if affine_matrix is None:
        affine_matrix = np.eye(4)
    return _get_image(np.random.random(data_shape), affine_matrix)


@pytest.fixture
def random_image():
    return _get_random_image(TEST_IMAGE_DATA_SHAPE)


def test_assert_nifti_equal(tmp_path, random_image):
    from clinica.utils.testing_utils import assert_nifti_equal

    nib.save(random_image, tmp_path / "image.nii.gz")
    nib.save(random_image, tmp_path / "image2.nii.gz")

    assert_nifti_equal(tmp_path / "image.nii.gz", tmp_path / "image.nii.gz")
    assert_nifti_equal(tmp_path / "image.nii.gz", tmp_path / "image2.nii.gz")


def test_assert_nifti_equal_different_dataobj_failure(tmp_path, random_image):
    from clinica.utils.testing_utils import assert_nifti_equal

    nib.save(random_image, tmp_path / "image.nii.gz")
    nib.save(_get_random_image(TEST_IMAGE_DATA_SHAPE), tmp_path / "image2.nii.gz")

    with pytest.raises(AssertionError):
        assert_nifti_equal(tmp_path / "image.nii.gz", tmp_path / "image2.nii.gz")


def test_assert_nifti_equal_different_affine_failure(tmp_path):
    from clinica.utils.testing_utils import assert_nifti_equal

    data = np.random.random((16, 10, 6, 3))
    nib.save(_get_image(data), tmp_path / "image.nii.gz")
    nib.save(_get_image(data, np.random.random((4, 4))), tmp_path / "image2.nii.gz")

    with pytest.raises(AssertionError):
        assert_nifti_equal(tmp_path / "image.nii.gz", tmp_path / "image2.nii.gz")


def test_assert_nifti_almost_equal(tmp_path):
    from clinica.utils.testing_utils import assert_nifti_almost_equal

    data = np.random.random((16, 10, 6, 3))
    noisy_data = data + 10**-7 * np.random.random(data.shape)
    nib.save(_get_image(data), tmp_path / "image.nii.gz")
    nib.save(_get_image(noisy_data), tmp_path / "image2.nii.gz")

    assert_nifti_almost_equal(tmp_path / "image.nii.gz", tmp_path / "image2.nii.gz")


def test_assert_nifti_almost_equal_fail(tmp_path):
    from clinica.utils.testing_utils import assert_nifti_almost_equal

    data = np.random.random((16, 10, 6, 3))
    noisy_data = data + 10**-4 * np.random.random(data.shape)
    nib.save(_get_image(data), tmp_path / "image.nii.gz")
    nib.save(_get_image(noisy_data), tmp_path / "image2.nii.gz")

    with pytest.raises(AssertionError):
        assert_nifti_almost_equal(tmp_path / "image.nii.gz", tmp_path / "image2.nii.gz")


def test_assert_large_image_dataobj_almost_equal():
    from clinica.utils.testing_utils import _assert_large_image_dataobj_almost_equal

    data = np.random.random((60, 60, 60, 20))
    noisy_data = data + 10**-7 * data
    img1 = _get_image(data)
    img2 = _get_image(noisy_data)

    _assert_large_image_dataobj_almost_equal(img1, img2)
    _assert_large_image_dataobj_almost_equal(img1, img2, n_samples=5)
    _assert_large_image_dataobj_almost_equal(img1, img2, n_samples=10, verbose=True)

    with pytest.raises(AssertionError):
        _assert_large_image_dataobj_almost_equal(img1, img2, decimal=9)
