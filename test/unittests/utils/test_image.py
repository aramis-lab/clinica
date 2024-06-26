import re

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal


@pytest.mark.parametrize(
    "aggregator,expected_value",
    [
        (np.average, 5.25),
        (np.mean, 5.25),
        (np.median, 5.5),
        (np.min, 1.0),
        (np.max, 9.0),
    ],
)
def test_compute_aggregated_volume(tmp_path, aggregator, expected_value):
    from clinica.utils.image import compute_aggregated_volume

    img_data = np.zeros((5, 5, 5, 8))
    img_data[2:4, 2:4, 2:4, 0:2] = 1.0
    img_data[2:4, 2:4, 2:4, 2:4] = 2.0
    img_data[2:4, 2:4, 2:4, 4:6] = 3.0
    img_data[2:4, 2:4, 2:4, 4:8] = 9.0
    img = nib.Nifti1Image(img_data, affine=np.eye(4))
    nib.save(img, tmp_path / "foo.nii.gz")
    expected = np.zeros((5, 5, 5))
    expected[2:4, 2:4, 2:4] = expected_value
    assert_array_equal(
        compute_aggregated_volume(tmp_path / "foo.nii.gz", aggregator), expected
    )


def test_get_new_image_like(tmp_path):
    from clinica.utils.image import get_new_image_like

    img_data = np.zeros((5, 5, 5, 8))
    img_data[2:4, 2:4, 2:4, 0:2] = 1.0
    img = nib.Nifti1Image(img_data, affine=np.eye(4))
    nib.save(img, tmp_path / "foo.nii.gz")
    img2 = get_new_image_like(tmp_path / "foo.nii.gz", np.ones((2, 6, 8, 10)))
    assert isinstance(img2, nib.nifti1.Nifti1Image)
    assert_array_equal(img.affine, img2.affine)


def test_merge_nifti_images_in_time_dimension_wrong_file_number(tmp_path):
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    with pytest.raises(
        ValueError,
        match="At least 2 files are required.",
    ):
        merge_nifti_images_in_time_dimension(())
        merge_nifti_images_in_time_dimension((tmp_path / "vol1.nii.gz"))


def test_merge_nifti_images_in_time_dimension_non_existing_file(tmp_path):
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    with pytest.raises(
        FileNotFoundError,
        match="the following",
    ):
        merge_nifti_images_in_time_dimension(
            (tmp_path / "vol1.nii.gz", tmp_path / "vol2.nii.gz")
        )


def test_merge_nifti_images_in_time_dimension_wrong_dimension(tmp_path):
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    for i, input_data in enumerate((np.zeros((5, 5, 5, 10)), np.zeros((5, 5)))):
        img = nib.Nifti1Image(input_data, affine=np.eye(4))
        nib.save(img, tmp_path / f"foo{i}.nii.gz")

    with pytest.raises(
        ValueError,
        match="Only 3D or 4D images can be concatenated.",
    ):
        merge_nifti_images_in_time_dimension(
            tuple(tmp_path / f"foo{i}.nii.gz" for i in range(2))
        )


@pytest.mark.parametrize("nb_images", [2, 3, 6])
def test_merge_nifti_images_in_time_dimension(tmp_path, nb_images):
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    input_data = []
    for i in range(nb_images):
        img_data = np.zeros((5, 5, 5, 10))
        img_data[2:4, 2:4, 2:4, 0:2] = i
        input_data.append(img_data)
        img = nib.Nifti1Image(img_data, affine=np.eye(4))
        nib.save(img, tmp_path / f"foo{i}.nii.gz")
    out_file = merge_nifti_images_in_time_dimension(
        tuple(tmp_path / f"foo{i}.nii.gz" for i in range(nb_images)),
        out_file=tmp_path / "merged_files.nii.gz",
    )
    out_img = nib.load(out_file)
    assert out_file == tmp_path / "merged_files.nii.gz"
    assert_array_equal(out_img.affine, np.eye(4))
    assert_array_equal(out_img.get_fdata(), np.concatenate(input_data, axis=-1))


@pytest.mark.parametrize(
    "shape1,shape2,expected_shape",
    [
        (
            (5, 5, 5, 10),
            (5, 5, 5),
            (5, 5, 5, 11),
        ),
        (
            (5, 5, 5),
            (5, 5, 5, 10),
            (5, 5, 5, 11),
        ),
        (
            (5, 5, 5),
            (5, 5, 5),
            (5, 5, 5, 2),
        ),
    ],
)
def test_merge_nifti_images_in_time_dimension_3d_and_4d(
    tmp_path, shape1, shape2, expected_shape
):
    """Checks that 3D and 4D images are merged correctly."""
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    input_data = (np.zeros(shape1), np.zeros(shape2))
    for i, data in enumerate(input_data):
        img = nib.Nifti1Image(data, affine=np.eye(4))
        nib.save(img, tmp_path / f"foo{i}.nii.gz")

    out_file = merge_nifti_images_in_time_dimension(
        tuple(tmp_path / f"foo{i}.nii.gz" for i in (0, 1)),
        out_file=tmp_path / "merged_files.nii.gz",
    )
    out_img = nib.load(out_file)
    assert out_file == tmp_path / "merged_files.nii.gz"
    assert_array_equal(out_img.affine, np.eye(4))
    expected_data = np.zeros(expected_shape)
    assert_array_equal(out_img.get_fdata(), expected_data)


def test_remove_dummy_dimension_from_image(tmp_path):
    from clinica.utils.image import remove_dummy_dimension_from_image

    input_data = np.random.randint(low=0, high=10, size=(16, 10, 6, 3))
    input_image = nib.Nifti1Image(input_data.astype(np.int32), affine=np.eye(4))
    nib.save(input_image, tmp_path / "input_image.nii.gz")

    result = remove_dummy_dimension_from_image(
        str(tmp_path / "input_image.nii.gz"),
        str(tmp_path / "output_image.nii.gz"),
    )
    result_image = nib.load(result)

    assert result == str(tmp_path / "output_image.nii.gz")
    assert_array_equal(result_image.affine, np.eye(4))
    assert_array_equal(result_image.get_fdata(), input_data)


def test_slice_error():
    from clinica.utils.image import Slice  # noqa

    with pytest.raises(
        ValueError,
        match=re.escape(
            "Slice instance has a start value (100) larger than the end value (10)."
        ),
    ):
        Slice(100, 10)


def test_slice():
    from clinica.utils.image import Slice  # noqa

    s = Slice(3, 16)

    assert s.start == 3
    assert s.end == 16
    assert s.get_slice() == slice(3, 16)


def test_mni_cropped_bbox():
    from clinica.utils.image import MNI_CROP_BBOX  # noqa

    assert MNI_CROP_BBOX.x_slice.start == 12
    assert MNI_CROP_BBOX.x_slice.end == 181
    assert MNI_CROP_BBOX.y_slice.start == 13
    assert MNI_CROP_BBOX.y_slice.end == 221
    assert MNI_CROP_BBOX.z_slice.start == 0
    assert MNI_CROP_BBOX.z_slice.end == 179


def test_load_mni_cropped_template(tmp_path, mocker):
    from clinica.utils.image import _load_mni_cropped_template  # noqa

    expected_affine = np.array(
        [
            [1.0, 0.0, 0.0, -84.0],
            [0.0, 1.0, 0.0, -119.0],
            [0.0, 0.0, 1.0, -78.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    expected_shape = (169, 208, 179)
    mocked = nib.Nifti1Image(np.zeros(expected_shape), affine=expected_affine)
    mocked.to_filename(tmp_path / "mocked.nii.gz")
    mocker.patch(
        "clinica.utils.image._get_file_locally_or_download",
        return_value=tmp_path / "mocked.nii.gz",
    )

    img = _load_mni_cropped_template()

    assert_array_equal(img.affine, expected_affine)
    assert img.shape == expected_shape


def test_load_mni_template(tmp_path, mocker):
    from clinica.utils.image import _load_mni_template  # noqa

    expected_affine = np.array(
        [
            [1.0, 0.0, 0.0, -96.0],
            [0.0, 1.0, 0.0, -132.0],
            [0.0, 0.0, 1.0, -78.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    )
    expected_shape = (193, 229, 193)
    mocked = nib.Nifti1Image(np.zeros(expected_shape), affine=expected_affine)
    mocked.to_filename(tmp_path / "mocked.nii.gz")
    mocker.patch(
        "clinica.utils.image._get_file_locally_or_download",
        return_value=tmp_path / "mocked.nii.gz",
    )

    img = _load_mni_template()

    assert_array_equal(img.affine, expected_affine)
    assert img.shape == expected_shape


def test_crop_nifti_error(tmp_path):
    from clinica.utils.image import crop_nifti

    nib.Nifti1Image(np.random.random((10, 10, 10, 10)), np.eye(4)).to_filename(
        tmp_path / "test.nii.gz"
    )

    with pytest.raises(
        ValueError,
        match=re.escape(
            "The function crop_nifti is implemented for anatomical 3D images. "
            "You provided an image of shape (10, 10, 10, 10)."
        ),
    ):
        crop_nifti(tmp_path / "test.nii.gz")


def test_crop_nifti(tmp_path):
    from clinica.utils.image import (
        _load_mni_cropped_template,  # noqa
        _load_mni_template,  # noqa
        crop_nifti,
    )

    _load_mni_template().to_filename(tmp_path / "mni.nii.gz")
    crop_nifti(tmp_path / "mni.nii.gz", output_dir=tmp_path)

    assert (tmp_path / "mni_cropped.nii.gz").exists()
    cropped = nib.load(tmp_path / "mni_cropped.nii.gz")
    assert_array_equal(cropped.affine, _load_mni_cropped_template().affine)
    assert_array_equal(cropped.get_fdata(), _load_mni_cropped_template().get_fdata())
