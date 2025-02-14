import re

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from clinica.utils.exceptions import ClinicaImageDimensionError
from clinica.utils.testing_utils import (
    assert_nifti_equal,
    build_test_image_cubic_object,
)


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


@pytest.mark.parametrize(
    "array,coords,expected",
    [
        (np.ones((10, 10, 10)), (3, 6, 0, 10, 6, 6), True),
        (np.ones((10, 10, 10)), (3, 11, 0, 10, 6, 6), False),
        (np.ones((10, 10, 10)), (3, 6, -1, 10, 6, 6), False),
        (np.ones((10, 10, 10)), (1, 1, 1, 1, 1, 1), True),
        (np.ones((10, 10, 10)), (11, 11, 11, 11, 11, 11), False),
    ],
)
def test_is_bbox_within_array(
    array: np.ndarray, coords: tuple[int, int, int, int, int, int], expected: bool
):
    from clinica.utils.image import Bbox3D, _is_bbox_within_array  # noqa

    assert _is_bbox_within_array(array, Bbox3D.from_coordinates(*coords)) is expected


def test_mni_cropped_bbox():
    from clinica.utils.image import MNI_CROP_BBOX  # noqa

    assert MNI_CROP_BBOX.x_slice.start == 12
    assert MNI_CROP_BBOX.x_slice.end == 181
    assert MNI_CROP_BBOX.y_slice.start == 13
    assert MNI_CROP_BBOX.y_slice.end == 221
    assert MNI_CROP_BBOX.z_slice.start == 0
    assert MNI_CROP_BBOX.z_slice.end == 179


def test_get_mni_cropped_template(tmp_path, mocker):
    from clinica.utils.image import get_mni_cropped_template

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

    img = nib.load(get_mni_cropped_template())

    assert_array_equal(img.affine, expected_affine)
    assert img.shape == expected_shape


def test_get_mni_template(tmp_path, mocker):
    from clinica.utils.image import get_mni_template

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

    img = nib.load(get_mni_template("t1"))

    assert_array_equal(img.affine, expected_affine)
    assert img.shape == expected_shape


def test_crop_nifti_input_image_not_3d_error(tmp_path):
    from clinica.utils.image import crop_nifti

    nib.Nifti1Image(np.random.random((10, 10, 10, 10)), np.eye(4)).to_filename(
        tmp_path / "test.nii.gz"
    )

    with pytest.raises(
        ClinicaImageDimensionError,
        match=f"The image in {tmp_path / 'test.nii.gz'} is not 3D.",
    ):
        crop_nifti(tmp_path / "test.nii.gz")


def test_crop_nifti_no_cropping(tmp_path):
    from clinica.utils.image import crop_nifti

    input_image = nib.Nifti1Image(np.random.random((10, 10, 10)), np.eye(4))
    input_image.to_filename(tmp_path / "test.nii.gz")

    cropped = crop_nifti(tmp_path / "test.nii.gz", output_dir=tmp_path)

    assert_nifti_equal(tmp_path / "test.nii.gz", tmp_path / "test_cropped.nii.gz")


def test_crop_nifti_using_t1_mni_template_with_resampling(tmp_path):
    from clinica.utils.image import (
        crop_nifti_using_t1_mni_template,
        get_mni_cropped_template,
        get_mni_template,
    )

    nib.load(get_mni_template("flair")).to_filename(tmp_path / "mni.nii.gz")
    with pytest.warns(
        UserWarning,
        match=re.escape(
            f"The image {tmp_path / 'mni.nii.gz'} has dimensions (182, 218, 182) and cannot "
            "be cropped using the bounding box ( ( 12, 181 ), ( 13, 221 ), ( 0, 179 ) ). "
            "The `crop_nifti` function will try to resample the input image to the reference template"
        ),
    ):
        cropped = crop_nifti_using_t1_mni_template(
            tmp_path / "mni.nii.gz", output_dir=tmp_path
        )
        assert nib.load(cropped).shape == nib.load(get_mni_cropped_template()).shape


def test_crop_nifti_using_t1_mni_template(tmp_path):
    from clinica.utils.image import (
        crop_nifti_using_t1_mni_template,
        get_mni_cropped_template,
        get_mni_template,
    )

    nib.load(get_mni_template("t1")).to_filename(tmp_path / "mni.nii.gz")
    crop_nifti_using_t1_mni_template(tmp_path / "mni.nii.gz", output_dir=tmp_path)

    assert (tmp_path / "mni_cropped.nii.gz").exists()
    cropped = nib.load(tmp_path / "mni_cropped.nii.gz")
    assert_array_equal(cropped.affine, nib.load(get_mni_cropped_template()).affine)
    assert_array_equal(
        cropped.get_fdata(), nib.load(get_mni_cropped_template()).get_fdata()
    )


@pytest.mark.parametrize(
    "input_background,input_object,threshold_low,threshold_high,expected_background,expected_object",
    [
        (-0.01, 1.0, 0.0, None, 0.0, 1.0),
        (1.1, 1.0, 1.1, None, 1.1, 1.1),
        (0.0, 0.0, 1.0, None, 1.0, 1.0),
        (0.0, 0.0, 0.0, None, 0.0, 0.0),
        (0.0, 0.0, -0.1, None, 0.0, 0.0),
        (10.0, 1.0, 5.0, None, 10.0, 5.0),
        (-0.01, 1.0, None, 0.0, -0.01, 0.0),
        (-0.01, 1.0, None, 1.0, -0.01, 1.0),
        (10.0, 1.0, None, 2.0, 2.0, 1.0),
        (10.0, 1.0, 2.0, 4.0, 4.0, 2.0),
        (0.0, 0.0, 2.0, 4.0, 2.0, 2.0),
        (10.0, 1.0, 1.0, 10.0, 10.0, 1.0),
    ],
)
def test_clip_nifti(
    tmp_path,
    input_background: float,
    input_object: float,
    threshold_low: float,
    threshold_high: float,
    expected_background: float,
    expected_object: float,
):
    from clinica.utils.image import clip_nifti

    build_test_image_cubic_object(
        shape=(10, 10, 10),
        background_value=input_background,
        object_value=input_object,
        object_size=2,
    ).to_filename(tmp_path / "input.nii.gz")

    clipped = clip_nifti(
        tmp_path / "input.nii.gz", low=threshold_low, high=threshold_high
    )

    build_test_image_cubic_object(
        shape=(10, 10, 10),
        background_value=expected_background,
        object_value=expected_object,
        object_size=2,
    ).to_filename(tmp_path / "expected.nii.gz")

    assert_nifti_equal(clipped, tmp_path / "expected.nii.gz")


def test_nifti_image_file_not_exist_error(tmp_path):
    from clinica.utils.image import NiftiImage

    with pytest.raises(
        FileNotFoundError,
        match=f"File {tmp_path / 'image.nii.gz'} does not exist.",
    ):
        NiftiImage(tmp_path / "image.nii.gz")


def test_nifti_image_file_corrupted_error(tmp_path):
    from clinica.utils.image import NiftiImage

    (tmp_path / "image.nii.gz").touch()

    with pytest.raises(
        IOError,
        match=f"File {tmp_path / 'image.nii.gz'} is not a nifti image or is corrupted.",
    ):
        NiftiImage(tmp_path / "image.nii.gz")


def test_nifti_image_get_filename(tmp_path):
    from clinica.utils.image import NiftiImage

    image = nib.Nifti1Image(np.random.random((10, 10, 10)), np.eye(4))
    nib.save(image, tmp_path / "image.nii.gz")
    image_instance = NiftiImage(tmp_path / "image.nii.gz")

    assert image_instance.get_filename() == "image.nii.gz"
    assert image_instance.get_filename(with_extension=False) == "image"


def test_nifti_image_3d_error(tmp_path):
    from clinica.utils.exceptions import ClinicaImageDimensionError
    from clinica.utils.image import NiftiImage3D

    image = nib.Nifti1Image(np.random.random((10, 10, 10, 10)), np.eye(4))
    nib.save(image, tmp_path / "image.nii.gz")

    with pytest.raises(
        ClinicaImageDimensionError,
        match=f"The image in {tmp_path / 'image.nii.gz'} is not 3D.",
    ):
        NiftiImage3D(tmp_path / "image.nii.gz")
