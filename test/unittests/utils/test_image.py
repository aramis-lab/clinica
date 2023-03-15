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
