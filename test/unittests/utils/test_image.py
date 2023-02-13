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


@pytest.mark.parametrize("axis", [-2, 5, 2.34])
def test_merge_volumes_errors(tmp_path, axis):
    from clinica.utils.image import merge_volumes

    with pytest.raises(
        ValueError,
        match="Axis should be an integer",
    ):
        merge_volumes(tmp_path / "vol1.nii.gz", tmp_path / "vol2.nii.gz", axis)


@pytest.mark.parametrize("axis", [0, 1, 2, 3, -1])
def test_merge_volumes(tmp_path, axis):
    from clinica.utils.image import merge_volumes

    input_data = []
    for i in (1, 2):
        img_data = np.zeros((5, 5, 5, 10))
        img_data[2:4, 2:4, 2:4, 0:2] = i
        input_data.append(img_data)
        img = nib.Nifti1Image(img_data, affine=np.eye(4))
        nib.save(img, tmp_path / f"foo{i}.nii.gz")
    out_file = merge_volumes(
        tmp_path / "foo1.nii.gz", tmp_path / "foo2.nii.gz", axis=axis
    )
    out_img = nib.load(out_file)
    assert out_file == tmp_path / "merged_files.nii.gz"
    assert_array_equal(out_img.affine, img.affine)
    assert_array_equal(out_img.get_fdata(), np.concatenate(input_data, axis=axis))
    (tmp_path / "bar").mkdir()
    nib.save(img, tmp_path / "bar" / f"foo2.nii.gz")
    with pytest.warns(
        UserWarning,
        match="Merging volumes",
    ):
        out_file = merge_volumes(
            tmp_path / "bar" / "foo2.nii.gz", tmp_path / "foo1.nii.gz", axis=axis
        )
    out_img = nib.load(out_file)
    assert out_file == tmp_path / "bar" / "merged_files.nii.gz"
    assert_array_equal(out_img.affine, img.affine)
    assert_array_equal(out_img.get_fdata(), np.concatenate(input_data[::-1], axis=axis))
    assert (
        merge_volumes(
            tmp_path / "foo1.nii.gz",
            tmp_path / "foo2.nii.gz",
            axis=axis,
            out_file=tmp_path / "foo3.nii.gz",
        )
        == tmp_path / "foo3.nii.gz"
    )
