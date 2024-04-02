import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal


def test_extract_sub_ses_folder_name():
    from clinica.pipelines.dwi.preprocessing.t1.utils import extract_sub_ses_folder_name

    assert (
        extract_sub_ses_folder_name(
            "/localdrive10TB/wd/dwi-preprocessing-using-t1/epi_pipeline/"
            "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/MergeDWIs/"
            "Jacobian_image_maths_thresh_merged.nii.gz"
        )
        == "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea"
    )


def test_change_itk_transform_type_error(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import change_itk_transform_type

    with pytest.raises(FileNotFoundError):
        change_itk_transform_type(tmp_path / "affine.txt")


def test_change_itk_transform_type_empty_file(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import change_itk_transform_type

    (tmp_path / "affine.txt").touch()
    assert (
        change_itk_transform_type(tmp_path / "affine.txt")
        == tmp_path / "updated_affine.txt"
    )
    assert (tmp_path / "updated_affine.txt").read_text() == ""


def test_change_itk_transform_type_no_transform_line(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import change_itk_transform_type

    origin_affine = tmp_path / "affine.txt"
    origin_content = (
        "fooo\nbar \n bazzzz\n transform\nTransform: AffineTransform_double_3_3"
    )
    origin_affine.write_text(origin_content)

    assert change_itk_transform_type(origin_affine) == tmp_path / "updated_affine.txt"
    assert (tmp_path / "updated_affine.txt").read_text() == origin_content


def test_change_itk_transform_type(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import change_itk_transform_type

    origin_affine = tmp_path / "affine.txt"
    origin_content = "fooo\nbar \n bazzzz\n transform\nTransform: MatrixOffsetTransformBase_double_3_3"
    expected_content = (
        "fooo\nbar \n bazzzz\n transform\nTransform: AffineTransform_double_3_3"
    )
    origin_affine.write_text(origin_content)

    assert change_itk_transform_type(origin_affine) == tmp_path / "updated_affine.txt"
    assert (tmp_path / "updated_affine.txt").read_text() == expected_content


def test_broadcast_matrix_filename_to_match_b_vector_length_file_not_found_error(
    tmp_path,
):
    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        broadcast_matrix_filename_to_match_b_vector_length,
    )

    with pytest.raises(FileNotFoundError):
        broadcast_matrix_filename_to_match_b_vector_length(
            Path("matrix.mat"), tmp_path / "foo.txt"
        )


def test_broadcast_matrix_filename_to_match_b_vector_length_file_format_error(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        broadcast_matrix_filename_to_match_b_vector_length,
    )

    filename = tmp_path / "vectors.bvec"
    filename.write_text("foooo\nbarrrr\nbazzzz")

    with pytest.raises(
        ValueError,
        match="could not convert string 'foooo' to float64",
    ):
        broadcast_matrix_filename_to_match_b_vector_length(Path("matrix.mat"), filename)


def test_broadcast_matrix_filename_to_match_b_vector_length(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        broadcast_matrix_filename_to_match_b_vector_length,
    )

    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(6, 10))

    assert (
        broadcast_matrix_filename_to_match_b_vector_length(
            Path("matrix.mat"), tmp_path / "vectors.bvec"
        )
        == [Path("matrix.mat")] * 10
    )


def test_rotate_b_vectors_size_error(tmp_path):
    """Test that a RuntimeError is raised when the number of b-vectors
    is different from the number of matrices in the list.
    """
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(6, 10))

    with pytest.raises(RuntimeError, match="Number of b-vectors"):
        rotate_b_vectors(
            tmp_path / "vectors.bvec",
            ["mat1", "mat2", "mat3"],
            output_dir=tmp_path,
        )


def test_rotate_b_vectors_missing_b_vectors_file(tmp_path):
    """Test that a FileNotFoundError is raised when the b-vectors file is missing."""
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    with pytest.raises(FileNotFoundError):
        rotate_b_vectors(tmp_path / "foo.bvec", ["mat1", "mat2"], output_dir=tmp_path)


def test_rotate_b_vectors_missing_rotation_matrix_file(tmp_path):
    """Test that a FileNotFoundError is raised when one of the matrix file is missing."""
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(3, 2))
    np.savetxt(tmp_path / "mat1", np.eye(4))

    with pytest.raises(
        FileNotFoundError,
        match="mat2 not found",
    ):
        rotate_b_vectors(
            tmp_path / "vectors.bvec",
            [str(tmp_path / "mat1"), str(tmp_path / "mat2")],
            output_dir=tmp_path,
        )


def test_rotate_b_vectors_all_zeros(tmp_path):
    """Test that B-vectors composed only of zeros aren't rotated."""
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    np.savetxt(tmp_path / "vectors.bvec", np.zeros((3, 2)))
    for i in range(1, 3):
        np.savetxt(tmp_path / f"mat{i}", np.random.rand(3, 3))

    rotated_b_vectors_file = rotate_b_vectors(
        tmp_path / "vectors.bvec",
        [str(tmp_path / "mat1"), str(tmp_path / "mat2")],
        output_dir=tmp_path,
    )

    assert rotated_b_vectors_file == tmp_path / "vectors_rotated.bvec"
    assert_array_equal(np.loadtxt(rotated_b_vectors_file), np.zeros((3, 2)))


def test_rotate_b_vectors_with_identity(tmp_path):
    """Test that b-vectors rotated with identity matrix are identical."""
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    origin_b_vectors = np.random.rand(3, 2)
    normalized_origin_b_vectors = np.array(
        [b / np.linalg.norm(b) for b in origin_b_vectors.T]
    ).T
    np.savetxt(tmp_path / "vectors.bvec", normalized_origin_b_vectors)
    for i in range(1, 3):
        np.savetxt(tmp_path / f"mat{i}", np.eye(4))

    rotated_b_vectors_file = rotate_b_vectors(
        tmp_path / "vectors.bvec",
        [str(tmp_path / "mat1"), str(tmp_path / "mat2")],
        output_dir=tmp_path,
    )

    assert rotated_b_vectors_file == tmp_path / "vectors_rotated.bvec"
    assert_array_almost_equal(
        np.loadtxt(rotated_b_vectors_file), normalized_origin_b_vectors
    )


def test_rotate_b_vectors_wrong_content_in_vector_file(tmp_path):
    """Test that a ValueError is raised when vector file isn't valid."""
    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    b_vectors_file = tmp_path / "foo.bvec"
    b_vectors_file.write_text("fooo\nbar")

    with pytest.raises(ValueError, match="could not convert string 'fooo' to float64"):
        rotate_b_vectors(b_vectors_file, ["mat1"], output_dir=tmp_path)


def test_configure_working_directory(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import configure_working_directory

    dwi_file = tmp_path / "sub-foo_ses-bar_dwi.nii.gz"
    file_hash = hashlib.md5(str(dwi_file).encode()).hexdigest()
    dwi_file.touch()
    work_dir = configure_working_directory(dwi_file, tmp_path)
    assert work_dir.exists()
    assert work_dir == tmp_path / file_hash
    work_dir = configure_working_directory(dwi_file)
    assert work_dir.exists()
    assert work_dir.stem == file_hash


def test_compute_reference_b0_errors(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import compute_reference_b0

    (tmp_path / "b0.nii.gz").touch()
    np.savetxt(tmp_path / "foo.bval", [1000] * 9)
    with pytest.raises(
        ValueError,
        match="The number of b0s should be strictly positive",
    ):
        compute_reference_b0(tmp_path / "b0.nii.gz", tmp_path / "foo.bval", 5, tmp_path)


def test_compute_reference_b0_with_single_b0(tmp_path):
    """Test the function compute_reference_b0 when there is a single B0 volume.

    In this case, compute_reference_b0 simply return the provided image.
    """
    from clinica.pipelines.dwi.preprocessing.t1.utils import compute_reference_b0

    b0_data = 4.0 * np.ones((5, 5, 5, 1))
    b0_img = nib.Nifti1Image(b0_data, affine=np.eye(4))
    nib.save(b0_img, tmp_path / "b0.nii.gz")
    np.savetxt(tmp_path / "foo.bval", [[1000] * 5 + [0] + [1000] * 3])
    working_directory = tmp_path / "tmp"
    working_directory.mkdir()

    ref = compute_reference_b0(
        tmp_path / "b0.nii.gz", tmp_path / "foo.bval", 5, tmp_path
    )
    assert ref == tmp_path / "b0.nii.gz"
    ref_img = nib.load(ref)
    assert_array_equal(b0_data, ref_img.get_fdata())
    assert_array_equal(np.eye(4), ref_img.affine)


@pytest.mark.parametrize("clean_working_dir", [False, True])
def test_compute_reference_b0_with_multiple_b0(tmp_path, mocker, clean_working_dir):
    """Tests the function compute_reference_b0.

    `compute_reference_b0` calls `register_b0` which is a simple wrapper around the
    FSL pipeline implemented in `b0_flirt_pipeline`.
    This pipeline is tested as part of the non regression tests suite which runs
    on a much longer timescale and on dedicated machines having both FSL installed
    and access to proper testing data.
    This unit test simply checks the additional logic brought by `compute_reference_b0`,
    and should not depend on the installation of a third party software or on specific
    testing data. For these reasons we mock the `register_b0 function`.
    """
    from clinica.pipelines.dwi.preprocessing.t1.utils import compute_reference_b0

    working_directory = tmp_path / "tmp"
    working_directory.mkdir()
    mocked_output = (
        working_directory
        / "b0_coregistration"
        / "concat_ref_moving"
        / "merged_files.nii.gz"
    )
    assert not mocked_output.exists()
    mocked_output.parent.mkdir(parents=True)
    mocker.patch(
        "clinica.pipelines.dwi.preprocessing.t1.utils.register_b0",
        return_value=mocked_output,
    )
    # Build a dataset with 3 B0 volumes
    b0_data = 4.0 * np.ones((5, 5, 5, 3))
    b0_img = nib.Nifti1Image(b0_data, affine=np.eye(4))
    nib.save(b0_img, tmp_path / "b0.nii.gz")
    np.savetxt(tmp_path / "foo.bval", [[0] + [1000] * 3 + [0, 0] + [1000] * 3])

    # Build the data resulting from the call to register_b0
    # This is part of the mocking strategy
    co_registered_b0_data = 5.0 * np.ones((5, 5, 5, 3))
    co_registered_b0_img = nib.Nifti1Image(co_registered_b0_data, affine=np.eye(4))
    nib.save(co_registered_b0_img, mocked_output)
    assert mocked_output.exists()
    ref = compute_reference_b0(
        tmp_path / "b0.nii.gz",
        tmp_path / "foo.bval",
        b_value_threshold=5.0,
        working_directory=tmp_path / "tmp",
        clean_working_dir=clean_working_dir,
    )
    assert ref == tmp_path / "reference_b0_volume.nii.gz"
    if clean_working_dir:
        assert not mocked_output.exists()
        assert not working_directory.exists()
        assert sorted([p.name for p in tmp_path.iterdir()]) == [
            "b0.nii.gz",
            "foo.bval",
            "reference_b0_volume.nii.gz",
        ]
    else:
        assert mocked_output.exists()
        assert working_directory.exists()
        assert sorted([p.name for p in tmp_path.iterdir()]) == [
            "b0.nii.gz",
            "foo.bval",
            "reference_b0_volume.nii.gz",
            "tmp",
        ]
    ref_img = nib.load(ref)
    assert ref_img.shape == (5, 5, 5, 1)
    assert_array_equal(ref_img.get_fdata(), 5.0 * np.ones((5, 5, 5, 1)))


def test_prepare_reference_b0(tmp_path, mocker):
    from clinica.pipelines.dwi.preprocessing.t1.utils import prepare_reference_b0
    from clinica.pipelines.dwi.utils import DWIDataset

    working_dir = tmp_path / "working_dir"
    working_dir.mkdir()
    mocked_output = working_dir / "reference_b0_volume.nii.gz"
    mocker.patch(
        "clinica.pipelines.dwi.preprocessing.t1.utils.compute_reference_b0",
        return_value=mocked_output,
    )
    bids_folder = tmp_path / "bids"
    bids_folder.mkdir()
    ref_b0_data = 5.0 * np.ones((5, 5, 5, 1))
    ref_b0_img = nib.Nifti1Image(ref_b0_data, affine=np.eye(4))
    nib.save(ref_b0_img, mocked_output)

    dwi_data = 4.0 * np.ones((5, 5, 5, 9))
    dwi_img = nib.Nifti1Image(dwi_data, affine=np.eye(4))
    nib.save(dwi_img, bids_folder / "sub-foo_ses-bar_dwi.nii.gz")
    np.savetxt(
        bids_folder / "sub-foo_ses-bar_dwi.bval",
        [[0] + [1000] * 3 + [0, 0] + [1000] * 3],
    )
    bvecs_data = np.random.random((3, 9))
    np.savetxt(bids_folder / "sub-foo_ses-bar_dwi.bvec", bvecs_data)

    dataset = DWIDataset(
        dwi=bids_folder / "sub-foo_ses-bar_dwi.nii.gz",
        b_values=bids_folder / "sub-foo_ses-bar_dwi.bval",
        b_vectors=bids_folder / "sub-foo_ses-bar_dwi.bvec",
    )
    reference_b0, reference_dataset = prepare_reference_b0(
        dataset,
        b_value_threshold=5.0,
        working_directory=working_dir,
    )
    dwi_filename_hash = hashlib.md5(str(dataset.dwi).encode()).hexdigest()
    assert reference_b0 == working_dir / "reference_b0_volume.nii.gz"
    assert (
        reference_dataset.dwi
        == working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.nii.gz"
    )
    assert (
        reference_dataset.b_values
        == working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.bval"
    )
    assert (
        reference_dataset.b_vectors
        == working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.bvec"
    )
    assert sorted([p.name for p in tmp_path.iterdir()]) == ["bids", "working_dir"]
    assert sorted([p.name for p in bids_folder.iterdir()]) == [
        "sub-foo_ses-bar_dwi.bval",
        "sub-foo_ses-bar_dwi.bvec",
        "sub-foo_ses-bar_dwi.nii.gz",
    ]
    assert sorted([p.name for p in working_dir.iterdir()]) == [
        dwi_filename_hash,
        "reference_b0_volume.nii.gz",
    ]
    assert (
        sorted([p.name for p in (working_dir / dwi_filename_hash).iterdir()])
        == [
            "sub-foo_ses-bar_dwi_large_b.bval",  # bvalue file corresponding to DWI volumes with b>low_b
            "sub-foo_ses-bar_dwi_large_b.bvec",  # bvectors file corresponding to DWI volumes with b>low_b
            "sub-foo_ses-bar_dwi_large_b.nii.gz",  # DWI image file holding volumes for which b>low_b
            "sub-foo_ses-bar_dwi_merged.bval",  # bvalue file corresponding to the merged DWI dataset
            "sub-foo_ses-bar_dwi_merged.bvec",  # bvectors file corresponding to the merged DWI dataset
            "sub-foo_ses-bar_dwi_merged.nii.gz",  # image file holding the merged DWI volumes
            "sub-foo_ses-bar_dwi_small_b.bval",  # bvalue file corresponding to the volumes for which b<=low_b
            "sub-foo_ses-bar_dwi_small_b.bvec",  # bvectors file corresponding to the volumes for which b<=low_b
            "sub-foo_ses-bar_dwi_small_b.nii.gz",  # DWI image file holding volumes for which b<=low_b
        ]
    )
    ref_b0_volume = nib.load(working_dir / "reference_b0_volume.nii.gz")
    assert_array_equal(ref_b0_volume.get_fdata(), 5.0 * np.ones((5, 5, 5, 1)))
    large_b_image = nib.load(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_large_b.nii.gz"
    )
    assert_array_equal(large_b_image.get_fdata(), 4.0 * np.ones((5, 5, 5, 6)))
    large_b_values = np.loadtxt(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_large_b.bval"
    )
    assert_array_equal(large_b_values, np.array([1000] * 6))
    large_b_vectors = np.loadtxt(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_large_b.bvec"
    )
    assert large_b_vectors.shape == (3, 6)
    small_b_image = nib.load(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_small_b.nii.gz"
    )
    assert_array_equal(small_b_image.get_fdata(), 4.0 * np.ones((5, 5, 5, 3)))
    merged_image = nib.load(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.nii.gz"
    )
    expected_merged_image_data = np.concatenate(
        (
            5.0 * np.ones((5, 5, 5, 1)),
            4.0 * np.ones((5, 5, 5, 6)),
        ),
        axis=-1,
    )
    assert_array_equal(merged_image.get_fdata(), expected_merged_image_data)
    merged_b_values = np.loadtxt(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.bval"
    )
    assert_array_equal(merged_b_values, np.array([0] + [1000] * 6))
    merged_b_vectors = np.loadtxt(
        working_dir / dwi_filename_hash / "sub-foo_ses-bar_dwi_merged.bvec"
    )
    assert merged_b_vectors.shape == (3, 7)


@pytest.mark.parametrize(
    "threshold,expected", [(None, 4), (-1, 0), (500, 4), (1000, 8), (1001, 8)]
)
def test_count_b0s(tmp_path, threshold, expected):
    from clinica.pipelines.dwi.preprocessing.t1.utils import _count_b0s  # noqa

    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    kwargs = {"b_value_threshold": threshold} if threshold else {}
    assert _count_b0s(tmp_path / "foo.bval", **kwargs) == expected


@pytest.fixture
def dwi_dataset(tmp_path):
    from clinica.pipelines.dwi.utils import DWIDataset

    return DWIDataset(
        dwi=str(tmp_path / "foo.nii.gz"),
        b_values=tmp_path / "foo.bval",
        b_vectors=tmp_path / "foo.bvec",
    )


def test_split_dwi_dataset_with_b_values_errors(tmp_path, dwi_dataset):
    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        split_dwi_dataset_with_b_values,
    )

    for filename in ("foo.nii.gz", "foo.bval", "foo.bvec"):
        (tmp_path / filename).touch()

    with pytest.raises(
        ValueError,
        match="b_value_threshold should be >=0. You provided -1.",
    ):
        split_dwi_dataset_with_b_values(dwi_dataset, b_value_threshold=-1)


@pytest.mark.parametrize("extension", ["nii", "nii.gz"])
def test_split_dwi_dataset_with_b_values(tmp_path, extension):
    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        split_dwi_dataset_with_b_values,
    )
    from clinica.pipelines.dwi.utils import DWIDataset

    b_vectors_data = np.random.random((3, 8))
    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    np.savetxt(tmp_path / "foo.bvec", b_vectors_data)
    img_data = np.zeros((5, 5, 5, 8))
    img_data[2:4, 2:4, 2:4, 0:4] = 1.0
    img_data[2:4, 2:4, 2:4, 4:8] = 2.0
    img = nib.Nifti1Image(img_data, affine=np.eye(4))
    nib.save(img, tmp_path / f"foo.{extension}")

    dwi_dataset = DWIDataset(
        dwi=tmp_path / f"foo.{extension}",
        b_values=tmp_path / "foo.bval",
        b_vectors=tmp_path / "foo.bvec",
    )
    small_b_dataset, large_b_dataset = split_dwi_dataset_with_b_values(dwi_dataset)
    assert small_b_dataset.dwi == tmp_path / f"foo_small_b.{extension}"
    assert small_b_dataset.b_values == tmp_path / "foo_small_b.bval"
    assert small_b_dataset.b_vectors == tmp_path / "foo_small_b.bvec"
    assert large_b_dataset.dwi == tmp_path / f"foo_large_b.{extension}"
    assert large_b_dataset.b_values == tmp_path / "foo_large_b.bval"
    assert large_b_dataset.b_vectors == tmp_path / "foo_large_b.bvec"
    b0_img = nib.load(small_b_dataset.dwi)
    expected = np.zeros((5, 5, 5, 4))
    expected[2:4, 2:4, 2:4, 0:2] = 1.0
    expected[2:4, 2:4, 2:4, 2:4] = 2.0
    assert_array_equal(b0_img.get_fdata(), expected)
    dwi = nib.load(large_b_dataset.dwi)
    assert_array_equal(dwi.get_fdata(), expected)
    b_values = np.loadtxt(large_b_dataset.b_values)
    assert_array_equal(b_values, np.array([1000] * 4))
    b_vectors = np.loadtxt(large_b_dataset.b_vectors)
    assert_array_almost_equal(
        b_vectors, b_vectors_data[:, np.array([0, 1, 5, 6])], decimal=5
    )


def test_insert_b0_into_dwi(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.utils import insert_b0_into_dwi
    from clinica.utils.testing_utils import build_dwi_dataset

    b0_data = 6.0 * np.ones((5, 5, 5, 1))
    b0_img = nib.Nifti1Image(b0_data, affine=np.eye(4))
    nib.save(b0_img, tmp_path / "b0.nii.gz")
    dwi_dataset = build_dwi_dataset(tmp_path, 9, 9, 9)

    out_dataset = insert_b0_into_dwi(tmp_path / "b0.nii.gz", dwi_dataset)
    assert out_dataset.dwi == tmp_path / "foo_merged.nii.gz"
    assert out_dataset.b_values == tmp_path / "foo_merged.bval"
    assert out_dataset.b_vectors == tmp_path / "foo_merged.bvec"
    dwi = nib.load(out_dataset.dwi)
    dwi_img = nib.load(tmp_path / "foo.nii.gz")
    assert_array_equal(dwi.affine, dwi_img.affine)
    expected = 4.0 * np.ones((5, 5, 5, 10))
    expected[..., 0] += 2.0
    assert_array_equal(dwi.get_fdata(), expected)
    bvals = np.loadtxt(out_dataset.b_values)
    assert_array_equal(bvals, np.array([0] + [1000] * 9))
    bvecs = np.loadtxt(out_dataset.b_vectors)
    bvecs_data = np.loadtxt(tmp_path / "foo.bvec")
    expected = np.insert(bvecs_data, 0, 0.0, axis=1)
    assert_array_almost_equal(bvecs, expected, decimal=5)
