import os
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import (
    broadcast_matrix_filename_to_match_b_vector_length,
    change_itk_transform_type,
    delete_temp_dirs,
    extract_sub_ses_folder_name,
    rotate_b_vectors,
)


def test_extract_sub_ses_folder_name():
    assert (
        extract_sub_ses_folder_name(
            "/localdrive10TB/wd/dwi-preprocessing-using-t1/epi_pipeline/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/MergeDWIs/Jacobian_image_maths_thresh_merged.nii.gz"
        )
        == "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea"
    )


def test_change_itk_transform_type_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        change_itk_transform_type(str(tmp_path / "affine.txt"))


def test_change_itk_transform_type_empty_file(tmp_path):
    (tmp_path / "affine.txt").touch()

    assert change_itk_transform_type(str(tmp_path / "affine.txt")) == str(
        tmp_path / "updated_affine.txt"
    )
    assert (tmp_path / "updated_affine.txt").read_text() == ""


def test_change_itk_transform_type_no_transform_line(tmp_path):
    origin_affine = tmp_path / "affine.txt"
    origin_content = (
        "fooo\nbar \n bazzzz\n transform\nTransform: AffineTransform_double_3_3"
    )
    origin_affine.write_text(origin_content)

    assert change_itk_transform_type(str(origin_affine)) == str(
        tmp_path / "updated_affine.txt"
    )
    assert (tmp_path / "updated_affine.txt").read_text() == origin_content


def test_change_itk_transform_type(tmp_path):
    origin_affine = tmp_path / "affine.txt"
    origin_content = "fooo\nbar \n bazzzz\n transform\nTransform: MatrixOffsetTransformBase_double_3_3"
    expected_content = (
        "fooo\nbar \n bazzzz\n transform\nTransform: AffineTransform_double_3_3"
    )
    origin_affine.write_text(origin_content)

    assert change_itk_transform_type(str(origin_affine)) == str(
        tmp_path / "updated_affine.txt"
    )
    assert (tmp_path / "updated_affine.txt").read_text() == expected_content


def test_broadcast_matrix_filename_to_match_b_vector_length_file_not_found_error(
    tmp_path,
):
    with pytest.raises(FileNotFoundError):
        broadcast_matrix_filename_to_match_b_vector_length(
            "matrix.mat", str(tmp_path / "foo.txt")
        )


def test_broadcast_matrix_filename_to_match_b_vector_length_file_format_error(tmp_path):
    filename = tmp_path / "vectors.bvec"
    filename.write_text("foooo\nbarrrr\nbazzzz")

    with pytest.raises(
        ValueError,
        match="could not convert string 'foooo' to float64",
    ):
        broadcast_matrix_filename_to_match_b_vector_length("matrix.mat", str(filename))


def test_broadcast_matrix_filename_to_match_b_vector_length(tmp_path):
    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(6, 10))

    assert (
        broadcast_matrix_filename_to_match_b_vector_length(
            "matrix.mat", str(tmp_path / "vectors.bvec")
        )
        == ["matrix.mat"] * 10
    )


def test_rotate_b_vectors_size_error(tmp_path):
    """Test that a RuntimeError is raised when the number of b-vectors
    is different from the number of matrices in the list.
    """
    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(6, 10))

    with pytest.raises(RuntimeError, match="Number of b-vectors"):
        rotate_b_vectors(str(tmp_path / "vectors.bvec"), ["mat1", "mat2", "mat3"])


def test_rotate_b_vectors_missing_b_vectors_file(tmp_path):
    """Test that a FileNotFoundError is raised when the b-vectors file is missing."""
    with pytest.raises(FileNotFoundError):
        rotate_b_vectors(str(tmp_path / "foo.bvec"), ["mat1", "mat2"])


def test_rotate_b_vectors_missing_rotation_matrix_file(tmp_path):
    """Test that a FileNotFoundError is raised when one of the matrix file is missing."""
    np.savetxt(tmp_path / "vectors.bvec", np.random.rand(3, 2))
    np.savetxt(tmp_path / "mat1", np.eye(4))

    with pytest.raises(
        FileNotFoundError,
        match="mat2 not found",
    ):
        rotate_b_vectors(
            str(tmp_path / "vectors.bvec"),
            [str(tmp_path / "mat1"), str(tmp_path / "mat2")],
        )


def test_rotate_b_vectors_all_zeros(tmp_path):
    """Test that B-vectors composed only of zeros aren't rotated."""
    np.savetxt(tmp_path / "vectors.bvec", np.zeros((3, 2)))
    for i in range(1, 3):
        np.savetxt(tmp_path / f"mat{i}", np.random.rand(3, 3))

    rotated_b_vectors_file = rotate_b_vectors(
        str(tmp_path / "vectors.bvec"), [str(tmp_path / "mat1"), str(tmp_path / "mat2")]
    )

    assert rotated_b_vectors_file == str(tmp_path / "vectors_rotated.bvec")
    assert_array_equal(np.loadtxt(rotated_b_vectors_file), np.zeros((3, 2)))


def test_rotate_b_vectors_with_identity(tmp_path):
    """Test that b-vectors rotated with identity matrix are identical."""
    origin_b_vectors = np.random.rand(3, 2)
    normalized_origin_b_vectors = np.array(
        [b / np.linalg.norm(b) for b in origin_b_vectors.T]
    ).T
    np.savetxt(tmp_path / "vectors.bvec", normalized_origin_b_vectors)
    for i in range(1, 3):
        np.savetxt(tmp_path / f"mat{i}", np.eye(4))

    rotated_b_vectors_file = rotate_b_vectors(
        str(tmp_path / "vectors.bvec"), [str(tmp_path / "mat1"), str(tmp_path / "mat2")]
    )

    assert rotated_b_vectors_file == str(tmp_path / "vectors_rotated.bvec")
    assert_array_almost_equal(
        np.loadtxt(rotated_b_vectors_file), normalized_origin_b_vectors
    )


def test_rotate_b_vectors_wrong_content_in_vector_file(tmp_path):
    """Test that a ValueError is raised when vector file isn't valid."""
    b_vectors_file = tmp_path / "foo.bvec"
    b_vectors_file.write_text("fooo\nbar")

    with pytest.raises(ValueError, match="could not convert string 'fooo' to float64"):
        rotate_b_vectors(b_vectors_file, ["mat1"])


@pytest.mark.parametrize(
    "checkpoint, dir_to_del, expected",
    [
        (
            Path(
                "ed42125b6abc244af649ff5eff1df41b/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [True, True],
        ),
        (
            Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [False, True],
        ),
        (
            Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name", "node_2_name"],
            [False, False],
        ),
    ],
)
def test_delete_temp_dirs(tmp_path, checkpoint, dir_to_del, expected):
    checkpoint = tmp_path / checkpoint
    base_dir = tmp_path / Path("wd")
    dirs = [
        tmp_path
        / Path(f"wd/pipeline_name/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/{name}")
        for name in ["node_1_name", "node_2_name"]
    ]
    for dir in dirs:
        if not dir.is_dir():
            os.makedirs(dir)
    delete_temp_dirs(checkpoint, dir_to_del, base_dir)
    assert [dir.is_dir() for dir in dirs] == expected
