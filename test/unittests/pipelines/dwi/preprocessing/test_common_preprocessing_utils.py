import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal


@pytest.mark.parametrize("phase", ["a", "X", "x+", "-x", "i", "j", "k", "foo"])
def test_generate_acq_file_errors(tmp_path, phase):
    from clinica.pipelines.dwi.preprocessing.utils import generate_acq_file

    data = 5.0 * np.ones((5, 5, 5, 10))
    dwi = nib.Nifti1Image(data, affine=np.eye(4))
    nib.save(dwi, tmp_path / "dwi.nii.gz")
    with pytest.raises(RuntimeError, match="FSL PhaseEncodingDirection"):
        generate_acq_file(tmp_path / "dwi.nii.gz", phase, "16")


@pytest.mark.parametrize("image_id", [None, "foo"])
@pytest.mark.parametrize(
    "phase,expected",
    [
        ("y", ["0 1 0 16.000000\n"]),
        ("y-", ["0 -1 0 16.000000\n"]),
        ("x", ["1 0 0 16.000000\n"]),
        ("x-", ["-1 0 0 16.000000\n"]),
        ("z", ["0 0 1 16.000000\n"]),
        ("z-", ["0 0 -1 16.000000\n"]),
    ],
)
def test_generate_acq_file(tmp_path, image_id, phase, expected):
    from clinica.pipelines.dwi.preprocessing.utils import generate_acq_file

    data = 5.0 * np.ones((5, 5, 5, 10))
    dwi = nib.Nifti1Image(data, affine=np.eye(4))
    nib.save(dwi, tmp_path / "dwi.nii.gz")
    acq_file = generate_acq_file(
        tmp_path / "dwi.nii.gz", phase, "16", image_id=image_id
    )
    filename = f"{image_id}_acq.txt" if image_id else "acq.txt"
    assert acq_file == tmp_path / filename
    with acq_file.open() as f:
        lines = f.readlines()
    assert lines == expected


def test_generate_index_file_bvalue_file_error(tmp_path):
    from clinica.pipelines.dwi.preprocessing.utils import generate_index_file

    with pytest.raises(
        FileNotFoundError,
        match="Unable to find b-values file",
    ):
        generate_index_file(tmp_path / "foo.txt")


@pytest.mark.parametrize("image_id", [None, "foo", "foo_bar"])
def test_generate_index_file(tmp_path, image_id):
    from clinica.pipelines.dwi.preprocessing.utils import generate_index_file

    np.savetxt(tmp_path / "foo.bval", [0] + [1000] * 7)
    index_file = generate_index_file(tmp_path / "foo.bval", image_id=image_id)
    if image_id:
        assert index_file == tmp_path / f"{image_id}_index.txt"
    else:
        assert index_file == tmp_path / "index.txt"
    index = np.loadtxt(index_file)
    assert_array_equal(index, np.ones(8))


def test_get_b0_filter_error(tmp_path):
    from clinica.pipelines.dwi.preprocessing.utils import get_b0_filter

    with pytest.raises(
        FileNotFoundError,
        match="File not found",
    ):
        get_b0_filter(tmp_path / "foo.bval")


@pytest.mark.parametrize(
    "threshold,expected",
    [
        (None, np.array([2, 3, 4, 7])),
        (-1, np.array([])),
        (500, np.array([2, 3, 4, 7])),
        (1000, np.arange(8)),
        (1001, np.arange(8)),
    ],
)
def test_get_b0_filter(tmp_path, threshold, expected):
    from clinica.pipelines.dwi.preprocessing.utils import get_b0_filter

    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    kwargs = {"b_value_threshold": threshold} if threshold else {}

    assert_array_equal(get_b0_filter(tmp_path / "foo.bval", **kwargs), expected)


@pytest.mark.parametrize("extension", ["nii", "nii.gz"])
def test_compute_average_b0(tmp_path, extension):
    from clinica.pipelines.dwi.preprocessing.utils import compute_average_b0

    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    img_data = np.zeros((5, 5, 5, 8))
    img_data[2:4, 2:4, 2:4, 0:4] = 1.0
    img_data[2:4, 2:4, 2:4, 4:8] = 2.0
    img = nib.Nifti1Image(img_data, affine=np.eye(4))
    nib.save(img, tmp_path / f"foo.{extension}")

    # No filtering with the bvalues
    out_file = compute_average_b0(tmp_path / f"foo.{extension}")
    assert out_file == tmp_path / f"foo_avg_b0.{extension}"
    result = nib.load(out_file)
    assert result.shape == (5, 5, 5, 1)
    expected = np.zeros((5, 5, 5, 1))
    expected[2:4, 2:4, 2:4] = 1.5
    assert_array_equal(result.get_fdata(), expected)

    # With filtering
    out_file = compute_average_b0(
        tmp_path / f"foo.{extension}",
        tmp_path / "foo.bval",
    )
    assert out_file == tmp_path / f"foo_avg_b0.{extension}"
    result = nib.load(out_file)
    assert result.shape == (5, 5, 5, 1)
    expected = np.zeros((5, 5, 5, 1))
    expected[2:4, 2:4, 2:4] = 1.5
    assert_array_equal(result.get_fdata(), expected)


@pytest.fixture
def dwi_dataset(tmp_path):
    from clinica.pipelines.dwi.utils import DWIDataset

    return DWIDataset(
        dwi=str(tmp_path / "foo.nii.gz"),
        b_values=tmp_path / "foo.bval",
        b_vectors=tmp_path / "foo.bvec",
    )


def test_check_dwi_dataset(tmp_path, dwi_dataset):
    from clinica.pipelines.dwi.preprocessing.utils import check_dwi_dataset

    for filename in ("foo.nii.gz", "foo.bval", "foo.bvec"):
        with pytest.raises(
            FileNotFoundError,
            match="File not found",
        ):
            check_dwi_dataset(dwi_dataset)
        (tmp_path / filename).touch()
    dwi_dataset_checked = check_dwi_dataset(dwi_dataset)
    assert len(dwi_dataset_checked) == 3
    assert dwi_dataset_checked.dwi == tmp_path / "foo.nii.gz"
    assert dwi_dataset_checked.b_values == tmp_path / "foo.bval"
    assert dwi_dataset_checked.b_vectors == tmp_path / "foo.bvec"


@pytest.mark.parametrize(
    "n_dwi,n_b_values,n_b_vectors", [(10, 9, 9), (9, 10, 9), (9, 9, 10), (8, 9, 10)]
)
def test_check_dwi_volume_errors(tmp_path, n_dwi, n_b_values, n_b_vectors):
    from clinica.pipelines.dwi.preprocessing.utils import check_dwi_volume
    from clinica.utils.testing_utils import build_dwi_dataset

    dwi_dataset = build_dwi_dataset(tmp_path, n_dwi, n_b_values, n_b_vectors)
    with pytest.raises(
        IOError,
        match="Number of DWIs, b-vals and b-vecs mismatch",
    ):
        check_dwi_volume(dwi_dataset)


def test_check_dwi_volume(tmp_path, dwi_dataset):
    from clinica.pipelines.dwi.preprocessing.utils import check_dwi_volume
    from clinica.utils.testing_utils import build_dwi_dataset

    check_dwi_volume(build_dwi_dataset(tmp_path, 9, 9, 9))
