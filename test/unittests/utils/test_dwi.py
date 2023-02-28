import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal


def test_rename_files_errors(tmp_path):
    from clinica.utils.dwi import rename_files

    with pytest.raises(
        ValueError,
        match=(
            "Could not extract the BIDS identifier "
            "from the DWI input filename foo.txt"
        ),
    ):
        rename_files("foo.txt", {})
    (tmp_path / "foo/bar").mkdir(parents=True)
    with pytest.raises(
        Exception,
        match=(
            "The 'in_file' trait of a RenameInputSpec instance must be a "
            "pathlike object or string representing an existing file"
        ),
    ):
        rename_files(
            "sub-01_ses-M000_dwi_space-b0_preproc.bval",
            {
                str(tmp_path / "foo/bar/baz.nii.gz"): "_space-b0_preproc.nii.gz",
            },
        )


def test_rename_files(tmp_path):
    from clinica.utils.dwi import rename_files

    assert rename_files("sub-01_ses-M000_dwi_space-b0_preproc.bval", {}) == ()
    (tmp_path / "foo/bar").mkdir(parents=True)
    (tmp_path / "foo/bar/baz.nii.gz").touch()
    (tmp_path / "foo/bar/sub-02_ses-M666_dwi_foo.nii.gz").touch()
    assert rename_files(
        "sub-01_ses-M000_dwi_space-b0_preproc.bval",
        {
            str(tmp_path / "foo/bar/baz.nii.gz"): "_space-b0_preproc.nii.gz",
            str(
                tmp_path / "foo/bar/sub-02_ses-M666_dwi_foo.nii.gz"
            ): "_space-b0_fwhm-4_fmap.nii.gz",
        },
    ) == (
        str(tmp_path / "foo/bar/sub-01_ses-M000_dwi_space-b0_preproc.nii.gz"),
        str(tmp_path / "foo/bar/sub-01_ses-M000_dwi_space-b0_fwhm-4_fmap.nii.gz"),
    )
    assert (tmp_path / "foo/bar/sub-01_ses-M000_dwi_space-b0_preproc.nii.gz").exists()
    assert (
        tmp_path / "foo/bar/sub-01_ses-M000_dwi_space-b0_fwhm-4_fmap.nii.gz"
    ).exists()
    assert (tmp_path / "foo/bar/baz.nii.gz").exists()
    assert (tmp_path / "foo/bar/sub-02_ses-M666_dwi_foo.nii.gz").exists()


@pytest.mark.parametrize("phase", ["a", "X", "x+", "-x", "i", "j", "k", "foo"])
def test_generate_acq_file_errors(tmp_path, phase):
    from clinica.utils.dwi import generate_acq_file

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
    from clinica.utils.dwi import generate_acq_file

    data = 5.0 * np.ones((5, 5, 5, 10))
    dwi = nib.Nifti1Image(data, affine=np.eye(4))
    nib.save(dwi, tmp_path / "dwi.nii.gz")
    acq_file = generate_acq_file(
        tmp_path / "dwi.nii.gz", phase, "16", image_id=image_id
    )
    if image_id:
        assert acq_file == str(tmp_path / f"{image_id}_acq.txt")
    else:
        assert acq_file == str(tmp_path / "acq.txt")
    with open(acq_file, "r") as fp:
        lines = fp.readlines()
    assert lines == expected


def test_generate_index_file_bvalue_file_error(tmp_path):
    from clinica.utils.dwi import generate_index_file

    with pytest.raises(
        FileNotFoundError,
        match="Unable to find b-values file",
    ):
        generate_index_file(str(tmp_path / "foo.txt"))


@pytest.mark.parametrize("image_id", [None, "foo", "foo_bar"])
def test_generate_index_file(tmp_path, image_id):
    from clinica.utils.dwi import generate_index_file

    np.savetxt(tmp_path / "foo.bval", [0] + [1000] * 7)
    index_file = generate_index_file(str(tmp_path / "foo.bval"), image_id=image_id)
    if image_id:
        assert index_file == str(tmp_path / f"{image_id}_index.txt")
    else:
        assert index_file == str(tmp_path / "index.txt")
    index = np.loadtxt(index_file)
    assert_array_equal(index, np.ones(8))


def test_get_b0_filter_error(tmp_path):
    from clinica.utils.dwi import get_b0_filter

    with pytest.raises(
        FileNotFoundError,
        match="Cannot find bval file",
    ):
        get_b0_filter(tmp_path / "foo.bval")


def test_get_b0_filter(tmp_path):
    from clinica.utils.dwi import get_b0_filter

    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    assert_array_equal(get_b0_filter(tmp_path / "foo.bval"), np.array([2, 3, 4, 7]))
    assert_array_equal(get_b0_filter(tmp_path / "foo.bval", low_bval=-1), np.array([]))
    assert_array_equal(
        get_b0_filter(tmp_path / "foo.bval", low_bval=500), np.array([2, 3, 4, 7])
    )
    assert_array_equal(
        get_b0_filter(tmp_path / "foo.bval", low_bval=1000), np.arange(8)
    )
    assert_array_equal(
        get_b0_filter(tmp_path / "foo.bval", low_bval=1001), np.arange(8)
    )


def test_count_b0s(tmp_path):
    from clinica.utils.dwi import count_b0s

    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    assert count_b0s(tmp_path / "foo.bval") == 4
    assert count_b0s(tmp_path / "foo.bval", low_bval=-1) == 0
    assert count_b0s(tmp_path / "foo.bval", low_bval=500) == 4
    assert count_b0s(tmp_path / "foo.bval", low_bval=1000) == 8
    assert count_b0s(tmp_path / "foo.bval", low_bval=1001) == 8


@pytest.mark.parametrize("extension", ["nii", "nii.gz"])
def test_compute_average_b0(tmp_path, extension):
    from clinica.utils.dwi import compute_average_b0

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
    assert result.shape == (5, 5, 5)
    expected = np.zeros((5, 5, 5))
    expected[2:4, 2:4, 2:4] = 1.5
    assert_array_equal(result.get_fdata(), expected)

    # With filtering
    out_file = compute_average_b0(
        tmp_path / f"foo.{extension}",
        tmp_path / "foo.bval",
    )
    assert out_file == tmp_path / f"foo_avg_b0.{extension}"
    result = nib.load(out_file)
    assert result.shape == (5, 5, 5)
    expected = np.zeros((5, 5, 5))
    expected[2:4, 2:4, 2:4] = 1.5
    assert_array_equal(result.get_fdata(), expected)


def test_b0_dwi_split_errors(tmp_path):
    from clinica.utils.dwi import DWIDataset, b0_dwi_split

    dwi_dataset = DWIDataset(
        dwi=tmp_path / "foo.nii.gz",
        b_values=tmp_path / "foo.bval",
        b_vectors=tmp_path / "foo.bvec",
    )
    for filename in ("foo.nii.gz", "foo.bval", "foo.bvec"):
        with pytest.raises(
            FileNotFoundError,
            match=f"File {tmp_path / filename} could not be found.",
        ):
            b0_dwi_split(dwi_dataset)
        (tmp_path / filename).touch()
    with pytest.raises(
        ValueError,
        match="low_bval should be >=0. You provided -1.",
    ):
        b0_dwi_split(dwi_dataset, low_bval=-1)


@pytest.mark.parametrize("extension", ["nii", "nii.gz"])
def test_b0_dwi_split(tmp_path, extension):
    from clinica.utils.dwi import DWIDataset, b0_dwi_split

    bvecs_data = np.random.random((3, 8))
    np.savetxt(tmp_path / "foo.bval", [1000, 1000, 0, 0, 0, 1000, 1000, 0])
    np.savetxt(tmp_path / "foo.bvec", bvecs_data)
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
    small_b_dataset, large_b_dataset = b0_dwi_split(dwi_dataset)
    assert small_b_dataset.dwi == tmp_path / f"foo_b0.{extension}"
    assert small_b_dataset.b_values is None
    assert small_b_dataset.b_vectors is None
    assert large_b_dataset.dwi == tmp_path / f"dwi.{extension}"
    assert large_b_dataset.b_values == tmp_path / "bvals"
    assert large_b_dataset.b_vectors == tmp_path / "bvecs"
    b0_img = nib.load(small_b_dataset.dwi)
    expected = np.zeros((5, 5, 5, 4))
    expected[2:4, 2:4, 2:4, 0:2] = 1.0
    expected[2:4, 2:4, 2:4, 2:4] = 2.0
    assert_array_equal(b0_img.get_fdata(), expected)
    dwi = nib.load(large_b_dataset.dwi)
    assert_array_equal(dwi.get_fdata(), expected)
    bvals = np.loadtxt(large_b_dataset.b_values)
    assert_array_equal(bvals, np.array([1000] * 4))
    bvecs = np.loadtxt(large_b_dataset.b_vectors)
    assert_array_almost_equal(bvecs, bvecs_data[:, np.array([0, 1, 5, 6])], decimal=5)


def test_insert_b0_into_dwi(tmp_path):
    from clinica.utils.dwi import DWIDataset, insert_b0_into_dwi

    b0_data = 6.0 * np.ones((5, 5, 5, 1))
    b0_img = nib.Nifti1Image(b0_data, affine=np.eye(4))
    nib.save(b0_img, tmp_path / "b0.nii.gz")
    dwi_data = 4.0 * np.ones((5, 5, 5, 9))
    dwi_img = nib.Nifti1Image(dwi_data, affine=np.eye(4))
    nib.save(dwi_img, tmp_path / "dwi.nii.gz")
    np.savetxt(tmp_path / "bvals", [1000] * 9)
    bvecs_data = np.random.random((3, 9))
    np.savetxt(tmp_path / "bvecs", bvecs_data)
    dwi_dataset = DWIDataset(
        dwi=tmp_path / "dwi.nii.gz",
        b_values=tmp_path / "bvals",
        b_vectors=tmp_path / "bvecs",
    )
    out_dataset = insert_b0_into_dwi(tmp_path / "b0.nii.gz", dwi_dataset)
    assert out_dataset.dwi == tmp_path / "merged_files.nii.gz"
    assert out_dataset.b_values == tmp_path / "bvals"
    assert out_dataset.b_vectors == tmp_path / "bvecs"
    dwi = nib.load(out_dataset.dwi)
    assert_array_equal(dwi.affine, dwi_img.affine)
    expected = 4.0 * np.ones((5, 5, 5, 10))
    expected[..., 0] += 2.0
    assert_array_equal(dwi.get_fdata(), expected)
    bvals = np.loadtxt(out_dataset.b_values)
    assert_array_equal(bvals, np.array([0] + [1000] * 9))
    bvecs = np.loadtxt(out_dataset.b_vectors)
    expected = np.insert(bvecs_data, 0, 0.0, axis=1)
    assert_array_almost_equal(bvecs, expected, decimal=5)
