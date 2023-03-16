import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal


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


def test_generate_index_file_no_b0_error(tmp_path):
    from clinica.utils.dwi import generate_index_file

    np.savetxt(tmp_path / "foo.bval", [1000] * 8)
    with pytest.raises(
        ValueError,
        match="Could not find b-value <= 5.0 in bval file",
    ):
        generate_index_file(str(tmp_path / "foo.bval"))


@pytest.mark.parametrize("image_id", [None, "foo", "foo_bar"])
def test_generate_index_file(tmp_path, image_id):
    from clinica.utils.dwi import generate_index_file

    #  This will work because there is a single B0 volume at index 0
    #  But a more general input like [1000, 1000, 0, 0, 0, 1000, 1000, 0]
    #  would fail...
    np.savetxt(tmp_path / "foo.bval", [0] + [1000] * 7)
    index_file = generate_index_file(str(tmp_path / "foo.bval"), image_id=image_id)
    if image_id:
        assert index_file == str(tmp_path / f"{image_id}_index.txt")
    else:
        assert index_file == str(tmp_path / "index.txt")
    index = np.loadtxt(index_file)
    # This is wrong and needs to be fixed
    # the array should be length 8...
    # assert_array_equal(index, np.array([1.0, 2.0, 3.0, 3.0, 3.0, 4.0]))
    assert_array_equal(index, np.ones(8))
