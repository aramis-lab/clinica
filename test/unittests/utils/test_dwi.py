import nibabel as nib
import numpy as np
import pytest


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
