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
