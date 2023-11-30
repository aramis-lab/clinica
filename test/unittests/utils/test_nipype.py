import pytest


@pytest.mark.parametrize(
    "input_path,expected_container",
    [
        (
            "/path/to/bids/sub-CLNC01/ses-M000/anat/sub-CLNC01_ses-M000_T1w.nii.gz",
            "subjects/sub-CLNC01/ses-M000",
        ),
        (
            "caps/subjects/sub-CLNC01/ses-M000/dwi/preprocessing/sub-CLNC01_ses-M000_preproc.nii",
            "subjects/sub-CLNC01/ses-M000",
        ),
        (
            "foo/bar/sub-01/ses-M327/foo/bar/baz/foo.txt",
            "subjects/sub-01/ses-M327",
        ),
        (
            "foo/bar/sub-666/ses-M45678/foo/bar/baz/foo.txt",
            "subjects/sub-666/ses-M45678",
        ),
        (
            "foo/bar/sub-baz/ses-foo/foo/bar/baz/foo.txt",
            "subjects/sub-baz/ses-foo",
        ),
    ],
)
def test_container_from_filename(input_path, expected_container):
    from clinica.utils.nipype import container_from_filename

    assert container_from_filename(input_path) == expected_container


@pytest.mark.parametrize(
    "input_path", ["foo", "sub-01", "ses-M000", "sub-01\\ses-M000", "100"]
)
def test_container_from_filename_errors(input_path):
    from clinica.utils.nipype import container_from_filename

    with pytest.raises(
        ValueError,
        match="Input filename",
    ):
        container_from_filename(input_path)
