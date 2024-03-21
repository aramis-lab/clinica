from pathlib import Path

import pytest


def test_get_luts(mocker):
    from clinica.pipelines.dwi.connectome.utils import get_luts

    mocked_freesurfer_home = Path("/Applications/freesurfer/7.2.0")
    mocker.patch(
        "clinica.utils.check_dependency.get_freesurfer_home",
        return_value=mocked_freesurfer_home,
    )
    assert get_luts() == [f"{mocked_freesurfer_home}/FreeSurferColorLUT.txt"] * 2


@pytest.mark.parametrize(
    "filename,expected_checksum",
    [
        (
            "fs_default.txt",
            "a8d561694887a1ca8d9df223aa5ef861b6c79d43ce9ed93835b9ce8aadc331b1",
        ),
        (
            "fs_a2009s.txt",
            "40b0d4d77bde7e1d265439347af5b30cc973748c1a88d203d7044cb35b3863e1",
        ),
    ],
)
def test_get_checksum_for_filename(filename, expected_checksum):
    from clinica.pipelines.dwi.connectome.utils import _get_checksum_for_filename  # noqa

    assert _get_checksum_for_filename(filename) == expected_checksum


def test_get_checksum_for_filename_error():
    from clinica.pipelines.dwi.connectome.utils import _get_checksum_for_filename  # noqa

    with pytest.raises(ValueError, match="File name foo.txt is not supported."):
        _get_checksum_for_filename("foo.txt")


@pytest.mark.parametrize(
    "filename,expected_length", [("fs_default.txt", 112), ("fs_a2009s.txt", 192)]
)
def test_download_mrtrix3_file(tmp_path, filename, expected_length):
    """Atm this test needs an internet connection to download the files.

    TODO: Use mocking in the fetch_file function to remove this necessity.
    """
    from clinica.pipelines.dwi.connectome.utils import _download_mrtrix3_file  # noqa

    _download_mrtrix3_file(filename, tmp_path)

    assert [f.name for f in tmp_path.iterdir()] == [filename]
    assert len((tmp_path / filename).read_text().split("\n")) == expected_length


def test_download_mrtrix3_file_error(tmp_path, mocker):
    import re

    from clinica.pipelines.dwi.connectome.utils import _download_mrtrix3_file  # noqa

    mocker.patch(
        "clinica.pipelines.dwi.connectome.utils._get_checksum_for_filename",
        return_value="foo",
    )
    mocker.patch("clinica.utils.inputs.fetch_file", side_effect=IOError)

    with pytest.raises(
        IOError,
        match=re.escape(
            "Unable to download required MRTRIX mapping (foo.txt) for processing"
        ),
    ):
        _download_mrtrix3_file("foo.txt", tmp_path)


def test_get_conversion_luts():
    from pathlib import Path

    from clinica.pipelines.dwi.connectome.utils import get_conversion_luts

    luts = [Path(_) for _ in get_conversion_luts()]

    assert [p.name for p in luts] == ["fs_default.txt", "fs_a2009s.txt"]
    assert all([p.is_file() for p in luts])


@pytest.mark.parametrize(
    "filename",
    [
        "foo.txt",
        "dwi.nii.gz",
        "sub-01_ses-M000_dwi.nii.gz",
        "sub-01_ses-M000_preproc.nii.gz",
        "sub-01_ses-M000_space-T1w_preproc.nii.gz",
        "sub-01_ses-M000_space-b0_preproc.nii.gz",
    ],
)
def test_get_caps_filenames_error(tmp_path, filename):
    from clinica.pipelines.dwi.connectome.utils import get_caps_filenames

    with pytest.raises(ValueError, match="is not in a CAPS compliant format."):
        get_caps_filenames(str(tmp_path / filename))


def test_get_caps_filenames(tmp_path):
    from clinica.pipelines.dwi.connectome.utils import get_caps_filenames

    dwi_caps = tmp_path / "dwi" / "preprocessing"
    dwi_caps.mkdir(parents=True)

    assert get_caps_filenames(
        str(dwi_caps / "sub-01_ses-M000_space-b0_desc-preproc_dwi.nii.gz")
    ) == (
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_responseFunction.txt",
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_diffmodel.nii.gz",
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_tractography.tck",
        [
            "sub-01_ses-M000_space-b0_desc-preproc_atlas-desikan_parcellation.nii.gz",
            "sub-01_ses-M000_space-b0_desc-preproc_atlas-destrieux_parcellation.nii.gz",
        ],
        [
            "sub-01_ses-M000_model-CSD_atlas-desikan_connectivity.tsv",
            "sub-01_ses-M000_model-CSD_atlas-destrieux_connectivity.tsv",
        ],
    )


@pytest.mark.parametrize(
    "subjects,sessions,expected",
    [
        ([], [], []),
        (["foo"], ["bar"], ["subjects/foo/bar/dwi"]),
        (
            ["sub-01", "sub-02"],
            ["ses-M000", "ses-M006"],
            ["subjects/sub-01/ses-M000/dwi", "subjects/sub-02/ses-M006/dwi"],
        ),
    ],
)
def test_get_containers(subjects, sessions, expected):
    from clinica.pipelines.dwi.connectome.utils import get_containers

    assert get_containers(subjects, sessions) == expected
