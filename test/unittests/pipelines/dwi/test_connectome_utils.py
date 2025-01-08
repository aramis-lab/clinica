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
            "bfebee26de22dc4cd03d5ee3f26524b046cce232679e1ba1bc26f18180d491f1",
        ),
        (
            "fs_a2009s.txt",
            "ae9660f2a9fb44b7d828dcf1f390ce81ed600471810af89042ba011c7a2a675f",
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


def _create_files_helper(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    flirt_matrix = tmp_path / "flirt.mat"
    flirt_matrix.touch()

    source_image = tmp_path / "src.nii.gz"
    source_image.touch()

    reference_image = tmp_path / "ref.nii.gz"
    reference_image.touch()

    mrtrix_matrix = tmp_path / "mrtrix_matrix.mat"
    mrtrix_matrix.touch()

    return source_image, reference_image, flirt_matrix, mrtrix_matrix


def test_get_transformation_cmd(tmp_path):
    from clinica.pipelines.dwi.connectome.utils import _get_transformation_cmd

    source_image, reference_image, flirt_matrix, mrtrix_matrix = _create_files_helper(
        tmp_path
    )

    assert (
        _get_transformation_cmd(
            source_image, reference_image, flirt_matrix, mrtrix_matrix
        )
        == f"transformconvert {flirt_matrix} {source_image} {reference_image} flirt_import {mrtrix_matrix}"
    )


def test_get_transformation_cmd_error(tmp_path):
    from clinica.pipelines.dwi.connectome.utils import _get_transformation_cmd

    with pytest.raises(
        FileNotFoundError, match=f"The file {tmp_path / 'src.nii.gz'} was not found"
    ):
        _get_transformation_cmd(
            tmp_path / "src.nii.gz",
            tmp_path / "ref.nii.gz",
            tmp_path / "flirt.mat",
            tmp_path / "mrtx.mat",
        )


def test_convert_flirt_to_mrtrix_transformation(tmp_path, mocker):
    from clinica.pipelines.dwi.connectome.utils import (
        convert_flirt_to_mrtrix_transformation,
    )

    mocker.patch("clinica.utils.check_dependency.check_software", return_value=None)
    mocker.patch("os.system", return_value=None)

    source_image, reference_image, flirt_matrix, output_matrix = _create_files_helper(
        tmp_path
    )

    assert (
        convert_flirt_to_mrtrix_transformation(
            source_image, reference_image, flirt_matrix, output_matrix.name
        )
        == Path(output_matrix.name).resolve()
    )
    assert (
        convert_flirt_to_mrtrix_transformation(
            source_image, reference_image, flirt_matrix, None
        )
        == Path("mrtrix_matrix.mat").resolve()
    )
