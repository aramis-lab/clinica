from os import PathLike
from pathlib import Path
from typing import Iterable, Optional, Union
from unittest.mock import patch

import pytest

from clinica.utils.exceptions import ClinicaBIDSError, ClinicaExistingDatasetError


def build_bids_folder(tmp_path: Path) -> Path:
    bids_path = tmp_path / "BIDS"
    bids_path.mkdir()

    (bids_path / "dataset_description.json").touch()

    (bids_path / "sub-1" / "ses-1").mkdir(parents=True)
    (bids_path / "sub-1" / "ses-1" / "sub-1_ses-1_T1w.nii.gz").touch()
    (bids_path / "sub-2" / "ses-1").mkdir(parents=True)
    (bids_path / "sub-2" / "ses-1" / "sub-2_ses-1_T1w.nii.gz").touch()
    (bids_path / "sub-2" / "ses-2").mkdir(parents=True)
    (bids_path / "sub-2" / "ses-2" / "sub-2_ses-2_pet.nii.gz").touch()
    (bids_path / "sub-2" / "ses-2" / "sub-2_ses-2_T1toto.nii.gz").touch()

    return bids_path


def test_handle_output_existing_files(tmp_path):
    from clinica.iotools.utils.data_handling._centering import (
        _handle_output_existing_files,
    )

    output_path = tmp_path / "COPY"
    assert _handle_output_existing_files(output_path) == output_path

    output_path = build_bids_folder(tmp_path)
    with pytest.raises(ClinicaExistingDatasetError):
        _handle_output_existing_files(output_path)

    _handle_output_existing_files(output_path, overwrite_existing_files=True)
    assert not output_path.exists()


def test_validate_bids_and_output_dir_equal_error(tmp_path):
    from clinica.iotools.utils.data_handling._centering import (
        _validate_bids_and_output_dir,
    )

    with pytest.raises(ClinicaBIDSError):
        _validate_bids_and_output_dir(tmp_path / "BIDS", tmp_path / "BIDS")


def test_validate_bids_and_output_dir_not_bids_error(tmp_path):
    from clinica.iotools.utils.data_handling._centering import (
        _validate_bids_and_output_dir,
    )

    (tmp_path / "BIDS").mkdir()
    with pytest.raises(ClinicaBIDSError):
        _validate_bids_and_output_dir(tmp_path / "BIDS", tmp_path / "COPY")


def test_validate_bids_and_output_dir_not_bids_success(tmp_path):
    from clinica.iotools.utils.data_handling._centering import (
        _validate_bids_and_output_dir,
    )

    bids_path = build_bids_folder(tmp_path)
    output_path = tmp_path / "COPY"
    assert bids_path, output_path == _validate_bids_and_output_dir(
        bids_path, output_path
    )


@pytest.mark.parametrize(
    "modalities, expected",
    [
        (
            None,
            {
                "sub-1_ses-1_T1w.nii.gz",
                "sub-2_ses-1_T1w.nii.gz",
                "sub-2_ses-2_pet.nii.gz",
                "sub-2_ses-2_T1toto.nii.gz",
            },
        ),
        (("T1w",), {"sub-1_ses-1_T1w.nii.gz", "sub-2_ses-1_T1w.nii.gz"}),
        (
            ("T1",),
            {
                "sub-1_ses-1_T1w.nii.gz",
                "sub-2_ses-1_T1w.nii.gz",
                "sub-2_ses-2_T1toto.nii.gz",
            },
        ),
        (
            ("t1w", "pet"),
            {
                "sub-1_ses-1_T1w.nii.gz",
                "sub-2_ses-1_T1w.nii.gz",
                "sub-2_ses-2_pet.nii.gz",
            },
        ),
    ],
)
def test_find_files_with_modality(tmp_path, modalities, expected):
    from clinica.iotools.utils.data_handling._centering import _find_files_with_modality

    bids_path = build_bids_folder(tmp_path)
    result = _find_files_with_modality(bids_path, modalities)
    assert set(r.name for r in result) == expected


def test_center_nifti_error(tmp_path):
    from clinica.iotools.utils import center_nifti

    out_path = tmp_path / "out"
    out_path.mkdir()
    (out_path / "foo.txt").touch()

    bids_path = build_bids_folder(tmp_path)

    with pytest.raises(
        ClinicaExistingDatasetError,
        match=f"Dataset located at {out_path} already contain some files.",
    ):
        center_nifti(bids_path, out_path)


def center_all_nifti_mock(
    bids_dir: Union[str, PathLike],
    output_dir: Union[str, PathLike],
    modalities: Optional[Iterable[str]] = None,
    centering_threshold: int = 50,
    overwrite_existing_file: bool = False,
) -> list[Path]:
    files = ("foo.txt", "bar.tsv", "baz.nii.gz")
    if centering_threshold > 50:
        files = files[:-1]
    if centering_threshold < 50:
        files += ("foobar.nii.gz",)
    result = []
    for filename in files:
        (Path(output_dir) / filename).touch()
        result.append(Path(output_dir) / filename)
    return result


def test_center_nifti_log_file_creation_error(tmp_path, mocker):
    from clinica.iotools.utils import center_nifti

    bids_path = build_bids_folder(tmp_path)

    out_path = tmp_path / "out"
    out_path.mkdir()

    mocker.patch(
        "clinica.iotools.utils.data_handling.write_list_of_files", return_value=False
    )
    with patch(
        "clinica.iotools.utils.data_handling.center_all_nifti",
        wraps=center_all_nifti_mock,
    ) as mock:
        with pytest.raises(
            IOError,
            match="Could not create log file",
        ):
            center_nifti(bids_path, out_path)
        mock.assert_called_once_with(bids_path, out_path, None, 50, False)


@pytest.mark.parametrize(
    "centering_threshold, expected_files",
    [
        (50, ("foo.txt", "bar.tsv", "baz.nii.gz")),
        (20, ("foo.txt", "bar.tsv", "baz.nii.gz", "foobar.nii.gz")),
        (80, ("foo.txt", "bar.tsv")),
    ],
)
def test_center_nifti(tmp_path, centering_threshold, expected_files):
    from clinica.iotools.utils import center_nifti

    bids_path = build_bids_folder(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with patch(
        "clinica.iotools.utils.data_handling.center_all_nifti",
        wraps=center_all_nifti_mock,
    ) as mock:
        center_nifti(bids_path, output_dir, ("anat", "pet", "dwi"), centering_threshold)
        mock.assert_called_once_with(
            bids_path, output_dir, ("anat", "pet", "dwi"), centering_threshold, False
        )
        assert (
            len([f for f in output_dir.iterdir()]) == len(expected_files) + 1
        )  # 3 files written by mock and 1 log file
        logs = [f for f in output_dir.glob("centered_nifti_list_*.txt")]
        assert len(logs) == 1
        assert logs[0].read_text() == "\n".join(
            [str(output_dir / f) for f in expected_files]
        )
