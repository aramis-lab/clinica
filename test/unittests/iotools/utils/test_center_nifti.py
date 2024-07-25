from os import PathLike
from pathlib import Path
from typing import Iterable, List, Optional
from unittest.mock import patch

import pytest


def test_center_nifti_error(tmp_path):
    from clinica.iotools.utils import center_nifti
    from clinica.utils.exceptions import ClinicaExistingDatasetError

    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "foo.txt").touch()

    with pytest.raises(
        ClinicaExistingDatasetError,
        match=f"Dataset located at {tmp_path / 'out'} already contain some files.",
    ):
        center_nifti(tmp_path / "bids", tmp_path / "out")


def center_all_nifti_mock(
    bids_dir: PathLike,
    output_dir: PathLike,
    modalities: Optional[Iterable[str]] = None,
    center_all_files: bool = False,
) -> List[Path]:
    result = []
    for filename in ("foo.txt", "bar.tsv", "baz.nii.gz"):
        (Path(output_dir) / filename).touch()
        result.append(Path(output_dir) / filename)
    return result


def test_center_nifti_log_file_creation_error(tmp_path, mocker):
    from clinica.iotools.utils import center_nifti

    (tmp_path / "bids").mkdir()
    (tmp_path / "out").mkdir()
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
            center_nifti(tmp_path / "bids", tmp_path / "out", center_all_files=True)
        mock.assert_called_once_with(tmp_path / "bids", tmp_path / "out", None, True)


def test_center_nifti(tmp_path):
    from clinica.iotools.utils import center_nifti

    (tmp_path / "bids").mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with patch(
        "clinica.iotools.utils.data_handling.center_all_nifti",
        wraps=center_all_nifti_mock,
    ) as mock:
        center_nifti(tmp_path / "bids", tmp_path / "out", ("anat", "pet", "dwi"))
        mock.assert_called_once_with(
            tmp_path / "bids", tmp_path / "out", ("anat", "pet", "dwi"), False
        )
        assert (
            len([f for f in output_dir.iterdir()]) == 4
        )  # 3 files written by mock and 1 log file
        logs = [f for f in output_dir.glob("centered_nifti_list_*.txt")]
        assert len(logs) == 1
        assert logs[0].read_text() == "\n".join(
            [str(tmp_path / "out" / f) for f in ("foo.txt", "bar.tsv", "baz.nii.gz")]
        )
