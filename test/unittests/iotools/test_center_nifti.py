import json
from os import PathLike
from pathlib import Path
from typing import Iterable, Optional, Union
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pytest

from clinica.utils.exceptions import ClinicaBIDSError, ClinicaExistingDatasetError


def build_bids_folder(tmp_path: Path) -> Path:
    bids_path = tmp_path / "BIDS"
    bids_path.mkdir()
    with open(bids_path / "dataset_description.json", "w") as fp:
        json.dump({"Name": "Test", "DatasetType": "raw", "BIDSVersion": "1.7.0"}, fp)
    (bids_path / "sub-1" / "ses-1").mkdir(parents=True)
    (bids_path / "sub-1" / "ses-1" / "sub-1_ses-1_T1w.nii.gz").touch()
    (bids_path / "sub-2" / "ses-1").mkdir(parents=True)
    (bids_path / "sub-2" / "ses-1" / "sub-2_ses-1_T1w.nii.gz").touch()
    (bids_path / "sub-2" / "ses-2").mkdir(parents=True)
    (bids_path / "sub-2" / "ses-2" / "sub-2_ses-2_pet.nii.gz").touch()
    (bids_path / "sub-2" / "ses-2" / "sub-2_ses-2_T1toto.nii.gz").touch()

    return bids_path


def test_handle_output_existing_files(tmp_path):
    from clinica.iotools.data_handling._centering import (
        _handle_output_existing_files,
    )

    output_path = tmp_path / "COPY"
    assert _handle_output_existing_files(output_path) == output_path

    output_path = build_bids_folder(tmp_path)
    with pytest.raises(ClinicaExistingDatasetError):
        _handle_output_existing_files(output_path)

    _handle_output_existing_files(output_path, overwrite_existing_files=True)
    assert not output_path.exists()


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
def test_find_files_with_modality(tmp_path, modalities, expected: set[str]):
    from clinica.iotools.data_handling._centering import _find_files_with_modality

    bids_path = build_bids_folder(tmp_path)
    result = _find_files_with_modality(bids_path, modalities)

    assert set(r.name for r in result) == expected


def test_center_nifti_error(tmp_path):
    from clinica.iotools import center_nifti

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
    from clinica.iotools import center_nifti

    bids_path = build_bids_folder(tmp_path)

    out_path = tmp_path / "out"
    out_path.mkdir()

    mocker.patch(
        "clinica.iotools.data_handling.write_list_of_files", return_value=False
    )
    with patch(
        "clinica.iotools.data_handling.center_all_nifti",
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
        (50, ("bar.tsv", "baz.nii.gz", "foo.txt")),
        (20, ("bar.tsv", "baz.nii.gz", "foo.txt", "foobar.nii.gz")),
        (80, ("bar.tsv", "foo.txt")),
    ],
)
def test_center_nifti(tmp_path, centering_threshold, expected_files):
    from clinica.iotools import center_nifti

    bids_path = build_bids_folder(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with patch(
        "clinica.iotools.data_handling.center_all_nifti",
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
        assert logs[0].read_text() == "\n".join([f for f in expected_files])


@pytest.mark.parametrize(
    "quatern, pixdim, qoffset",
    [
        # 90 deg rotation around X with anisotropic voxel sizes
        ((np.sqrt(0.5), 0.0, 0.0), (1.0, 1.0, 1.0, 2.0), (-128.0, 160.0, -128.0)),
        # General 3D rotation with qfac = -1 and anisotropic voxel sizes
        ((0.3, -0.4, 0.5), (-1.0, 1.2, 0.8, 1.5), (10.0, -20.0, 30.0)),
    ],
)
def test_qform_rotation_and_world_center_match_nibabel(
    tmp_path, quatern, pixdim, qoffset
):
    from clinica.iotools.data_handling._centering import (
        _get_rotation_matrix,
        _get_world_coordinate_of_center,
        _scale_with_rotation_matrix,
    )

    header = nib.Nifti1Header()
    header.set_data_shape((256, 256, 160))
    header["qform_code"] = 1
    header["sform_code"] = 0
    header["pixdim"][0:4] = pixdim
    header["quatern_b"], header["quatern_c"], header["quatern_d"] = quatern
    header["qoffset_x"], header["qoffset_y"], header["qoffset_z"] = qoffset

    rot = _get_rotation_matrix(header)
    np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-6)
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-6)

    qform = header.get_qform()
    center_vox = np.array([128.0, 128.0, 80.0])
    expected_world = (qform @ np.array([*center_vox, 1.0]))[:3]

    np.testing.assert_allclose(
        _scale_with_rotation_matrix(center_vox, header),
        expected_world,
        atol=1e-4,
    )

    img = nib.Nifti1Image(np.zeros((256, 256, 160), dtype=np.uint8), None, header)
    img.set_sform(np.eye(4), code=0)
    img.set_qform(qform, code=1)
    nii_path = tmp_path / "rotated.nii.gz"
    nib.save(img, nii_path)

    np.testing.assert_allclose(
        _get_world_coordinate_of_center(nii_path),
        expected_world,
        atol=1e-4,
    )
