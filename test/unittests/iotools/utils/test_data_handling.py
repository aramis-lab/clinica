from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_array_equal


def test_scale_coordinates_by_pixdim():
    """Test function `_scale_coordinates_by_pixdim`."""
    from nibabel.nifti1 import Nifti1Header

    from clinica.iotools.utils.data_handling._centering import (
        _scale_coordinates_by_pixdim,
    )

    head = Nifti1Header()
    assert not head["qform_code"] > 0
    assert head["sform_code"] == 0
    assert np.all(head["pixdim"] == 1)
    center_coordinates = np.array([0.5] * 3)
    assert_array_equal(
        _scale_coordinates_by_pixdim(coordinates_vol=center_coordinates, header=head),
        center_coordinates,
    )
    head["pixdim"][1] = 2
    assert_array_equal(
        _scale_coordinates_by_pixdim(coordinates_vol=center_coordinates, header=head),
        np.array([1.0, 0.5, 0.5]),
    )


def test_check_relative_volume_location_in_world_coordinate_system(tmp_path, mocker):
    from clinica.iotools.utils.data_handling import (
        check_relative_volume_location_in_world_coordinate_system,
    )

    mocker.patch(
        "clinica.iotools.utils.data_handling._centering._get_problematic_pairs_with_l2_norm",
        return_value=[
            ("foo.nii.gz", "bar.nii.gz", 81.2),
            ("baz.nii", "goo.nii", 113.3),
        ],
    )
    with pytest.warns(
        UserWarning,
        match=(
            "It appears that 2 pairs of files have an important relative offset. "
            "SPM co-registration has a high probability to fail on these files:\n\n"
            "    foo_label  gloo_label  Relative distance\n"
            "0  foo.nii.gz  bar.nii.gz               81.2\n"
            "1     baz.nii     goo.nii              113.3\n"
            "Clinica provides a tool to counter this problem by replacing the center "
            "of the volume at the origin of the world coordinates.\n"
            "Use the following command line to correct the header of the faulty NIFTI "
            "volumes in a new folder:\n\n`"
        ),
    ):
        check_relative_volume_location_in_world_coordinate_system(
            "foo_label",
            [tmp_path / "foo.nii.gz", tmp_path / "bar.nii.gz", tmp_path / "baz.nii.gz"],
            "gloo_label",
            [
                tmp_path / "gloo" / "foo.nii.gz",
                tmp_path / "gloo" / "bar.nii.gz",
                tmp_path / "gloo" / "baz.nii.gz",
            ],
            tmp_path / "bids",
            "T1W",
            skip_question=True,
        )


def test_get_problematic_pairs_empty(tmp_path, mocker):
    from clinica.iotools.utils.data_handling._centering import (
        _get_problematic_pairs_with_l2_norm,
    )

    mocker.patch(
        "clinica.iotools.utils.data_handling._centering._get_world_coordinate_of_center",
        return_value=np.zeros((3, 3)),
    )
    assert (
        _get_problematic_pairs_with_l2_norm(
            [
                (tmp_path / "foo.nii.gz", tmp_path / "baz.nii.gz"),
                (tmp_path / "bar.nii.gz", tmp_path / "foo_bar.nii.gz"),
            ]
        )
        == []
    )


def test_get_problematic_pairs(tmp_path, mocker):
    from clinica.iotools.utils.data_handling._centering import (
        _get_problematic_pairs_with_l2_norm,
    )

    mocker.patch(
        "clinica.iotools.utils.data_handling._centering._compute_l2_norm",
        return_value=[81.0, 79.9],
    )
    assert _get_problematic_pairs_with_l2_norm(
        [
            (tmp_path / "foo.nii.gz", tmp_path / "baz.nii.gz"),
            (tmp_path / "bar.nii.gz", tmp_path / "foo_bar.nii.gz"),
        ]
    ) == [("foo.nii.gz", "baz.nii.gz", 81.0)]


def test_build_warning_message(tmp_path):
    from clinica.iotools.utils.data_handling._centering import _build_warning_message

    assert (
        "It appears that 3 files have a center way out of the origin of the world coordinate system. "
        "SPM has a high probability to fail on these files (for co-registration or segmentation):\n\n"
        "      File Coordinate of center  Distance to origin\n"
        "0  foo.nii            [1, 2, 3]                 0.2\n"
        "1  bar.nii            [4, 5, 6]                12.7\n"
        "2  baz.nii            [7, 8, 9]                26.6\n"
        "If you are trying to launch the t1-freesurfer pipeline, you can ignore this message if you "
        "do not want to run the pet-surface pipeline afterward.\n"
        "Clinica provides a tool to counter this problem by replacing the center of the volume at the "
        "origin of the world coordinates.\nUse the following command line to correct the header of the "
        "faulty NIFTI volumes in a new folder:\n"
    ) in _build_warning_message(
        [tmp_path / "foo.nii", tmp_path / "bar.nii", tmp_path / "baz.nii"],
        [np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])],
        [0.2, 12.7, 26.6],
        tmp_path / "bids",
        "T1W",
    )


def test_validate_output_tsv_path(tmp_path):
    from clinica.iotools.utils.data_handling._merging import _validate_output_tsv_path

    (tmp_path / "foo.tsv").touch()

    assert _validate_output_tsv_path(tmp_path) == tmp_path / "merge.tsv"
    assert _validate_output_tsv_path(tmp_path / "foo.tsv") == tmp_path / "foo.tsv"
    assert (
        _validate_output_tsv_path(tmp_path / "bids" / "foo.tsv")
        == tmp_path / "bids" / "foo.tsv"
    )
    assert (tmp_path / "bids").exists()
    assert not (tmp_path / "bids" / "foo.tsv").exists()


def test_validate_output_tsv_path_error(tmp_path):
    from clinica.iotools.utils.data_handling._merging import _validate_output_tsv_path

    (tmp_path / "bar.txt").touch()

    with pytest.raises(
        TypeError,
        match="Output path extension must be tsv.",
    ):
        _validate_output_tsv_path(tmp_path / "bar.txt")


def test_get_groups(tmp_path):
    from clinica.iotools.utils.data_handling._missing import _get_groups

    assert _get_groups(tmp_path) == []
    (tmp_path / "groups").mkdir()
    assert _get_groups(tmp_path) == []
    for group in ("foo", "bar", "baz"):
        (tmp_path / "groups" / group).mkdir()
    assert sorted(_get_groups(tmp_path)) == ["bar", "baz", "foo"]


def create_bids_dataset(folder: Path, write_tsv_files: bool = False) -> None:
    from clinica.utils.testing_utils import build_bids_directory

    folder.mkdir()
    build_bids_directory(
        folder,
        {
            "sub-01": ["ses-M000", "ses-M006"],
            "sub-02": ["ses-M000"],
            "sub-03": ["ses-M000", "ses-M012", "ses-M024"],
        },
        modalities={
            "anat": {"T1w", "flair"},
            "pet": {"trc-18FFDG_pet"},
        },
        write_tsv_files=write_tsv_files,
    )


def test_find_mods_and_sess(tmp_path):
    from clinica.iotools.utils.data_handling._missing import _find_mods_and_sess

    create_bids_dataset(tmp_path / "bids")
    ses_mod = _find_mods_and_sess(tmp_path / "bids")
    assert len(ses_mod) == 3
    assert ses_mod["sessions"] == {"ses-M000", "ses-M006", "ses-M012", "ses-M024"}
    assert ses_mod["anat"] == {"t1w", "flair"}
    assert ses_mod["pet"] == {"pet_trc-18FFDG"}


@pytest.fixture
def expected_tsv_content() -> str:
    return (
        "participant_id\tsession_id\n"
        "sub-01\tses-M000\n"
        "sub-01\tses-M006\n"
        "sub-02\tses-M000\n"
        "sub-03\tses-M000\n"
        "sub-03\tses-M012\n"
        "sub-03\tses-M024\n"
    )


def test_create_subs_sess_list_as_text(tmp_path, expected_tsv_content):
    from clinica.iotools.utils.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    create_bids_dataset(tmp_path / "bids")
    assert (
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=False)
        == expected_tsv_content
    )


def test_create_subs_sess_list_as_text_using_tsv(tmp_path, expected_tsv_content):
    from clinica.iotools.utils.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)
    assert (
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=True)
        == expected_tsv_content
    )


def test_create_subs_sess_list_as_text_errors(tmp_path):
    from clinica.iotools.utils.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    with pytest.raises(IOError, match="Dataset empty or not BIDS/CAPS compliant."):
        _create_subs_sess_list_as_text(tmp_path, use_session_tsv=False)
    create_bids_dataset(tmp_path / "bids")
    with pytest.raises(
        ValueError, match="there is no session TSV file for subject sub-01"
    ):
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=True)


@pytest.mark.parametrize("use_tsv", [True, False])
def test_create_subs_sess_list(tmp_path, use_tsv, expected_tsv_content):
    from clinica.iotools.utils.data_handling import create_subs_sess_list

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)
    create_subs_sess_list(
        tmp_path / "bids",
        tmp_path / "output",
        file_name="foo.tsv",
        use_session_tsv=use_tsv,
    )
    assert (tmp_path / "output" / "foo.tsv").exists()
    assert (tmp_path / "output" / "foo.tsv").read_text() == expected_tsv_content


def test_write_list_of_files_errors(tmp_path):
    from clinica.iotools.utils.data_handling import write_list_of_files

    with pytest.raises(
        TypeError,
        match="`file_list` argument must be a list of paths. Instead <class 'int'> was provided.",
    ):
        write_list_of_files(10, tmp_path / "foo.txt")
    (tmp_path / "foo.txt").touch()
    with pytest.raises(
        IOError,
        match=f"Output file {tmp_path / 'foo.txt'} already exists.",
    ):
        write_list_of_files([], tmp_path / "foo.txt")


def test_write_list_of_files(tmp_path):
    from clinica.iotools.utils.data_handling import write_list_of_files

    file_list = [
        tmp_path / "foo.csv",
        tmp_path / "bar" / "baz.jpg",
        tmp_path / "boo.nii",
    ]
    write_list_of_files(file_list, tmp_path / "foo.txt")
    assert (tmp_path / "foo.txt").exists()
    assert (tmp_path / "foo.txt").read_text() == "\n".join(
        [str(f) for f in file_list]
    ) + "\n"
