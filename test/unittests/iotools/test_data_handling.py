from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_equal
from pandas.testing import assert_frame_equal


def test_scale_coordinates_by_pixdim():
    """Test function `_scale_coordinates_by_pixdim`."""
    from nibabel.nifti1 import Nifti1Header

    from clinica.iotools.data_handling._centering import (
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
    from clinica.iotools.data_handling import (
        check_relative_volume_location_in_world_coordinate_system,
    )

    mocker.patch(
        "clinica.iotools.data_handling._centering._get_problematic_pairs_with_l2_norm",
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
        )


def test_get_problematic_pairs_empty(tmp_path, mocker):
    from clinica.iotools.data_handling._centering import (
        _get_problematic_pairs_with_l2_norm,
    )

    mocker.patch(
        "clinica.iotools.data_handling._centering._get_world_coordinate_of_center",
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
    from clinica.iotools.data_handling._centering import (
        _get_problematic_pairs_with_l2_norm,
    )

    mocker.patch(
        "clinica.iotools.data_handling._centering._compute_l2_norm",
        return_value=[81.0, 79.9],
    )
    assert _get_problematic_pairs_with_l2_norm(
        [
            (tmp_path / "foo.nii.gz", tmp_path / "baz.nii.gz"),
            (tmp_path / "bar.nii.gz", tmp_path / "foo_bar.nii.gz"),
        ]
    ) == [("foo.nii.gz", "baz.nii.gz", 81.0)]


def test_build_warning_message(tmp_path):
    from clinica.iotools.data_handling._centering import _build_warning_message

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
    from clinica.iotools.data_handling._merging import _validate_output_tsv_path

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
    from clinica.iotools.data_handling._merging import _validate_output_tsv_path

    (tmp_path / "bar.txt").touch()

    with pytest.raises(
        TypeError,
        match="Output path extension must be tsv.",
    ):
        _validate_output_tsv_path(tmp_path / "bar.txt")


def test_get_groups(tmp_path):
    from clinica.iotools.data_handling._missing import _get_groups

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
            "anat": ("T1w", "flair"),
            "pet": ("trc-18FFDG_pet",),
        },
        write_tsv_files=write_tsv_files,
    )


def test_find_mods_and_sess(tmp_path):
    from clinica.iotools.data_handling._missing import _find_mods_and_sess

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


def test_create_subs_sess_list_as_text(tmp_path, expected_tsv_content: str):
    from clinica.iotools.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    create_bids_dataset(tmp_path / "bids")

    assert (
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=False)
        == expected_tsv_content
    )


def test_create_subs_sess_list_as_text_using_tsv(tmp_path, expected_tsv_content: str):
    from clinica.iotools.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)

    assert (
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=True)
        == expected_tsv_content
    )


def test_create_subs_sess_list_as_text_errors(tmp_path):
    from clinica.iotools.data_handling._files import (
        _create_subs_sess_list_as_text,
    )

    with pytest.raises(IOError, match="Dataset empty or not BIDS/CAPS compliant."):
        _create_subs_sess_list_as_text(tmp_path, use_session_tsv=False)

    create_bids_dataset(tmp_path / "bids")

    with pytest.raises(
        ValueError,
        match="there is no session TSV file for subject sub-01",
    ):
        _create_subs_sess_list_as_text(tmp_path / "bids", use_session_tsv=True)


@pytest.mark.parametrize("use_tsv", [True, False])
def test_create_subs_sess_list(tmp_path, use_tsv: bool, expected_tsv_content: str):
    from clinica.iotools.data_handling import create_subs_sess_list

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
    from clinica.iotools.data_handling import write_list_of_files

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
    from clinica.iotools.data_handling import write_list_of_files

    file_list = [
        tmp_path / "foo.csv",
        tmp_path / "bar" / "baz.jpg",
        tmp_path / "boo.nii",
    ]
    write_list_of_files(file_list, tmp_path / "foo.txt")

    assert (tmp_path / "foo.txt").exists()
    assert (tmp_path / "foo.txt").read_text() == "\n".join([str(f) for f in file_list])


def test_get_participants_and_subjects_sessions_df(tmp_path):
    from clinica.iotools.data_handling._merging import (
        _get_participants_and_subjects_sessions_df,
    )

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)
    participants, sessions = _get_participants_and_subjects_sessions_df(
        tmp_path / "bids"
    )
    assert_frame_equal(
        participants, pd.DataFrame({"participant_id": ["sub-01", "sub-02", "sub-03"]})
    )

    assert len(sessions) == 6


@pytest.mark.parametrize("ignore_sessions", (True, False))
def test_create_merge_file_from_bids_ignore_sessions_and_scans(
    tmp_path, ignore_sessions: bool
):
    from clinica.iotools.data_handling._merging import (
        _create_merge_file_from_bids,
        _get_participants_and_subjects_sessions_df,
    )

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)
    participants, sessions = _get_participants_and_subjects_sessions_df(
        tmp_path / "bids"
    )
    df = _create_merge_file_from_bids(
        tmp_path / "bids",
        sessions,
        participants,
        ignore_sessions_files=ignore_sessions,
        ignore_scan_files=True,
    )
    expected = pd.DataFrame(
        {
            "participant_id": [
                "sub-01",
                "sub-01",
                "sub-02",
                "sub-03",
                "sub-03",
                "sub-03",
            ],
            "session_id": [
                "ses-M000",
                "ses-M006",
                "ses-M000",
                "ses-M000",
                "ses-M012",
                "ses-M024",
            ],
        }
    )

    assert_frame_equal(df.reset_index(drop=True), expected.reset_index(drop=True))


def test_create_merge_file_from_bids(tmp_path):
    from clinica.iotools.data_handling._merging import (
        _create_merge_file_from_bids,
        _get_participants_and_subjects_sessions_df,
    )

    create_bids_dataset(tmp_path / "bids", write_tsv_files=True)
    participants, sessions = _get_participants_and_subjects_sessions_df(
        tmp_path / "bids"
    )
    df = _create_merge_file_from_bids(
        tmp_path / "bids",
        sessions,
        participants,
        ignore_sessions_files=False,
        ignore_scan_files=False,
    )
    expected = pd.DataFrame(
        {
            "participant_id": [
                "sub-01",
                "sub-01",
                "sub-02",
                "sub-03",
                "sub-03",
                "sub-03",
            ],
            "session_id": [
                "ses-M000",
                "ses-M006",
                "ses-M000",
                "ses-M000",
                "ses-M012",
                "ses-M024",
            ],
            "T1w_scan_id": ["670488", "777573", "234054", "107474", "935519", "619177"],
            "flair_scan_id": [
                "116740",
                "288390",
                "146317",
                "709571",
                "571859",
                "442418",
            ],
            "trc-18FFDG_pet_scan_id": [
                "26226",
                "256788",
                "772247",
                "776647",
                "91162",
                "33327",
            ],
        }
    )
    assert_frame_equal(
        df.reset_index(drop=True),
        expected.reset_index(drop=True),
        check_dtype=False,
        check_column_type=False,
    )


def test_post_process_merge_file_from_bids_column_reordering():
    from clinica.iotools.data_handling._merging import (
        _post_process_merge_file_from_bids,
    )

    df = pd.DataFrame(columns=["foo", "session_id", "bar", "participant_id", "baz"])

    assert_frame_equal(
        _post_process_merge_file_from_bids(df),
        pd.DataFrame(columns=["participant_id", "session_id", "foo", "bar", "baz"]),
    )


def test_post_process_merge_file_from_bids_number_rounding():
    from clinica.iotools.data_handling._merging import (
        _post_process_merge_file_from_bids,
    )

    df = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M006"],
            "participant_id": ["sub-01", "sub-01"],
            "foo": [1.0123456789, 100.987654321],
        }
    )
    expected = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01"],
            "session_id": ["ses-M000", "ses-M006"],
            "foo": [1.012346, 100.987654],
        }
    )
    assert_frame_equal(_post_process_merge_file_from_bids(df), expected)


def test_add_data_to_merge_file_from_caps_wrong_handler(tmp_path):
    from clinica.iotools.data_handling._merging import (
        _add_data_to_merge_file_from_caps,
    )

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid PipelineNameForMetricExtraction",
    ):
        _add_data_to_merge_file_from_caps(tmp_path, pd.DataFrame(), pipelines=["foo"])


def test_add_data_to_merge_file_from_caps_empty_pipeline(tmp_path):
    from clinica.iotools.data_handling._merging import (
        _add_data_to_merge_file_from_caps,
    )

    df = pd.DataFrame(
        {
            "participant_id": [
                "sub-01",
                "sub-01",
                "sub-02",
            ],
            "session_id": [
                "ses-M000",
                "ses-M006",
                "ses-M000",
            ],
        }
    )
    (tmp_path / "groups").mkdir()
    with pytest.raises(
        FileNotFoundError,
        match=(
            "No outputs were found for any pipeline in the CAPS folder. "
            "The output only contains BIDS information."
        ),
    ):
        _add_data_to_merge_file_from_caps(tmp_path, df, pipelines=[])


@pytest.mark.parametrize(
    "input_list,expected",
    [
        (
            ["ses-M000", "ses-M006", "ses-M012", "ses-M024", "ses-M048", "ses-M003"],
            ["ses-M000", "ses-M003", "ses-M006", "ses-M012", "ses-M024", "ses-M048"],
        ),
        (
            ["ses-M00", "ses-M06", "ses-M12", "ses-M24", "ses-M48", "ses-M03"],
            ["ses-M00", "ses-M03", "ses-M06", "ses-M12", "ses-M24", "ses-M48"],
        ),
        (
            ["ses-M0", "ses-M6", "ses-M12", "ses-M24", "ses-M48", "ses-M3"],
            ["ses-M0", "ses-M3", "ses-M6", "ses-M12", "ses-M24", "ses-M48"],
        ),
    ],
)
def test_sort_session_list(input_list, expected):
    """Test function `sort_session_list`."""
    from clinica.iotools.data_handling._missing_modality_tracker import (
        _sort_session_list,
    )

    assert _sort_session_list(input_list) == expected


@pytest.mark.parametrize(
    "sessions,modalities,expected",
    [
        ([], None, {}),
        (
            ["ses-M000", "ses-M006"],
            None,
            {
                k: {
                    "session": 0,
                    "dwi": 0,
                    "func": 0,
                    "fieldmap": 0,
                    "flair": 0,
                    "t1w": 0,
                }
                for k in ["ses-M000", "ses-M006"]
            },
        ),
        ([], ["fmri"], {}),
        (
            ["ses-M000", "ses-M006"],
            ["fmri", "dwi"],
            {k: {"session": 0, "dwi": 0, "fmri": 0} for k in ["ses-M000", "ses-M006"]},
        ),
    ],
)
def test_missing_mods_tracker_instantiation(sessions, modalities, expected):
    from clinica.iotools.data_handling._missing_modality_tracker import (
        MissingModsTracker,
    )

    tracker = MissingModsTracker(sessions, modalities)

    assert tracker.missing == expected
    assert tracker.ses == sessions
    assert tracker.get_missing_list() == expected


def test_missing_mods_tracker_add_errors():
    from clinica.iotools.data_handling._missing_modality_tracker import (
        MissingModsTracker,
    )

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])

    with pytest.raises(
        ValueError,
        match="Session foo was not provided to the MissingModsTracker constructor.",
    ):
        tracker.add_missing_mod("foo", "fmri")
        tracker.increase_missing_ses("foo")

    with pytest.raises(
        ValueError,
        match="Modality foo is not tracked by this instance of MissingModsTracker.",
    ):
        tracker.add_missing_mod("ses-M006", "foo")


def test_missing_mods_tracker_add():
    from clinica.iotools.data_handling._missing_modality_tracker import (
        MissingModsTracker,
    )

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])
    tracker.increase_missing_ses("ses-M000")
    assert tracker.missing["ses-M000"]["session"] == 1
    assert tracker.missing["ses-M006"]["session"] == 0
    tracker.add_missing_mod("ses-M006", "dwi")
    assert tracker.missing["ses-M006"]["dwi"] == 1
    assert tracker.missing["ses-M000"]["dwi"] == 0
    tracker.add_missing_mod("ses-M000", "dwi")
    assert tracker.missing["ses-M006"]["dwi"] == 1
    assert tracker.missing["ses-M000"]["dwi"] == 1


def test_compute_statistics():
    from clinica.iotools.data_handling._missing_modality_tracker import (
        MissingModsTracker,
        _compute_statistics,
    )

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])
    tracker.increase_missing_ses("ses-M000")
    tracker.add_missing_mod("ses-M006", "dwi")

    assert _compute_statistics(2, ["ses-M000", "ses-M006"], tracker) == (
        "**********************************************\n"
        "Number of subjects converted: 2\n"
        "Sessions available: ['ses-M000', 'ses-M006']\n\n"
        "Number of sessions ses-M000 found: 1 (50.0%)\n\n"
        "Number of sessions ses-M006 found: 2 (100.0%)\n\n"
        "**********************************************\n\n"
        "Number of missing modalities for each session:\n\n"
        "ses-M000\nfmri: 0 (0.0%) \ndwi: 0 (0.0%) \n\n"
        "ses-M006\nfmri: 0 (0.0%) \ndwi: 1 (50.0%) \n"
    )


def test_compute_table():
    from clinica.iotools.data_handling._missing_modality_tracker import _compute_table

    mods = {
        "foo": {
            "missing": 42,
            "n/a": 16,
            "bar": 23,
        },
        "bar": {
            "missing": 0,
            "baz": 7,
        },
    }

    assert _compute_table(mods) == (
        "\tbar\t| baz\t| missing\t| n/a\n------------------------------------------------\n"
        "foo\t23\t| 0\t| 42\t| 16\nbar\t0\t| 7\t| 0\t| 0\n"
    )


def test_compute_longitudinal_analysis(tmp_path):
    from clinica.iotools.data_handling._missing_modality_tracker import (
        _compute_longitudinal_analysis,
    )

    sessions = ["ses-M000", "ses-M006"]
    for subject in ("sub-01", "sub-03"):
        (tmp_path / "bids" / subject).mkdir(parents=True)
    for ses in sessions:
        df = pd.DataFrame(
            [["sub-01", 1, 0, 1], ["sub-03", 0, 1, 1]],
            columns=["participant_id", "foo", "bar", "baz"],
        )
        df.to_csv(tmp_path / f"missing_mods_{ses}.tsv", sep="\t", index=False)
    df = pd.DataFrame(
        [
            ["ses-M000", 18, "foooo", "CN", "ADNI"],
            ["ses-M006", 19, "foooo", "AD", "ADNI"],
        ],
        columns=["session_id", "age", "foobarbaz", "diagnosis", "study"],
    )
    df.to_csv(
        tmp_path / "bids" / "sub-01" / "sub-01_sessions.tsv", sep="\t", index=False
    )
    df = pd.DataFrame(
        [["ses-M000", 66, "bar", 42, "ADNI"], ["ses-M012", 67, "baz", 69, "ADNI"]],
        columns=["session_id", "age", "bar", "diagnosis", "study"],
    )
    df.to_csv(
        tmp_path / "bids" / "sub-03" / "sub-03_sessions.tsv", sep="\t", index=False
    )

    assert _compute_longitudinal_analysis(
        tmp_path / "bids", tmp_path, sessions, "missing_mods_"
    ) == (
        "**********************************************\n\n"
        "Number of present diagnoses and modalities "
        "for each session:\nses-M000\n\tCN\t| n/a\n"
        "--------------------------------\nfoo\t1\t| 0\nbar\t0\t| "
        "1\nbaz\t1\t| 1\n\n\nses-M006\n\tAD\t| missing\n"
        "--------------------------------\nfoo\t1\t| "
        "0\nbar\t0\t| 1\nbaz\t1\t| 1\n\n\n"
    )
