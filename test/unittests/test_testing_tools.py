import os
from os import PathLike
from pathlib import Path, PurePath
from test.nonregression.testing_tools import (
    Level,
    _load_participants_tsv,
    _load_scans_tsv,
    _load_sessions_tsv,
    compare_bids_tsv,
    compare_folders_structures,
    compare_folders_with_hashes,
    compare_niftis,
    create_list_hashes,
)
from typing import Callable, Optional

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


def test_likeliness_measure(tmp_path: PurePath):
    from test.nonregression.testing_tools import likeliness_measure

    rng = np.random.RandomState(42)
    img1 = nib.Nifti1Image(rng.random((2, 2, 2, 2)), affine=np.eye(4))
    img2 = nib.Nifti1Image(rng.random((2, 2, 2, 2)), affine=np.eye(4))
    img1.to_filename(str(tmp_path / "img1.nii"))
    img2.to_filename(str(tmp_path / "img2.nii"))
    assert not likeliness_measure(
        tmp_path / "img1.nii", tmp_path / "img2.nii", (1e-4, 1e-4), (1e-4, 1e-4)
    )
    assert likeliness_measure(
        tmp_path / "img1.nii", tmp_path / "img1.nii", (1e-4, 1e-4), (1e-4, 1e-4)
    )


def test_similarity_measure(tmp_path: PurePath):
    from test.nonregression.testing_tools import similarity_measure

    rng = np.random.RandomState(42)
    shape = (16, 16, 16, 16)

    img1 = nib.Nifti1Image(rng.random(shape), affine=np.eye(4))
    file1 = tmp_path / "img1.nii"
    img1.to_filename(os.fspath(file1))

    img2 = nib.Nifti1Image(rng.random(shape), affine=np.eye(4))
    file2 = tmp_path / "img2.nii"
    img2.to_filename(os.fspath(file2))

    assert similarity_measure(file1, file1, 0.8)
    assert similarity_measure(file2, file2, 0.8)
    assert not similarity_measure(file1, file2, 0.2)


def test_identical_subject_list(tmp_path: PurePath):
    from test.nonregression.testing_tools import compare_subject_session_tsv

    import pandas as pd

    df1 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M00", "ses-M006", "ses-M000"],
        }
    )
    df1.to_csv(tmp_path / "df1.tsv", sep="\t")
    df2 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-02"],
            "session_id": ["ses-M000", "ses-M000"],
        }
    )
    df2.to_csv(tmp_path / "df2.tsv", sep="\t")
    df3 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-03"],
            "session_id": ["ses-M000", "ses-M000"],
        }
    )
    df3.to_csv(tmp_path / "df3.tsv", sep="\t")
    df4 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M000", "ses-M012", "ses-M000"],
        }
    )
    df4.to_csv(tmp_path / "df4.tsv", sep="\t")

    assert compare_subject_session_tsv(tmp_path / "df1.tsv", tmp_path / "df1.tsv")
    assert not compare_subject_session_tsv(tmp_path / "df1.tsv", tmp_path / "df2.tsv")
    assert not compare_subject_session_tsv(tmp_path / "df2.tsv", tmp_path / "df3.tsv")
    assert not compare_subject_session_tsv(tmp_path / "df1.tsv", tmp_path / "df4.tsv")


def test_same_missing_modality_tsv(tmp_path: PurePath):
    from test.nonregression.testing_tools import compare_missing_modality_tsv

    import pandas as pd

    df1 = pd.DataFrame(
        {
            "participant_id": ["sub-02", "sub-01", "sub-03", "sub-02"],
            "pet_trc-18FAV45": ["b", "b", "a", "a"],
            "pet_trc-18FFDG": ["x", "x", "x", "x"],
            "t1w": ["a", "b", "a", "b"],
            "func_task-rest": ["a", "b", "b", "a"],
        }
    )
    df1.to_csv(tmp_path / "df1.tsv", sep="\t")
    df2 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-02", "sub-03", "sub-02"],
            "pet_trc-18FAV45": ["b", "b", "a", "a"],
            "pet_trc-18FFDG": ["x", "x", "x", "x"],
            "t1w": ["b", "a", "a", "b"],
            "func_task-rest": ["b", "a", "b", "a"],
        }
    )
    df2.to_csv(tmp_path / "df2.tsv", sep="\t")
    df3 = pd.DataFrame(
        {
            "participant_id": ["sub-03", "sub-02", "sub-01", "sub-02"],
            "pet_trc-18FAV45": ["b", "b", "a", "a"],
            "pet_trc-18FFDG": ["x", "x", "x", "x"],
            "t1w": ["b", "a", "a", "b"],
            "func_task-rest": ["b", "a", "b", "a"],
        }
    )
    df3.to_csv(tmp_path / "df3.tsv", sep="\t")

    assert compare_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df1.tsv")
    assert compare_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df2.tsv")
    assert not compare_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df3.tsv")


def test_tree(tmp_path):
    from test.nonregression.testing_tools import tree

    tree(tmp_path, tmp_path / "file_out.txt")
    assert (tmp_path / "file_out.txt").exists()

    content = (tmp_path / "file_out.txt").read_text()
    assert len(content) == 0

    (tmp_path / "subjects/sub-01/ses-M000/").mkdir(parents=True)
    tree(tmp_path, tmp_path / "file_out.txt")

    content = (tmp_path / "file_out.txt").read_text()
    assert "".join(content) == (
        "    + file_out.txt\n    + subjects\n        + sub-01\n"
    )


def _create_files(folder, list_of_filenames):
    for f in list_of_filenames:
        with open(folder / f, "w") as _:
            pass


def test_list_files_with_extensions(tmp_path: PurePath) -> None:
    from test.nonregression.testing_tools import list_files_with_extensions

    _create_files(tmp_path, ["foo.txt", "bar.png"])
    assert len(list_files_with_extensions(tmp_path, (".nii.gz", ".tsv"))) == 0
    assert list_files_with_extensions(tmp_path, (".txt",)) == [
        str(tmp_path / "foo.txt")
    ]
    assert set(list_files_with_extensions(tmp_path, (".txt", ".png"))) == {
        str(tmp_path / "foo.txt"),
        str(tmp_path / "bar.png"),
    }


def test_create_list_hashes(tmp_path):
    from test.nonregression.testing_tools import create_list_hashes

    _create_files(tmp_path, ["foo.nii.gz", "bar.tsv", "baz.json", "foo.txt", "bar.png"])
    hashes = create_list_hashes(tmp_path)
    assert set(hashes.keys()) == {"/foo.nii.gz", "/bar.tsv", "/baz.json"}

    # change content of "baz.json" and check that hash is different
    (tmp_path / "baz.json").write_text("data")
    hashes2 = create_list_hashes(tmp_path)
    assert hashes["/baz.json"] != hashes2["/baz.json"]
    for key in ["/foo.nii.gz", "/bar.tsv"]:
        assert hashes[key] == hashes2[key]


@pytest.mark.parametrize(
    "compare_func", [compare_folders_structures, compare_folders_with_hashes]
)
def test_compare_folders_structures(
    tmp_path, compare_func: Callable[[PathLike, PathLike], None]
):
    import pickle
    import shutil

    # Setup the test data
    for subject in ["sub-01", "sub-02"]:
        os.makedirs(tmp_path / subject)
        _create_files(
            tmp_path / subject,
            ["foo.nii.gz", "bar.tsv", "baz.json", "foo.txt", "bar.png"],
        )
    hashes = create_list_hashes(tmp_path)
    with open(tmp_path / "hashes.pl", "wb") as fp:
        pickle.dump(hashes, fp)

    # Basic check that structures and files match the ref
    compare_func(tmp_path, tmp_path / "hashes.pl")

    # Change the content of a file should not change the structure
    # but only the hashes
    (tmp_path / "sub-02/baz.json").write_text("data")
    if compare_func.__name__ == "compare_folders_structures":
        compare_func(tmp_path, tmp_path / "hashes.pl")
    else:
        with pytest.raises(
            ValueError, match="/sub-02/baz.json does not match the reference file !"
        ):
            compare_func(tmp_path, tmp_path / "hashes.pl")

    # Delete all files and folders for subject 2
    # This should change the structure compared to the reference
    shutil.rmtree(tmp_path / "sub-02")
    with pytest.raises(ValueError, match="/sub-02/bar.tsv not found !"):
        compare_func(tmp_path, tmp_path / "hashes.pl")


def build_bids_tsv(tmp_path: Path) -> Path:
    bids_path = tmp_path / "BIDS"
    bids_path.mkdir()
    prpc = pd.DataFrame({"participant_id": ["sub-002", "sub-001"], "age": [20, 26]})
    prpc.to_csv(bids_path / "participants.tsv", sep="\t", index=False)
    sub_path = bids_path / "sub-001"
    sub_path.mkdir()
    sess = pd.DataFrame({"session_id": ["ses-M012", "ses-M006"], "age": [20, 20]})
    sess.to_csv(sub_path / "sub-001_sessions.tsv", sep="\t", index=False)
    ses_path = sub_path / "ses-M016"
    ses_path.mkdir()
    scans = pd.DataFrame(
        {
            "filename": ["pet/foo.json", "anat/foo.json"],
            "acq_time": ["00:00:00", "00:00:00"],
        }
    )
    scans.to_csv(ses_path / "sub-001_ses-M016_scans.tsv", sep="\t", index=False)
    return bids_path


def test_loader_participants(tmp_path):
    bids_path = build_bids_tsv(tmp_path)

    assert_frame_equal(
        pd.DataFrame({"participant_id": ["sub-001", "sub-002"], "age": [26, 20]}),
        _load_participants_tsv(bids_path, Path("")),
    )


def test_loader_sessions(tmp_path):
    bids_path = build_bids_tsv(tmp_path)

    assert_frame_equal(
        pd.DataFrame({"session_id": ["ses-M006", "ses-M012"], "age": [20, 20]}),
        _load_sessions_tsv(bids_path, bids_path / "sub-001" / "sub-001_sessions.tsv"),
    )


def test_loader_scans(tmp_path):
    bids_path = build_bids_tsv(tmp_path)

    assert_frame_equal(
        pd.DataFrame(
            {
                "filename": ["anat/foo.json", "pet/foo.json"],
                "acq_time": ["00:00:00", "00:00:00"],
            }
        ),
        _load_scans_tsv(
            bids_path, bids_path / "sub-001" / "ses-M016" / "sub-001_ses-M016_scans.tsv"
        ),
    )


@pytest.mark.parametrize(
    "level, expected",
    [
        ("participants", _load_participants_tsv),
        (Level.SESSIONS, _load_sessions_tsv),
        (Level.SCANS, _load_scans_tsv),
    ],
)
def test_loader_factory(level, expected):
    from test.nonregression.testing_tools import _loader_factory

    assert expected == _loader_factory(level)


def test_loader_factory_error():
    from test.nonregression.testing_tools import _loader_factory

    with pytest.raises(ValueError):
        _loader_factory("foo")


def test_compare_bids_tsv_success(tmp_path):
    bids_path = build_bids_tsv(tmp_path)
    compare_bids_tsv(bids_path, bids_path)


@pytest.mark.parametrize(
    "modified_frame, frame_path, error_message",
    [
        (
            pd.DataFrame({"participant_id": ["sub-001"], "age": [26]}),
            "participants.tsv",
            "participants.tsv shape mismatch",
        ),
        (
            pd.DataFrame({"session_id": ["ses-M012", "ses-M006"], "age": [20, 25]}),
            "sub-001/sub-001_sessions.tsv",
            r"sub-001_sessions.tsv.* values are different",
        ),
        (
            pd.DataFrame(
                {
                    "filename": ["pet/foo.nii.gz", "anat/foo.json"],
                    "acq_time": ["00:00:00", "00:00:00"],
                }
            ),
            "sub-001/ses-M016/sub-001_ses-M016_scans.tsv",
            r"sub-001_ses-M016_scans.tsv.* values are different",
        ),
    ],
)
def test_compare_bids_tsv_error(tmp_path, modified_frame, frame_path, error_message):
    from shutil import copytree

    bids_path = build_bids_tsv(tmp_path)
    copy = tmp_path / "BIDS_copy"
    copytree(bids_path, copy)
    modified_frame.to_csv(copy / frame_path, sep="\t", index=False)
    with pytest.raises(AssertionError, match=error_message):
        compare_bids_tsv(bids_path, copy)


def _create_nifti_in_dir(dir_path: Path, rds) -> Path:
    dir_path.mkdir()
    nib.Nifti1Image(rds.random((16, 16, 16, 16)), affine=np.eye(4)).to_filename(
        dir_path / "nifti.nii.gz"
    )
    return dir_path


def test_compare_niftis_success(tmp_path):
    compare_niftis(
        _create_nifti_in_dir(tmp_path / "ref", np.random.RandomState(42)),
        _create_nifti_in_dir(tmp_path / "out", np.random.RandomState(42)),
    )


def test_compare_niftis_error(tmp_path):
    with pytest.raises(
        AssertionError,
        match="Following images do not meet the similarity criteria : \n\nnifti.nii.gz",
    ):
        rds = np.random.RandomState(42)
        compare_niftis(
            _create_nifti_in_dir(tmp_path / "ref", rds),
            _create_nifti_in_dir(tmp_path / "out", rds),
        )
