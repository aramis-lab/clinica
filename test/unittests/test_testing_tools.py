import os
from os import PathLike
from pathlib import PurePath
from test.nonregression.testing_tools import (
    compare_folders_structures,
    compare_folders_with_hashes,
    create_list_hashes,
)
from typing import Callable

import nibabel as nib
import numpy as np
import pytest


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
    shape = (7, 7, 7, 7)
    img1 = nib.Nifti1Image(rng.random(shape), affine=np.eye(4))
    img2 = nib.Nifti1Image(rng.random(shape), affine=np.eye(4))
    img1.to_filename(str(tmp_path / "img1.nii"))
    img2.to_filename(str(tmp_path / "img2.nii"))
    assert not similarity_measure(tmp_path / "img1.nii", tmp_path / "img2.nii", 0.8)
    assert similarity_measure(tmp_path / "img1.nii", tmp_path / "img1.nii", 0.8)
    assert similarity_measure(tmp_path / "img2.nii", tmp_path / "img2.nii", 0.8)
    assert similarity_measure(tmp_path / "img1.nii", tmp_path / "img2.nii", 0.2)


def test_identical_subject_list(tmp_path: PurePath):
    from test.nonregression.testing_tools import identical_subject_list

    import pandas as pd

    df1 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M00", "ses-M06", "ses-M00"],
        }
    )
    df1.to_csv(tmp_path / "df1.tsv", sep="\t")
    df2 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-02"],
            "session_id": ["ses-M00", "ses-M00"],
        }
    )
    df2.to_csv(tmp_path / "df2.tsv", sep="\t")
    df3 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-03"],
            "session_id": ["ses-M00", "ses-M00"],
        }
    )
    df3.to_csv(tmp_path / "df3.tsv", sep="\t")
    df4 = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M00", "ses-M12", "ses-M00"],
        }
    )
    df4.to_csv(tmp_path / "df4.tsv", sep="\t")
    assert identical_subject_list(tmp_path / "df1.tsv", tmp_path / "df1.tsv")
    assert not identical_subject_list(tmp_path / "df1.tsv", tmp_path / "df2.tsv")
    assert not identical_subject_list(tmp_path / "df2.tsv", tmp_path / "df3.tsv")
    assert not identical_subject_list(tmp_path / "df1.tsv", tmp_path / "df4.tsv")


def test_same_missing_modality_tsv(tmp_path: PurePath):
    from test.nonregression.testing_tools import same_missing_modality_tsv

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
    assert same_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df1.tsv")
    assert same_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df2.tsv")
    assert not same_missing_modality_tsv(tmp_path / "df1.tsv", tmp_path / "df3.tsv")


def test_tree(tmp_path: PurePath):
    from test.nonregression.testing_tools import tree

    tree(str(tmp_path), str(tmp_path / "file_out.txt"))
    assert os.path.exists(tmp_path / "file_out.txt")
    with open(tmp_path / "file_out.txt", "r") as fp:
        content = fp.readlines()
    assert len(content) == 0
    os.makedirs(tmp_path / "subjects/sub-01/ses-M00/")
    tree(str(tmp_path), str(tmp_path / "file_out.txt"))
    with open(tmp_path / "file_out.txt", "r") as fp:
        content = fp.readlines()
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
    assert len(list_files_with_extensions(str(tmp_path), (".nii.gz", ".tsv"))) == 0
    assert list_files_with_extensions(str(tmp_path), (".txt",)) == [
        str(tmp_path / "foo.txt")
    ]
    assert set(list_files_with_extensions(str(tmp_path), (".txt", ".png"))) == {
        str(tmp_path / "foo.txt"),
        str(tmp_path / "bar.png"),
    }


def test_create_list_hashes(tmp_path: PurePath):
    from test.nonregression.testing_tools import create_list_hashes

    _create_files(tmp_path, ["foo.nii.gz", "bar.tsv", "baz.json", "foo.txt", "bar.png"])
    hashes = create_list_hashes(str(tmp_path))
    assert set(hashes.keys()) == {"/foo.nii.gz", "/bar.tsv", "/baz.json"}
    # change content of "baz.json" and check that hash is different
    with open(tmp_path / "baz.json", "w") as fp:
        fp.write("data")
    hashes2 = create_list_hashes(str(tmp_path))
    assert hashes["/baz.json"] != hashes2["/baz.json"]
    for key in ["/foo.nii.gz", "/bar.tsv"]:
        assert hashes[key] == hashes2[key]


@pytest.mark.parametrize(
    "compare_func", [compare_folders_structures, compare_folders_with_hashes]
)
def test_compare_folders_structures(
    tmp_path: PurePath, compare_func: Callable[[PathLike, PathLike], None]
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
    hashes = create_list_hashes(str(tmp_path))
    with open(tmp_path / "hashes.pl", "wb") as fp:
        pickle.dump(hashes, fp)

    # Basic check that structures and files match the ref
    compare_func(str(tmp_path), str(tmp_path / "hashes.pl"))

    # Change the content of a file should not change the structure
    # but only the hashes
    with open(str(tmp_path / "sub-02/baz.json"), "w") as fp:
        fp.write("data")
    if compare_func.__name__ == "compare_folders_structures":
        compare_func(str(tmp_path), str(tmp_path / "hashes.pl"))
    else:
        with pytest.raises(
            ValueError, match="/sub-02/baz.json does not match the reference file !"
        ):
            compare_func(str(tmp_path), str(tmp_path / "hashes.pl"))

    # Delete all files and folders for subject 2
    # This should change the structure compared to the reference
    shutil.rmtree(tmp_path / "sub-02")
    with pytest.raises(ValueError, match="/sub-02/bar.tsv not found !"):
        compare_func(str(tmp_path), str(tmp_path / "hashes.pl"))
