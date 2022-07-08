import os
import pytest
import numpy as np
import nibabel as nib


def test_likeliness_measure():
    pass


def test_similarity_measure(tmp_path):
    from test.nonregression.testing_tools import similarity_measure
    rng = np.random.RandomState(42)
    img1 = nib.Nifti1Image(rng.random((2, 2, 2, 2)), affine=np.eye(4))
    img2 = nib.Nifti1Image(rng.random((2, 2, 2, 2)), affine=np.eye(4))
    img1.to_filename(str(tmp_path / 'img1.nii'))
    img2.to_filename(str(tmp_path / 'img2.nii'))
    assert not similarity_measure(tmp_path / 'img1.nii', tmp_path / 'img2.nii', .8)
    assert similarity_measure(tmp_path / 'img1.nii', tmp_path / 'img1.nii', .8)
    assert similarity_measure(tmp_path / 'img2.nii', tmp_path / 'img2.nii', .8)
    assert similarity_measure(tmp_path / 'img1.nii', tmp_path / 'img2.nii', .2)


def test_identical_subject_list():
    pass


def test_same_missing_modality_tsv():
    pass


def test_tree(tmp_path):
    from test.nonregression.testing_tools import tree
    tree(tmp_path, tmp_path / "file_out.txt")
    assert os.path.exists(tmp_path / "file_out.txt")
    with open(tmp_path / "file_out.txt", "r") as fp:
        content = fp.readlines()
    assert len(content) == 0
    os.makedirs(tmp_path / "subjects/sub-01/ses-M00/")
    tree(tmp_path, tmp_path / "file_out.txt")
    with open(tmp_path / "file_out.txt", "r") as fp:
        content = fp.readlines()
    assert (
        "".join(content) == (
            '    + file_out.txt\n    + subjects\n        + sub-01\n'
            )
    )


def _create_files(folder, list_of_filenames):
    for f in list_of_filenames:
        with open(folder / f, "w") as fp:
            pass


def test_list_files_with_extensions(tmp_path):
    from test.nonregression.testing_tools import list_files_with_extensions
    _create_files(tmp_path, ["foo.txt", "bar.png"])
    assert len(list_files_with_extensions(tmp_path, (".nii.gz", ".tsv"))) == 0
    assert list_files_with_extensions(tmp_path, (".txt")) == [str(tmp_path / "foo.txt")]
    assert (
        set(list_files_with_extensions(tmp_path, (".txt", ".png")))
        == set([str(tmp_path / "foo.txt"), str(tmp_path / "bar.png")])
    )


def test_create_list_hashes(tmp_path):
    from test.nonregression.testing_tools import create_list_hashes
    _create_files(tmp_path, ["foo.nii.gz", "bar.tsv", "baz.json", "foo.txt", "bar.png"])
    hashes = create_list_hashes(tmp_path)
    assert set(hashes.keys()) == set(["/foo.nii.gz", "/bar.tsv", "/baz.json"])
    # change content of "baz.json" and check that hash is different
    with open(tmp_path / "baz.json", "w") as fp:
        fp.write("data")
    hashes2 = create_list_hashes(tmp_path)
    assert hashes["/baz.json"] != hashes2["/baz.json"]
    for key in ["/foo.nii.gz", "/bar.tsv"]:
        assert hashes[key] == hashes2[key]


def test_compare_folders_structures(tmp_path):
    import pickle
    import shutil
    from test.nonregression.testing_tools import (
            compare_folders_structures, create_list_hashes
    )
    for subject in ["sub-01", "sub-02"]:
        os.makedirs(tmp_path / subject)
        _create_files(
                tmp_path / subject,
                ["foo.nii.gz", "bar.tsv", "baz.json", "foo.txt", "bar.png"]
        )
    hashes = create_list_hashes(tmp_path)
    with open(tmp_path / "hashes.pl", "wb") as fp:
        pickle.dump(hashes, fp)
    compare_folders_structures(tmp_path, tmp_path / "hashes.pl")
    shutil.rmtree(tmp_path / "sub-02")
    with pytest.raises(
        ValueError,
        match="/sub-02/bar.tsv not found !"
    ):
        compare_folders_structures(tmp_path, tmp_path / "hashes.pl")
