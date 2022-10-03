from pathlib import Path

import pytest


def rmtree(f: Path):
    if f.is_file():
        f.unlink()
    else:
        for child in f.iterdir():
            rmtree(child)
        f.rmdir()


def test_insensitive_glob(tmp_path):
    from clinica.utils.inputs import insensitive_glob

    files = [
        "foo.py",
        "Bar.txt",
        "bAZ.py",
        "Fooo.PY",
        "folder_1",
        "folder_2",
        "folder_1/foo1.py",
        "folder_2/BaR2.PY",
        "folder_1/baz.txt",
    ]
    for file in files:
        d = tmp_path / file
        d.mkdir()
    python_files = insensitive_glob(str(tmp_path / "*.py"))
    assert set([Path(f).name for f in python_files]) == {"foo.py", "bAZ.py", "Fooo.PY"}
    text_files = insensitive_glob(str(tmp_path / "*.txt"))
    assert set([Path(f).name for f in text_files]) == {"Bar.txt"}
    assert len(insensitive_glob(str(tmp_path / "*.json"))) == 0
    all_python_files = insensitive_glob(str(tmp_path / "**/*.py"), recursive=True)
    assert set([Path(f).name for f in all_python_files]) == {
        "foo.py",
        "bAZ.py",
        "Fooo.PY",
        "foo1.py",
        "BaR2.PY",
    }


def test_determine_caps_or_bids(tmp_path):
    from clinica.utils.inputs import determine_caps_or_bids

    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "subjects" / "sub-01").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / "foo.txt").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    rmtree(tmp_path / "subjects")
    (tmp_path / "sub-01").mkdir()
    (tmp_path / "sub-01" / "foo.txt").mkdir()
    assert determine_caps_or_bids(tmp_path)


@pytest.mark.parametrize("folder_type", ["BIDS", "CAPS"])
def test_common_checks(folder_type):
    from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError
    from clinica.utils.inputs import _common_checks

    with pytest.raises(
        ValueError,
        match="Argument you provided to ",
    ):
        _common_checks(1, folder_type)

    error = ClinicaBIDSError if folder_type == "BIDS" else ClinicaCAPSError

    with pytest.raises(
        error,
        match=f"The {folder_type} directory you gave is not a folder.",
    ):
        _common_checks("fooooo", folder_type)


def test_check_bids_folder(tmp_path):
    from clinica.utils.exceptions import ClinicaBIDSError
    from clinica.utils.inputs import check_bids_folder

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="The BIDS directory",
    ):
        check_bids_folder(tmp_path)
    rmtree(tmp_path / "subjects")
    (tmp_path / "data").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="The BIDS directory you provided is empty.",
    ):
        check_bids_folder(tmp_path / "data")
    (tmp_path / "data" / "foo").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="Your BIDS directory does not contains a single folder whose name",
    ):
        check_bids_folder(tmp_path / "data")
