import pytest


def test_listdir_nohidden(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import _listdir_nohidden

    (tmp_path / "file").touch()
    (tmp_path / ".hidden_file").touch()
    (tmp_path / "dir").mkdir()
    (tmp_path / ".hidden_dir").mkdir()

    assert _listdir_nohidden(tmp_path) == ["dir"]


def test_get_first_file_matching_pattern(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import (
        _get_first_file_matching_pattern,
    )

    (tmp_path / "foo_bar.nii.gz").touch()
    (tmp_path / "foo_bar_baz.nii").touch()

    assert (
        _get_first_file_matching_pattern(tmp_path, "foo*.nii*")
        == tmp_path / "foo_bar.nii.gz"
    )


@pytest.mark.parametrize(
    "pattern,msg",
    [("", "Pattern is not valid."), ("foo*.txt*", "No file matching pattern")],
)
def test_get_first_file_matching_pattern_error(tmp_path, pattern, msg):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import (
        _get_first_file_matching_pattern,
    )

    (tmp_path / "foo_bar.nii.gz").touch()
    (tmp_path / "foo_bar_baz.nii").touch()

    with pytest.raises(ValueError, match=msg):
        _get_first_file_matching_pattern(tmp_path, pattern)
