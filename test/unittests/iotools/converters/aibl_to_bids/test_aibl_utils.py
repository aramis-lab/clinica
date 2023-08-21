def test_listdir_nohidden(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import _listdir_nohidden

    (tmp_path / "file").touch()
    (tmp_path / ".hidden_file").touch()
    (tmp_path / "dir").mkdir()
    (tmp_path / ".hidden_dir").mkdir()

    assert _listdir_nohidden(tmp_path) == ["dir"]
