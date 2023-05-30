from clinica.iotools.converters.aibl_to_bids import aibl_utils


def test_listdir_nohidden():
    from pathlib import Path
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as tmpdir:
        path = Path(tmpdir)
        (path / "file").touch()
        (path / ".hidden_file").touch()
        (path / "dir").mkdir()
        (path / ".hidden_dir").mkdir()
        result = aibl_utils.listdir_nohidden(str(path))
        assert result == ["dir"]
