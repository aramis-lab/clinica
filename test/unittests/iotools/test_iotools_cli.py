from unittest.mock import patch

from click.testing import CliRunner

from clinica.iotools.cli import check_missing_modalities, check_missing_processing


@patch("clinica.iotools.data_handling.compute_missing_processing")
def test_check_missing_processing(mock_compute_missing_processing, tmp_path):
    (tmp_path / "bids").mkdir()
    (tmp_path / "caps").mkdir()
    runner = CliRunner()
    result = runner.invoke(
        check_missing_processing,
        [
            str(tmp_path / "bids"),
            str(tmp_path / "caps"),
            str(tmp_path / "outfile"),
        ],
    )
    assert result.exit_code == 0
    mock_compute_missing_processing.assert_called_once_with(
        str(tmp_path / "bids"),
        str(tmp_path / "caps"),
        str(tmp_path / "outfile"),
    )


@patch("clinica.iotools.data_handling.compute_missing_mods")
def test_check_missing_modalities(mock_compute_missing_mods, tmp_path):
    (tmp_path / "bids").mkdir()
    runner = CliRunner()
    result = runner.invoke(
        check_missing_modalities,
        [
            str(tmp_path / "bids"),
            str(tmp_path / "output"),
        ],
    )
    assert result.exit_code == 0
    mock_compute_missing_mods.assert_called_once_with(
        str(tmp_path / "bids"),
        str(tmp_path / "output"),
        "missing_mods",
    )
