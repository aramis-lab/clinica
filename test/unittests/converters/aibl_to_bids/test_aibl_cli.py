import pytest
from click.testing import CliRunner

from clinica.converters.aibl_to_bids.cli import cli


@pytest.fixture
def cli_paths(tmp_path):
    dataset = tmp_path / "imaging data"
    clinical = tmp_path / "clinical data"
    bids = tmp_path / "bids output"
    dataset.mkdir()
    clinical.mkdir()
    return [str(dataset), str(clinical), str(bids)]


@pytest.fixture
def mock_convert(mocker):
    mocker.patch("clinica.option.cpu_count", return_value=4)
    return mocker.patch("clinica.converters.aibl_to_bids._converter.convert")


@pytest.mark.parametrize(
    "options,clinical_data_only,overwrite",
    [
        ([], False, False),
        (["-c"], True, False),
        (["--clinical-data-only"], True, False),
        (["--overwrite"], False, True),
        (["--clinical-data-only", "--overwrite"], True, True),
    ],
)
def test_aibl_cli_forwards_options(
    cli_paths, mock_convert, options, clinical_data_only, overwrite
):
    result = CliRunner().invoke(cli, [*cli_paths, *options])

    assert result.exit_code == 0, result.output
    dataset, clinical, bids = cli_paths
    mock_convert.assert_called_once_with(
        dataset,
        bids,
        clinical,
        overwrite=overwrite,
        clinical_data_only=clinical_data_only,
        n_procs=3,
    )


@pytest.mark.parametrize("option", ["-np", "--n_procs"])
@pytest.mark.parametrize("n_procs", [1, 4])
def test_aibl_cli_forwards_n_procs(cli_paths, mock_convert, option, n_procs):
    result = CliRunner().invoke(cli, [*cli_paths, option, str(n_procs)])

    assert result.exit_code == 0, result.output
    dataset, clinical, bids = cli_paths
    mock_convert.assert_called_once_with(
        dataset,
        bids,
        clinical,
        overwrite=False,
        clinical_data_only=False,
        n_procs=n_procs,
    )


@pytest.mark.parametrize("n_procs", ["0", "-1", "5", "not-an-integer"])
def test_aibl_cli_rejects_invalid_n_procs(cli_paths, mock_convert, n_procs):
    result = CliRunner().invoke(cli, [*cli_paths, "--n_procs", n_procs])

    assert result.exit_code == 2
    assert "Invalid value" in result.output
    mock_convert.assert_not_called()


@pytest.mark.parametrize("input_index", [0, 1])
@pytest.mark.parametrize("existing_file", [False, True])
def test_aibl_cli_requires_input_directories(
    cli_paths, mock_convert, tmp_path, input_index, existing_file
):
    invalid_path = tmp_path / "invalid input"
    if existing_file:
        invalid_path.write_text("not a directory")
    cli_paths[input_index] = str(invalid_path)

    result = CliRunner().invoke(cli, cli_paths)

    assert result.exit_code == 2
    assert "Invalid value" in result.output
    mock_convert.assert_not_called()
