from pathlib import Path
from typing import Optional
from unittest.mock import patch

import pytest
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


@pytest.fixture
def input_argument_list(
    tmp_path: Path, option_name: Optional[str], option_value: Optional[str]
) -> list:
    input_argument_list = [
        str(tmp_path / "bids"),
        str(tmp_path / "output"),
    ]
    if option_name:
        input_argument_list.extend([option_name, option_value])
    return input_argument_list


@pytest.fixture
def expected_arguments_passed_to_mock(
    tmp_path: Path, option_name: Optional[str], option_value: Optional[str]
) -> list:
    expected_arguments_passed_to_mock = [
        str(tmp_path / "bids"),
        str(tmp_path / "output"),
    ]
    if option_name:
        expected_arguments_passed_to_mock.append(option_value)
    else:
        expected_arguments_passed_to_mock.append("missing_mods")
    return expected_arguments_passed_to_mock


@pytest.mark.parametrize(
    "option_name,option_value",
    [
        (None, None),
        ("--output_prefix", "foobar"),
        ("--output_prefix", ""),
    ],
)
@patch("clinica.iotools.data_handling.compute_missing_mods")
def test_check_missing_modalities(
    mock_compute_missing_mods,
    tmp_path,
    option_name,
    option_value,
    input_argument_list,
    expected_arguments_passed_to_mock,
):
    (tmp_path / "bids").mkdir()
    runner = CliRunner()
    result = runner.invoke(check_missing_modalities, input_argument_list)
    assert result.exit_code == 0
    mock_compute_missing_mods.assert_called_once_with(
        *expected_arguments_passed_to_mock
    )
