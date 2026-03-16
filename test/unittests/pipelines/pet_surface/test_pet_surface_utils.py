from pathlib import Path
from unittest import mock

import pytest


def _padding_to_3(number: int) -> str:
    number_string = str(number)
    return "0" * (3 - len(number_string)) + number_string


def _running_mris_expand_with_subprocess_mock(cmd: str):
    cmd_parts = cmd.split(" ")
    # cmd line contains : number of files to write -1 (-4), path for files (-3), base file name (-1)
    for n in range(0, int(cmd_parts[-4]) + 1):
        (Path(cmd_parts[-3]).parent / (cmd_parts[-1] + _padding_to_3(n))).touch()


def test_mris_expand(tmp_path):
    from clinica.pipelines.pet_surface.pet_surface_utils import mris_expand

    name = "lh.white"
    file = tmp_path / name
    file.touch()

    with mock.patch(
        "clinica.pipelines.pet_surface.pet_surface_utils._running_mris_expand_with_subprocess",
        wraps=_running_mris_expand_with_subprocess_mock,
    ) as mocked:
        mris_expand(str(file))
        mocked.assert_called_once()

    outputs = list(tmp_path.glob(rf"*{name}*"))
    assert len(outputs) == 7
    # todo : test which file was removed
