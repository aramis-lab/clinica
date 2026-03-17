import os
from pathlib import Path
from unittest import mock

import pytest


def _padding_to_3(number: int) -> str:
    number_string = str(number)
    return "0" * (3 - len(number_string)) + number_string


def _running_mris_expand_with_subprocess_mock(cmd: str):
    cmd_parts = cmd.split(" ")
    # cmd line contains : number of files to write -1 (-4), base file name (-1)
    for n in range(0, int(cmd_parts[-4]) + 1):
        (Path(os.getcwd()) / (cmd_parts[-1] + _padding_to_3(n))).touch()


@pytest.mark.parametrize("platform", ["darwin", "linux"])
def test_setting_mris_expand_cmd(tmp_path, platform):
    import sys

    from clinica.pipelines.pet_surface.pet_surface_utils import _setting_mris_expand_cmd

    in_surface = str(tmp_path / "lh.white")
    with mock.patch.object(sys, "platform", platform):
        result = _setting_mris_expand_cmd(in_surface)
        assert f"mris_expand -thickness -N 13 {in_surface} 0.65 lh.white_exp-" in result
        if platform == "darwin":
            assert "export" in result


def test_mris_expand(tmp_path):
    from clinica.pipelines.pet_surface.pet_surface_utils import mris_expand

    name = "lh.white"
    in_surface = tmp_path / name
    in_surface.touch()

    with mock.patch(
        "clinica.pipelines.pet_surface.pet_surface_utils._running_mris_expand_with_subprocess",
        wraps=_running_mris_expand_with_subprocess_mock,
    ) as mocked:
        mris_expand(str(in_surface))
        mocked.assert_called_once()

    outputs = list(Path(os.getcwd()).glob(rf"*{name}*"))
    assert len(outputs) == 7
    assert {_padding_to_3(x) for x in range(7, 14)} == set(
        map(lambda x: x.name.split("-")[-1], outputs)
    )
