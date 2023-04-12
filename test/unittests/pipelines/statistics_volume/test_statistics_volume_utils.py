import math
from pathlib import Path

import numpy as np
import pytest

from clinica.utils.exceptions import ClinicaException

TSV_COLUMN_NAMES = ("participant_id", "session_id", "group", "age", "sex")


def _get_tsv_header() -> str:
    return "\t".join(TSV_COLUMN_NAMES)


@pytest.fixture
def tsv_header():
    return _get_tsv_header()


@pytest.fixture
def tsv_data_two_classes():
    return _get_tsv_data_two_classes()


def _get_tsv_data_two_classes():
    data = (
        "sub-01\tses-M000\tA\t78,0\tFemale\n"
        "sub-02\tses-M000\tB\t82,0\tMale\n"
        "sub-03\tses-M000\tA\t79,0\tFemale\n"
        "sub-03\tses-M006\tB\t80,0\tFemale"
    )
    return f"{_get_tsv_header()}\n{data}"


def test_get_group_1_and_2_wrong_contrast(tmp_path, tsv_header):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        get_group_1_and_2,
    )

    (tmp_path / "groups.tsv").write_text(
        f"{tsv_header}\nsub-01\tses-M000\tCN\t78,0\tFemale"
    )

    with pytest.raises(
        ClinicaException,
        match="foo is not present",
    ):
        get_group_1_and_2(tmp_path / "groups.tsv", "foo")


@pytest.mark.parametrize(
    "data",
    [
        "sub-01\tses-M000\tA\t78,0\tFemale",
        (
            "sub-01\tses-M000\tA\t78,0\tFemale\n"
            "sub-02\tses-M000\tB\t82,0\tMale\n"
            "sub-03\tses-M000\tC\t79,0\tFemale"
        ),
    ],
    ids=("only one class", "three classes"),
)
def test_get_group_1_and_2_wrong_number_of_classes(tmp_path, tsv_header, data):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        get_group_1_and_2,
    )

    (tmp_path / "groups.tsv").write_text(f"{tsv_header}\n{data}")

    with pytest.raises(
        ClinicaException,
        match="It must exist only 2 classes in the column",
    ):
        get_group_1_and_2(tmp_path / "groups.tsv", "group")


def test_get_group_1_and_2(tmp_path, tsv_data_two_classes):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        get_group_1_and_2,
    )

    (tmp_path / "groups.tsv").write_text(tsv_data_two_classes)

    first_group_idx, second_group_idx, class_names = get_group_1_and_2(
        tmp_path / "groups.tsv", "group"
    )

    assert first_group_idx == [0, 2]
    assert second_group_idx == [1, 3]
    assert class_names == ["A", "B"]


def test_create_spm_output_folder(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        create_spm_output_folder,
    )

    (tmp_path / "spm" / "models").mkdir(parents=True)
    output_folder = Path(
        create_spm_output_folder(tmp_path / "spm" / "models" / "my_fake_model.m")
    )

    assert output_folder == tmp_path / "spm" / "2_sample_t_test"
    assert output_folder.exists()

    (output_folder / "foo.m").touch()

    output_folder = Path(
        create_spm_output_folder(tmp_path / "spm" / "models" / "my_fake_model.m")
    )

    assert output_folder == tmp_path / "spm" / "2_sample_t_test"
    assert output_folder.exists()
    assert not (output_folder / "foo.m").exists()


def test_set_output_and_groups(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        set_output_and_groups,
    )

    current_model_filename = tmp_path / "current_model.m"
    current_model_filename.write_text("foo bar baz\nfoo\nbar baz")
    template = "foo @OUTPUTDIR bar\nfoo bar @SCANS2\nfoo baz @SCANS1\nfoo bar baz"
    set_output_and_groups(
        str(tmp_path / "output_folder"),
        str(current_model_filename),
        [f"file_{i}" for i in range(1, 5)],
        [1, 3],
        [0, 2],
        template,
    )
    assert current_model_filename.exists()
    lines = current_model_filename.read_text().split("\n")
    assert len(lines) == 4
    for f in ("file_1", "file_3"):
        assert f in lines[1]
    for f in ("file_2", "file_4"):
        assert f in lines[2]
    for tag in ("@OUTPUTDIR", "@SCANS1", "@SCANS2"):
        assert tag not in current_model_filename.read_text()


@pytest.mark.parametrize(
    "input_data,expected",
    [
        (["1", "3", "6"], [1, 3, 6]),
        (["1.0", "2,6", "19"], [1.0, 2.6, 19]),
        (["foo", "bar", "foo", "baz", "baz"], [2, 0, 2, 1, 1]),
    ],
)
def test_convert_to_numeric(input_data, expected):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        convert_to_numeric,
    )

    assert convert_to_numeric(input_data) == expected


@pytest.mark.parametrize(
    "input_list,expected",
    [
        ([], "''"),
        (["a"], "'a'"),
        (["a", "b"], "'a','b'"),
        (["foo", "bar", "baz"], "'foo','bar','baz'"),
    ],
)
def test_unravel_list_for_matlab(input_list, expected):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        unravel_list_for_matlab,
    )

    assert unravel_list_for_matlab(input_list) == expected


@pytest.mark.parametrize(
    "input_number,expected",
    [
        ("", False),
        (0, True),
        (-1, True),
        (None, False),
        (True, True),
        (False, True),
        (3.14, True),
        (math.pi, True),
        (np.nan, True),  # float(nan) = nan and does not raise...
    ],
)
def test_is_number(input_number, expected):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import is_number

    assert is_number(input_number) == expected


@pytest.mark.parametrize(
    "input_data,expected",
    [
        ("", ""),
        ("foo bar baz", "foo bar baz"),
        ("foo bar baz\n", "foo bar baz\n"),
        ("foo bar\nbaz", "foo bar"),
        ("foo bar\nbaz\n", "foo bar"),
        (
            _get_tsv_data_two_classes(),
            "\n".join(_get_tsv_data_two_classes().split("\n")[:-1]),
        ),
    ],
    ids=(
        "empty content",
        "one line",
        "one line with new line",
        "two lines",
        "two lines with new line",
        "multiple lines",
    ),
)
def test_delete_last_line(tmp_path, input_data, expected):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        delete_last_line,
    )

    (tmp_path / "foo.tsv").write_text(input_data)

    delete_last_line(tmp_path / "foo.tsv")

    assert (tmp_path / "foo.tsv").read_text() == expected


def test_delete_last_line_empty(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        delete_last_line,
    )

    (tmp_path / "foo.tsv").touch()

    delete_last_line(tmp_path / "foo.tsv")

    assert (tmp_path / "foo.tsv").read_text() == ""


def test_write_covariate_lines_error(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        write_covariate_lines,
    )

    with pytest.raises(
        FileNotFoundError,
        match="Could not find file",
    ):
        write_covariate_lines(tmp_path / "foo.txt", 1, "age", [1.6, 2.3])


def test_write_covariate_lines(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        write_covariate_lines,
    )

    (tmp_path / "foo.m").write_text("Initial content\nfoo bar")

    write_covariate_lines(tmp_path / "foo.m", 1, "age", [1.6, 2.3])

    assert (tmp_path / "foo.m").read_text() == (
        "Initial content\n"
        "foo barmatlabbatch{1}.spm.stats.factorial_design.cov(1).c = [1.6 2.3];\n"
        "matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';\n"
        "matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;\n"
        "matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;\n"
    )


def test_run_m_script_file_not_found_error(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    with pytest.raises(
        FileNotFoundError,
        match="File",
    ):
        run_m_script(tmp_path / "script.m")


def test_run_m_script_file_extension_error(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    (tmp_path / "script.txt").touch()

    with pytest.raises(
        ValueError,
        match="is not a Matlab file",
    ):
        run_m_script(tmp_path / "script.txt")


def test_run_m_script_no_output_file_error(tmp_path, mocker):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    m_file = tmp_path / "script.m"
    m_file.touch()
    mocker.patch(
        "clinica.pipelines.statistics_volume.statistics_volume_utils.run_matlab_script_with_matlab",
        return_value=m_file,
    )

    with pytest.raises(
        RuntimeError,
        match="Output matrix",
    ):
        run_m_script(m_file)
