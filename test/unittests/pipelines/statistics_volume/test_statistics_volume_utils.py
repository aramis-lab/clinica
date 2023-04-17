import math
from functools import partial
from pathlib import Path
from typing import List

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
        _create_spm_output_folder,
    )

    (tmp_path / "spm" / "models").mkdir(parents=True)
    output_folder = Path(
        _create_spm_output_folder(tmp_path / "spm" / "models" / "my_fake_model.m")
    )

    assert output_folder == tmp_path / "spm" / "2_sample_t_test"
    assert output_folder.exists()

    (output_folder / "foo.m").touch()

    output_folder = Path(
        _create_spm_output_folder(tmp_path / "spm" / "models" / "my_fake_model.m")
    )

    assert output_folder == tmp_path / "spm" / "2_sample_t_test"
    assert output_folder.exists()
    assert not (output_folder / "foo.m").exists()


def test_set_output_and_groups(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _set_output_and_groups,
    )

    current_model_filename = tmp_path / "current_model.m"
    current_model_filename.write_text("foo bar baz\nfoo\nbar baz")
    template = "foo @OUTPUTDIR bar\nfoo bar @SCANS2\nfoo baz @SCANS1\nfoo bar baz"
    _set_output_and_groups(
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
        _convert_to_numeric,
    )

    assert _convert_to_numeric(input_data) == expected


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
        _unravel_list_for_matlab,
    )

    assert _unravel_list_for_matlab(input_list) == expected


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
    from clinica.pipelines.statistics_volume.statistics_volume_utils import _is_number

    assert _is_number(input_number) == expected


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
        _delete_last_line,
    )

    (tmp_path / "foo.tsv").write_text(input_data)

    _delete_last_line(tmp_path / "foo.tsv")

    assert (tmp_path / "foo.tsv").read_text() == expected


def test_delete_last_line_empty(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _delete_last_line,
    )

    (tmp_path / "foo.tsv").touch()

    _delete_last_line(tmp_path / "foo.tsv")

    assert (tmp_path / "foo.tsv").read_text() == ""


def test_write_covariate_lines_error(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _write_covariate_lines,
    )

    with pytest.raises(
        FileNotFoundError,
        match="Could not find file",
    ):
        _write_covariate_lines(tmp_path / "foo.txt", 1, "age", [1.6, 2.3])


def test_write_covariate_lines(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _write_covariate_lines,
    )

    (tmp_path / "foo.m").write_text("Initial content\nfoo bar")

    _write_covariate_lines(tmp_path / "foo.m", 1, "age", [1.6, 2.3])

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
        "clinica.pipelines.statistics_volume.statistics_volume_utils._run_matlab_script_with_matlab",
        return_value=None,
    )

    with pytest.raises(
        RuntimeError,
        match="Output matrix",
    ):
        run_m_script(m_file)


def test_run_m_script(tmp_path, mocker):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import run_m_script

    model_folder = tmp_path / "matlab" / "models"
    model_folder.mkdir(parents=True)
    m_file = model_folder / "script.m"
    m_file.touch()
    expected_output_folder = tmp_path / "matlab" / "2_sample_t_test"
    expected_output_folder.mkdir()
    expected_output_file = expected_output_folder / "SPM.mat"
    expected_output_file.touch()
    mocker.patch(
        "clinica.pipelines.statistics_volume.statistics_volume_utils._run_matlab_script_with_matlab",
        return_value=None,
    )

    assert run_m_script(m_file) == str(expected_output_file)


@pytest.mark.parametrize(
    "mocked_values,expected",
    [
        (
            [True, False],
            "cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch script.m",
        ),
        ([False, True], "$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch script.m"),
    ],
    ids=("running on MAC OS", "running on Linux"),
)
def test_get_matlab_standalone_command(mocker, mocked_values, expected):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _get_matlab_standalone_command,
    )

    for func, return_value in zip(
        ["_is_running_on_mac", "_is_running_on_linux"], mocked_values
    ):
        mocker.patch(
            f"clinica.pipelines.statistics_volume.statistics_volume_utils.{func}",
            return_value=return_value,
        )

    assert _get_matlab_standalone_command(Path("script.m")) == expected


def test_get_matlab_standalone_command_system_error(mocker):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _get_matlab_standalone_command,
    )

    for func, return_value in zip(
        ["_is_running_on_mac", "_is_running_on_linux"], [False, False]
    ):
        mocker.patch(
            f"clinica.pipelines.statistics_volume.statistics_volume_utils.{func}",
            return_value=return_value,
        )

    with pytest.raises(
        SystemError,
        match="Clinica only support Mac OS and Linux",
    ):
        _get_matlab_standalone_command(Path("script.m"))


@pytest.mark.parametrize(
    "files",
    [
        ["figure_1.txt", "figure_2.pdf"],
        ["figure_1.txt", "figure_2.pdf", "figure_3.png"],
    ],
    ids=("two files but no png", "only one png"),
)
def test_rename_spm_figures_error(tmp_path, files):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_figures,
    )

    create_spm_files(tmp_path, files)

    with pytest.raises(
        RuntimeError,
        match="Wrong number of SPM figures files",
    ):
        _rename_spm_figures(
            spm_dir=tmp_path / "spm",
            output_dir=tmp_path,
            parameters_for_new_filename_construction={"group_label": "group"},
        )


def create_spm_files(folder: Path, files: List[str]) -> List[str]:
    spm_folder = create_spm_folder(folder)
    for filename in files:
        (spm_folder / filename).touch()

    return files


def create_spm_folder(folder: Path) -> Path:
    spm_folder = folder / "spm"
    if not spm_folder.exists():
        spm_folder.mkdir()

    return spm_folder


create_spm_figures = partial(
    create_spm_files, files=["figure_001.png", "figure_013.png", "figure_666.png"]
)
create_spm_t_maps = partial(create_spm_files, files=["spmT_0001.nii", "spmT_0002.nii"])
create_other_spm_files = partial(
    create_spm_files, files=["ResMS.nii", "RPV.nii", "mask.nii"]
)
create_spm_contrast_files = partial(create_spm_files, files=["con_1.txt", "con_2.txt"])
create_spm_beta_files = partial(
    create_spm_files,
    files=[
        "beta_1.nii.gz",
        "beta_2.nii.gz",
        "beta_cov1.nii.gz",
        "beta_cov2.nii.gz",
        "beta_cov3.nii.gz",
    ],
)


def create_all_spm_files(folder: Path) -> List[str]:
    spm_file_creators = (
        create_spm_figures,
        create_spm_t_maps,
        create_spm_beta_files,
        create_spm_contrast_files,
        create_other_spm_files,
    )

    return [func(folder) for func in spm_file_creators]


def test_rename_spm_figures(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_figures,
    )

    _ = create_spm_figures(tmp_path)

    output_files = _rename_spm_figures(
        spm_dir=tmp_path / "spm",
        output_dir=tmp_path,
        parameters_for_new_filename_construction={"group_label": "foo"},
    )

    assert output_files == [
        str(tmp_path / f"group-foo_report-{i}.png") for i in (1, 13, 666)
    ]
    assert all([Path(f).is_file() for f in output_files])


@pytest.mark.parametrize(
    "files",
    [
        ["file1.txt", "file2.pdf"],
        ["spmT_map.nii.gz", "figure_2.pdf", "figure_3.png"],
        ["spmT_map.nii.gz", "figure_2.pdf", "spmT_map2.nii.gz", "spmT_map3.nii.gz"],
    ],
    ids=("two files but no spmT", "only one spmT", "three spmT"),
)
def test_rename_spm_t_maps_error(tmp_path, files):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_t_maps,
    )

    create_spm_files(tmp_path, files)
    params = {
        "fwhm": 8,
        "group_label": "group",
        "class_names": ["A", "B"],
        "measure": "measure",
    }

    with pytest.raises(
        RuntimeError,
        match="SPM t-map",
    ):
        _rename_spm_t_maps(
            spm_dir=tmp_path / "spm",
            output_dir=tmp_path,
            parameters_for_new_filename_construction=params,
        )


@pytest.mark.parametrize(
    "fwhm,expected1,expected2",
    [
        (
            8,
            "group-foo_A-lt-B_measure-measure_fwhm-8_TStatistics.nii",
            "group-foo_B-lt-A_measure-measure_fwhm-8_TStatistics.nii",
        ),
        (
            None,
            "group-foo_A-lt-B_measure-measure_TStatistics.nii",
            "group-foo_B-lt-A_measure-measure_TStatistics.nii",
        ),
    ],
)
def test_rename_spm_t_maps(tmp_path, fwhm, expected1, expected2):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_t_maps,
    )

    _ = create_spm_t_maps(tmp_path)
    params = {
        "fwhm": fwhm,
        "group_label": "foo",
        "class_names": ["A", "B"],
        "measure": "measure",
    }

    spm_t_maps = _rename_spm_t_maps(
        spm_dir=tmp_path / "spm",
        output_dir=tmp_path,
        parameters_for_new_filename_construction=params,
    )
    assert spm_t_maps == [str(tmp_path / expected1), str(tmp_path / expected2)]


@pytest.mark.parametrize(
    "files",
    [
        ["file1.txt", "file2.pdf"],
        ["con_1.txt", "figure_2.pdf", "figure_3.png"],
        ["con_1.txt", "figure_2.pdf", "con_2.txt", "con_3.txt"],
    ],
    ids=("two files but no contrast", "only one contrast", "three contrasts"),
)
def test_rename_spm_contrast_files_error(tmp_path, files):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_contrast_files,
    )

    create_spm_files(tmp_path, files)
    params = {"group_label": "foo", "class_names": ["A", "B"], "measure": "measure"}

    with pytest.raises(
        RuntimeError,
        match="Wrong number of SPM contrast files",
    ):
        _rename_spm_contrast_files(
            spm_dir=tmp_path / "spm",
            output_dir=tmp_path,
            parameters_for_new_filename_construction=params,
        )


def test_rename_spm_contrast_files(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_spm_contrast_files,
    )

    _ = create_spm_contrast_files(tmp_path)
    params = {"group_label": "foo", "class_names": ["A", "B"], "measure": "measure"}

    contrast_files = _rename_spm_contrast_files(
        spm_dir=tmp_path / "spm",
        output_dir=tmp_path,
        parameters_for_new_filename_construction=params,
    )

    assert len(contrast_files) == 2
    assert contrast_files[0] == str(
        tmp_path / "group-foo_A-lt-B_measure-measure_contrast.nii"
    )
    assert contrast_files[1] == str(
        tmp_path / "group-foo_B-lt-A_measure-measure_contrast.nii"
    )
    assert all([Path(contrast_file).is_file() for contrast_file in contrast_files])


@pytest.mark.parametrize(
    "covariates",
    [["cov1", "cov2"], ["cov1", "cov2", "cov3"]],
    ids=("two covariates", "three covariates"),
)
@pytest.mark.parametrize(
    "files",
    [
        ["file1.txt", "file2.pdf"],
        ["beta_1.nii.gz", "file_1.pdf", "file_2.png"],
        ["beta_1.nii.gz", "file2.pdf", "beta_2.nii.gz", "beta_3.nii.gz"],
    ],
    ids=("two files but no beta file", "only one beta file", "three beta files"),
)
def test_rename_beta_files_error(tmp_path, files, covariates):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_beta_files,
    )

    create_spm_files(tmp_path, files)
    params = {"covariates": covariates, "class_names": ["A", "B"]}

    with pytest.raises(
        RuntimeError,
        match="Wrong number of SPM betas files",
    ):
        _rename_beta_files(
            spm_dir=tmp_path / "spm",
            output_dir=tmp_path,
            parameters_for_new_filename_construction=params,
        )


def test_rename_beta_files(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        _rename_beta_files,
    )

    _ = create_spm_beta_files(tmp_path)
    params = {"covariates": ["cov1", "cov2", "cov3"], "class_names": ["A", "B"]}

    new_files = _rename_beta_files(
        spm_dir=tmp_path / "spm",
        output_dir=tmp_path,
        parameters_for_new_filename_construction=params,
    )

    assert set(new_files) == {
        str(tmp_path / filename)
        for filename in ("A.nii", "B.nii", "cov1.nii", "cov2.nii", "cov3.nii")
    }
    assert all([Path(new_file).is_file() for new_file in new_files])


def test_copy_files(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import _copy_files

    for folder in ("sources", "destinations"):
        (tmp_path / folder).mkdir()
    source_destination_mapping = {}
    sources = [tmp_path / "sources" / f"file_{i}" for i in (1, 3, 6, 42)]
    for source in sources:
        source.touch()
        source_destination_mapping[str(source)] = str(
            tmp_path / "destinations" / source.name
        )

    _copy_files(source_destination_mapping)

    for source, destination in source_destination_mapping.items():
        assert Path(source).is_file()
        assert Path(destination).is_file()


def test_copy_and_rename_spm_output_files_error(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        copy_and_rename_spm_output_files,
    )

    with pytest.raises(
        RuntimeError,
        match=" SPM matrix",
    ):
        copy_and_rename_spm_output_files(
            str(tmp_path), ["A", "B"], ["cov1", "cov2"], "group", 8, "measure"
        )

    with pytest.raises(
        RuntimeError,
        match="output folder",
    ):
        copy_and_rename_spm_output_files(
            str(tmp_path / "spm" / "matrix.mat"),
            ["A", "B"],
            ["cov1", "cov2"],
            "group",
            8,
            "measure",
        )


def test_copy_and_rename_spm_output_files(tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        copy_and_rename_spm_output_files,
    )

    _ = create_all_spm_files(tmp_path)

    matrix_filename = tmp_path / "spm" / "matrix.mat"
    matrix_filename.touch()

    res = copy_and_rename_spm_output_files(
        str(matrix_filename),
        ["A", "B"],
        ["cov1", "cov2", "cov3"],
        "group",
        8,
        "measure",
        output_dir=str(tmp_path),
    )

    assert res == [
        [
            str(tmp_path / "group-group_A-lt-B_measure-measure_fwhm-8_TStatistics.nii"),
            str(tmp_path / "group-group_B-lt-A_measure-measure_fwhm-8_TStatistics.nii"),
        ],
        [
            str(tmp_path / "group-group_report-1.png"),
            str(tmp_path / "group-group_report-13.png"),
            str(tmp_path / "group-group_report-666.png"),
        ],
        [
            str(tmp_path / "group-group_VarianceError.nii"),
            str(tmp_path / "resels_per_voxel.nii"),
            str(tmp_path / "included_voxel_mask.nii"),
        ],
        [
            str(tmp_path / "A.nii"),
            str(tmp_path / "B.nii"),
            str(tmp_path / "cov1.nii"),
            str(tmp_path / "cov2.nii"),
            str(tmp_path / "cov3.nii"),
        ],
        [
            str(tmp_path / "group-group_A-lt-B_measure-measure_contrast.nii"),
            str(tmp_path / "group-group_B-lt-A_measure-measure_contrast.nii"),
        ],
    ]
