import pytest
from click.testing import CliRunner

from clinica.cmdline import cli

# Test to ensure that the help string at the command line is invoked without errors

# Test for the first level at the command line
@pytest.fixture(
    params=["convert", "generate", "iotools", "run",]
)
def generate_cli_first_lv(request):
    task = request.param
    return task


def test_first_lv(generate_cli_first_lv):
    runner = CliRunner()
    cli_input = generate_cli_first_lv
    print(f"Testing input cli {cli_input}")
    result = runner.invoke(cli, cli_input)
    assert result.exit_code == 0


# Test for the converters cli (second level)
@pytest.fixture(
    params=["adni-to-bids", "aibl-to-bids", "nifd-to-bids", "oasis-to-bids",]
)
def generate_cli_second_lv_convert(request):
    task = request.param
    return task


def test_second_lv_convert(generate_cli_second_lv_convert):
    runner = CliRunner()
    cli_input = generate_cli_second_lv_convert
    print(f"Testing input cli convert {cli_input}")
    result = runner.invoke(cli, f"convert {cli_input} -h")
    assert result.exit_code == 0


# Test for the iotools cli (second level)
@pytest.fixture(
    params=[
        "center-nifti",
        "check-missing-modalities",
        "check-missing-processing",
        "create-subjects-visits",
        "merge-tsv",
    ]
)
def generate_cli_second_lv_iotools(request):
    task = request.param
    return task


def test_second_lv_iotools(generate_cli_second_lv_iotools):
    runner = CliRunner()
    cli_input = generate_cli_second_lv_iotools
    print(f"Testing input cli iotools {cli_input}")
    result = runner.invoke(cli, f"iotools {cli_input} -h")
    assert result.exit_code == 0


# Test for the pipelines cli (second level)
@pytest.fixture(
    params=[
        "t1-freesurfer",
        "t1-volume",
        "t1-freesurfer-longitudinal",
        "t1-linear",
        "dwi-preprocessing-using-phasediff-fmap",
        "dwi-preprocessing-using-t1",
        "dwi-dti",
        "dwi-connectome",
        "pet-linear",
        "pet-volume",
        "pet-surface",
        "pet-surface-longitudinal",
        "deeplearning-prepare-data",
        "machinelearning-prepare-spatial-svm",
        "machinelearning-classification",
        "statistics-surface",
        "statistics-volume",
        "statistics-volume-correction",
        "t1-volume-existing-template",
        "t1-volume-tissue-segmentation",
        "t1-volume-create-dartel",
        "t1-volume-register-dartel",
        "t1-volume-dartel2mni",
        "t1-volume-parcellation",
        "t1-freesurfer-template",
        "t1-freesurfer-longitudinal-correction",
    ]
)
def generate_cli_second_lv_run(request):
    task = request.param
    return task


def test_second_lv_run(generate_cli_second_lv_run):
    runner = CliRunner()
    cli_input = generate_cli_second_lv_run
    print(f"Testing input cli run {cli_input}")
    result = runner.invoke(cli, f"run {cli_input} -h")
    assert result.exit_code == 0
