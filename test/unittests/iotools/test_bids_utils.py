from pathlib import Path
from string import Template
from typing import Union

import pytest

from clinica.iotools.bids_utils import (
    StudyName,
    _write_bids_validator_config,
    _write_bidsignore,
    _write_readme,
)

MODALITY_AGNOSTIC_FILE_WRITERS = {
    #    "readme": _write_readme,
    "bids-validator": _write_bids_validator_config,
    "bidsignore": _write_bidsignore,
}


EXPECTED_MODALITY_AGNOSTIC_FILES = {
    "description": "dataset_description.json",
    "readme": "README",
    "bids-validator": ".bids-validator-config.json",
    "bidsignore": ".bidsignore",
}


EXPECTED_README_CONTENT = Template(
    (
        "This BIDS directory was generated with Clinica v$version.\n"
        "More information on $website\n"
        "\n"
        "Study: $study\n"
        "\n"
        "description\n\n"
        "Find more about it and about the data user agreement: link"
    )
)


def test_get_bids_subjs_list(tmp_path):
    from clinica.iotools.bids_utils import get_bids_subjs_list

    (tmp_path / "file").touch()
    (tmp_path / "sub-03").touch()
    (tmp_path / "folder").mkdir()

    for sub in ("sub-01", "sub-02", "sub-16"):
        (tmp_path / sub).mkdir()

    assert set(get_bids_subjs_list(tmp_path)) == {"sub-01", "sub-02", "sub-16"}


@pytest.mark.parametrize(
    "input_string,expected",
    [
        ("foo", "foo"),
        ("foo_bar", "foobar"),
        ("foo bar", "foobar"),
        ("foo bar_baz", "foobarbaz"),
        ("foo-ba_r baz", "foobarbaz"),
    ],
)
def test_remove_space_and_symbols(input_string, expected):
    from clinica.iotools.bids_utils import remove_space_and_symbols

    assert remove_space_and_symbols(input_string) == expected


@pytest.mark.parametrize("compress", [True, False])
@pytest.mark.parametrize("sidecar", [True, False])
def test_build_dcm2niix_command(tmp_path, compress, sidecar):
    from clinica.iotools.bids_utils import _build_dcm2niix_command

    compress_flag = "y" if compress else "n"
    sidecar_flag = "y" if sidecar else "n"
    expected = ["dcm2niix", "-w", "0", "-f", "fmt", "-o", str(tmp_path / "out")]
    if compress:
        expected += ["-9"]
    expected += ["-z", compress_flag, "-b", sidecar_flag]
    if sidecar:
        expected += ["-ba", "y"]
    expected += [str(tmp_path / "in")]
    assert (
        _build_dcm2niix_command(
            tmp_path / "in",
            tmp_path / "out",
            "fmt",
            compress=compress,
            bids_sidecar=sidecar,
        )
        == expected
    )


def _validate_file_and_content(file: Path, expected_content: str) -> None:
    assert file.exists()
    assert file.read_text() == expected_content


@pytest.fixture
def expected_description_content(
    study_name: StudyName,
    bids_version: Union[None, str],
) -> str:
    import json

    from clinica.utils.bids import BIDS_VERSION

    expected_version = BIDS_VERSION if bids_version is None else bids_version
    desc_dict = {
        "Name": study_name.value,
        "BIDSVersion": expected_version,
        "DatasetType": "raw",
    }
    return json.dumps(desc_dict, indent=4)


@pytest.mark.parametrize("study_name", StudyName)
@pytest.mark.parametrize("bids_version", [None, "1.6.0", "1.7.0"])
def test_write_bids_dataset_description(
    tmp_path,
    study_name,
    bids_version,
    expected_description_content,
):
    """Test function `_write_bids_dataset_description`.

    .. note::
        Tested independently for convenience since it takes
        a different set of input parameters.

    """
    from clinica.iotools.bids_utils import _write_bids_dataset_description

    _write_bids_dataset_description(study_name, tmp_path, bids_version=bids_version)
    _validate_file_and_content(
        tmp_path / EXPECTED_MODALITY_AGNOSTIC_FILES["description"],
        expected_description_content,
    )


def expected_validator_content() -> str:
    import json

    from clinica.iotools.bids_utils import BIDS_VALIDATOR_CONFIG

    return json.dumps(BIDS_VALIDATOR_CONFIG, indent=4)


@pytest.fixture
def expected_content(name: str, study_name: StudyName) -> str:
    if name == "readme":
        return get_expected_readme_content(study_name)
    elif name == "bids-validator":
        return expected_validator_content()
    return "\n".join(["swi/", "conversion_info/"])


@pytest.mark.parametrize("study_name", StudyName)
@pytest.mark.parametrize("name,writer", MODALITY_AGNOSTIC_FILE_WRITERS.items())
def test_modality_agnostic_file_writers(
    tmp_path, study_name, name, writer, expected_content
):
    """Test helper functions of the function `write_modality_agnostic_files`."""
    writer(tmp_path)
    _validate_file_and_content(
        tmp_path / EXPECTED_MODALITY_AGNOSTIC_FILES[name], expected_content
    )


def test_write_modality_agnostic_files(tmp_path):
    """Test function `write_modality_agnostic_files`."""
    import os

    from clinica.iotools.bids_utils import write_modality_agnostic_files

    data_dict = {"link": "", "desc": ""}
    assert len(os.listdir(tmp_path)) == 0
    write_modality_agnostic_files(StudyName.ADNI, data_dict, tmp_path)
    files = os.listdir(tmp_path)
    assert len(files) == 4
    for _, v in EXPECTED_MODALITY_AGNOSTIC_FILES.items():
        assert v in files


@pytest.mark.parametrize("study_name", StudyName)
@pytest.mark.parametrize("bids_version", ["1.7.0"])
def test_write_bids_readme(
    tmp_path,
    study_name: StudyName,
    bids_version,
):
    """Test function `_write_bids_readme`.

    .. note::
        Tested independently for convenience since it takes
        a different set of input parameters.

    """

    data_dict = {"link": "link", "desc": "description"}
    _write_readme(study_name=study_name, data_dict=data_dict, bids_dir=tmp_path)
    _validate_file_and_content(
        file=tmp_path / EXPECTED_MODALITY_AGNOSTIC_FILES["readme"],
        expected_content=get_expected_readme_content(study_name),
    )


def get_expected_readme_content(study_name: StudyName) -> str:
    import clinica

    return EXPECTED_README_CONTENT.safe_substitute(
        version=clinica.__version__,
        website="https://www.clinica.run",
        study=study_name.value,
    )
