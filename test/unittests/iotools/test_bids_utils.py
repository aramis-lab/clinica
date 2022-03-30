import pytest
from typing import Union
from pathlib import Path
from string import Template

from clinica.iotools.bids_utils import (
    _write_readme, _write_bids_validator_config, _write_bidsignore
)


MODALITY_AGNOSTIC_FILE_WRITERS = {
    "readme": _write_readme,
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
    ("This BIDS directory was generated with Clinica v$version.\n"
     "More information on $website\n")
)


def _validate_file_and_content(
        file: Path,
        expected_content: str
) -> None:
    import os
    assert os.path.exists(file)
    assert file.read_text() == expected_content


@pytest.fixture
def expected_description_content(
        study_name: str,
        bids_version: Union[None, str],
) -> str:
    import json
    from clinica.iotools.bids_dataset_description import BIDS_VERSION
    expected_version = BIDS_VERSION if bids_version is None else bids_version
    desc_dict = {
        "Name": study_name,
        "BIDSVersion": expected_version,
        "DatasetType": "raw",
    }
    return json.dumps(desc_dict, indent=4)


@pytest.mark.parametrize("study_name", ["ADNI", "foo"])
@pytest.mark.parametrize("bids_version", [None, "1.6.0", "1.7.0"])
def test_write_bids_dataset_description(
        tmp_path,
        study_name,
        bids_version,
        expected_description_content,
):
    """Test function `_write_bids_dataset_description`.

    .. note::
        Tested independantly for convenience since it takes
        a different set of input parameters.

    """
    from clinica.iotools.bids_utils import _write_bids_dataset_description
    _write_bids_dataset_description(
        study_name, tmp_path, bids_version=bids_version
    )
    _validate_file_and_content(
        tmp_path / EXPECTED_MODALITY_AGNOSTIC_FILES["description"],
        expected_description_content
    )


def expected_readme_content() -> str:
    import clinica
    return EXPECTED_README_CONTENT.safe_substitute(
        version=clinica.__version__,
        website="http://www.clinica.run"
    )


def expected_validator_content() -> str:
    import json
    from clinica.iotools.bids_utils import BIDS_VALIDATOR_CONFIG
    return json.dumps(BIDS_VALIDATOR_CONFIG, indent=4)


@pytest.fixture
def expected_content(name: str) -> str:
    if name == "readme":
        return expected_readme_content()
    elif name == "bids-validator":
        return expected_validator_content()
    return "\n".join(["pet/", "conversion_info/"])


@pytest.mark.parametrize("name,writer", MODALITY_AGNOSTIC_FILE_WRITERS.items())
def test_modality_agnostic_file_writers(tmp_path, name, writer, expected_content):
    """Test helper functions of the function `write_modality_agnostic_files`."""
    writer(tmp_path)
    _validate_file_and_content(
        tmp_path / EXPECTED_MODALITY_AGNOSTIC_FILES[name],
        expected_content
    )


def test_write_modality_agnostic_files(tmp_path):
    """Test function `write_modality_agnostic_files`."""
    import os
    from clinica.iotools.bids_utils import write_modality_agnostic_files
    assert len(os.listdir(tmp_path)) == 0
    write_modality_agnostic_files("ADNI", tmp_path)
    files = os.listdir(tmp_path)
    assert len(files) == 4
    for _, v in EXPECTED_MODALITY_AGNOSTIC_FILES.items():
        assert v in files
