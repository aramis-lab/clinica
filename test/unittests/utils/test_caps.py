import datetime
import json
from pathlib import Path

import pytest
from packaging.version import Version


def mock_processing_metadata(mocker):
    """Processing metadata is specific to the user and machine.

    This mock makes sure tests are reproducible by setting:
        - date to be 2024-08-06T16:30:00
        - user name to be 'John Doe'
        - machine name to be 'my machine'
    """
    mocker.patch(
        "clinica.utils.caps._get_current_timestamp",
        return_value=datetime.datetime(2024, 8, 6, 16, 30, 0),
    )
    mocker.patch("clinica.utils.caps._get_username", return_value="John Doe")
    mocker.patch("clinica.utils.caps._get_machine_name", return_value="my machine")
    return mocker


def test_caps_processing_description(tmp_path, mocker):
    from clinica.utils.caps import CAPSProcessingDescription

    mocker = mock_processing_metadata(mocker)
    desc = CAPSProcessingDescription.from_values("foo", str(tmp_path / "input"))

    assert desc.name == "foo"
    assert desc.date == datetime.datetime(2024, 8, 6, 16, 30)
    assert desc.author == "John Doe"
    assert desc.machine == "my machine"
    assert desc.input_path == str(tmp_path / "input")
    assert json.loads(str(desc)) == {
        "Name": "foo",
        "Date": "2024-08-06T16:30:00",
        "Author": "John Doe",
        "Machine": "my machine",
        "InputPath": f"{tmp_path}/input",
    }
    with open(tmp_path / "dataset_description.json", "w") as fp:
        desc.write(fp)
    desc2 = CAPSProcessingDescription.from_file(tmp_path / "dataset_description.json")
    assert desc == desc2


def test_caps_dataset_description(tmp_path, mocker):
    from clinica.utils.caps import CAPSDatasetDescription

    mocker = mock_processing_metadata(mocker)
    mocker.patch(
        "clinica.utils.caps._generate_random_name", return_value="my caps dataset"
    )

    desc = CAPSDatasetDescription.from_values()

    assert desc.name == "my caps dataset"
    assert desc.bids_version == Version("1.7.0")
    assert desc.caps_version == Version("1.0.0")
    assert desc.processing == []

    desc.add_processing("processing-1", str(tmp_path / "bids"))

    assert len(desc.processing) == 1
    assert desc.has_processing("processing-1")
    assert not desc.has_processing("processing-2")
    proc = desc.get_processing("processing-1")
    assert proc.name == "processing-1"
    assert proc.date == datetime.datetime(2024, 8, 6, 16, 30, 0)
    assert proc.author == "John Doe"
    assert proc.machine == "my machine"
    assert proc.input_path == str(tmp_path / "bids")
    assert json.loads(str(desc)) == {
        "Name": "my caps dataset",
        "BIDSVersion": "1.7.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
        "Processing": [
            {
                "Name": "processing-1",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            }
        ],
    }
    desc.add_processing("processing-2", str(tmp_path / "bids"))
    assert len(desc.processing) == 2
    assert json.loads(str(desc)) == {
        "Name": "my caps dataset",
        "BIDSVersion": "1.7.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
        "Processing": [
            {
                "Name": "processing-1",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            },
            {
                "Name": "processing-2",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            },
        ],
    }
    desc.delete_processing("processing-1")
    assert len(desc.processing) == 1
    desc.delete_processing("processing-2")
    assert len(desc.processing) == 0


def initialize_input_dir(folder: Path):
    desc = {"Name": "Input dataset", "BIDSVersion": "1.7.0"}
    folder.mkdir(exist_ok=True, parents=True)
    with open(folder / "dataset_description.json", "w") as fp:
        json.dump(desc, fp)


def test_write_caps_dataset_description(tmp_path, mocker):
    from clinica.utils.caps import write_caps_dataset_description

    mocker = mock_processing_metadata(mocker)
    (tmp_path / "caps").mkdir()
    initialize_input_dir(tmp_path / "bids")

    write_caps_dataset_description(
        tmp_path / "bids",
        tmp_path / "caps",
        processing_name="foo",
        dataset_name="my CAPS dataset",
    )

    files = [f for f in (tmp_path / "caps").iterdir()]
    assert len(files) == 1
    assert json.loads(files[0].read_text()) == {
        "Name": "my CAPS dataset",
        "BIDSVersion": "1.7.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
        "Processing": [
            {
                "Name": "foo",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            }
        ],
    }


def test_write_caps_dataset_description_specify_bids_and_caps_versions(tmp_path):
    from clinica.utils.caps import write_caps_dataset_description
    from clinica.utils.exceptions import ClinicaBIDSError

    (tmp_path / "caps").mkdir()
    initialize_input_dir(tmp_path / "bids")

    with pytest.raises(
        ClinicaBIDSError,
        match=(
            f"The input dataset {tmp_path}/bids has BIDS specifications following version 1.7.0, "
            "while the BIDS specifications version asked for the CAPS creation is 1.18.23. "
            "Please make sure the versions are the same before processing."
        ),
    ):
        write_caps_dataset_description(
            tmp_path / "bids",
            tmp_path / "caps",
            processing_name="foo",
            bids_version="1.18.23",
            caps_version="2.0.0",
        )


def test_read_caps_dataset_description(tmp_path, mocker):
    from clinica.utils.caps import (
        CAPSDatasetDescription,
        write_caps_dataset_description,
    )
    from clinica.utils.inputs import DatasetType

    mocker = mock_processing_metadata(mocker)
    initialize_input_dir(tmp_path / "bids")
    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()

    write_caps_dataset_description(
        tmp_path / "bids",
        caps_dir,
        processing_name="foo",
        dataset_name="my CAPS dataset",
    )

    desc = CAPSDatasetDescription.from_file(caps_dir / "dataset_description.json")

    assert desc.name == "my CAPS dataset"
    assert desc.bids_version == Version("1.7.0")
    assert desc.caps_version == Version("1.0.0")
    assert desc.dataset_type == DatasetType.DERIVATIVE
    assert len(desc.processing) == 1
    proc = desc.get_processing("foo")
    assert proc.name == "foo"
    assert proc.author == "John Doe"
    assert proc.date == datetime.datetime(2024, 8, 6, 16, 30)
    assert proc.machine == "my machine"


def test_write_caps_dataset_description_renaming_gives_warning(tmp_path):
    from clinica.utils.caps import write_caps_dataset_description

    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()
    initialize_input_dir(tmp_path / "bids")

    write_caps_dataset_description(
        tmp_path / "bids",
        tmp_path / "caps",
        processing_name="foo",
        dataset_name="my CAPS dataset",
    )
    with pytest.warns(
        UserWarning,
        match=(
            f"The existing CAPS dataset, located at {tmp_path}/caps has a name 'my CAPS dataset' "
            "different from the new name 'my CAPS dataset 2'. The old name will be kept."
        ),
    ):
        write_caps_dataset_description(
            tmp_path / "bids",
            tmp_path / "caps",
            processing_name="foo",
            dataset_name="my CAPS dataset 2",
        )


def test_write_caps_dataset_description_version_mismatch_error(tmp_path):
    from clinica.utils.caps import write_caps_dataset_description
    from clinica.utils.exceptions import ClinicaCAPSError

    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()
    initialize_input_dir(tmp_path / "bids")

    # Write a first processing named 'foo', with a CAPS version of 1.0.1
    write_caps_dataset_description(
        tmp_path / "bids",
        tmp_path / "caps",
        processing_name="foo",
        dataset_name="my CAPS dataset",
        caps_version="1.0.1",
    )
    # Now, write a second processing, named 'bar', but with a CAPS version of 1.0.2
    with pytest.raises(
        ClinicaCAPSError,
        match=(
            f"Impossible to write the 'dataset_description.json' file in {tmp_path}/caps "
            "because it already exists and it contains incompatible metadata."
        ),
    ):
        write_caps_dataset_description(
            tmp_path / "bids",
            tmp_path / "caps",
            processing_name="bar",
            dataset_name="my CAPS dataset",
            caps_version="1.0.2",
        )


def test_write_caps_dataset_description_multiple_processing(tmp_path, mocker):
    from clinica.utils.caps import write_caps_dataset_description

    mocker = mock_processing_metadata(mocker)
    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()
    initialize_input_dir(tmp_path / "bids")

    # Write a first processing named 'foo'
    write_caps_dataset_description(
        tmp_path / "bids",
        caps_dir,
        processing_name="foo",
        dataset_name="my CAPS dataset",
    )
    # Write a second processing, named 'bar'
    write_caps_dataset_description(
        tmp_path / "bids",
        caps_dir,
        processing_name="bar",
        dataset_name="my CAPS dataset",
    )
    files = [f for f in (tmp_path / "caps").iterdir()]
    assert len(files) == 1
    assert json.loads(files[0].read_text()) == {
        "Name": "my CAPS dataset",
        "BIDSVersion": "1.7.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
        "Processing": [
            {
                "Name": "foo",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            },
            {
                "Name": "bar",
                "Date": "2024-08-06T16:30:00",
                "Author": "John Doe",
                "Machine": "my machine",
                "InputPath": f"{tmp_path}/bids",
            },
        ],
    }


@pytest.mark.parametrize(
    "version1,version2,policy,expected",
    [
        ("1.2.3", "1.2.3", "strict", True),
        ("1.2.3", "1.2.3", "minor", True),
        ("1.2.3", "1.2.3", "major", True),
        ("1.2.3", "1.2.4", "strict", False),
        ("1.2.3", "1.2.4", "minor", True),
        ("1.2.3", "1.2.4", "major", True),
        ("1.2.3", "1.3.3", "strict", False),
        ("1.2.3", "1.3.3", "minor", False),
        ("1.2.3", "1.3.3", "major", True),
        ("1.2.3", "2.2.3", "strict", False),
        ("1.2.3", "2.2.3", "minor", False),
        ("1.2.3", "2.2.3", "major", False),
    ],
)
def test_are_versions_compatible(version1, version2, policy, expected):
    from clinica.utils.caps import are_versions_compatible

    assert are_versions_compatible(version1, version2, policy=policy) is expected
