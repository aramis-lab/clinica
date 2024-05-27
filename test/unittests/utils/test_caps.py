import json

import pytest


def test_write_caps_dataset_description(tmp_path):
    from clinica.utils.caps import write_caps_dataset_description

    (tmp_path / "caps").mkdir()

    write_caps_dataset_description("foo", tmp_path / "caps")

    files = [f for f in (tmp_path / "caps").iterdir()]
    assert len(files) == 1
    assert json.loads(files[0].read_text()) == {
        "Name": "foo",
        "BidsVersion": "1.7.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
    }


def test_write_caps_dataset_description_specify_bids_and_caps_versions(tmp_path):
    from clinica.utils.caps import write_caps_dataset_description

    (tmp_path / "caps").mkdir()

    write_caps_dataset_description(
        "foo", tmp_path / "caps", bids_version="foobar", caps_version="2.0.0"
    )

    files = [f for f in (tmp_path / "caps").iterdir()]
    assert len(files) == 1
    assert json.loads(files[0].read_text()) == {
        "Name": "foo",
        "BidsVersion": "foobar",
        "CAPSVersion": "2.0.0",
        "DatasetType": "derivative",
    }


def test_read_caps_dataset_description(tmp_path):
    from clinica.utils.caps import (
        CAPSDatasetDescription,
        DatasetType,
        write_caps_dataset_description,
    )

    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()

    write_caps_dataset_description(
        "foo", caps_dir, bids_version="1.7.0", caps_version="1.0.0"
    )
    desc = CAPSDatasetDescription.from_file(caps_dir / "dataset_description.json")

    assert desc.name == "foo"
    assert desc.bids_version == "1.7.0"
    assert desc.caps_version == "1.0.0"
    assert desc.dataset_type == DatasetType.DERIVATIVE


def test_write_caps_dataset_description_error(tmp_path):
    from clinica.utils.caps import (
        CAPSDatasetDescription,
        DatasetType,
        write_caps_dataset_description,
    )
    from clinica.utils.exceptions import ClinicaCAPSError

    caps_dir = tmp_path / "caps"
    caps_dir.mkdir()

    write_caps_dataset_description(
        "foo", caps_dir, bids_version="1.7.0", caps_version="1.0.0"
    )
    # Re-writing the same description works
    write_caps_dataset_description(
        "foo", caps_dir, bids_version="1.7.0", caps_version="1.0.0"
    )
    # Re-writing the same description with a different name works
    write_caps_dataset_description(
        "bar", caps_dir, bids_version="1.7.0", caps_version="1.0.0"
    )

    desc = CAPSDatasetDescription.from_file(caps_dir / "dataset_description.json")

    assert desc.name == "foo + bar"
    assert desc.bids_version == "1.7.0"
    assert desc.caps_version == "1.0.0"
    assert desc.dataset_type == DatasetType.DERIVATIVE

    # But re-writing a different description raises an error
    with pytest.raises(
        ClinicaCAPSError,
        match=(
            f"Impossible to write the dataset_description.json file in {caps_dir} "
            "because it already exists and it contains incompatible metadata."
        ),
    ):
        write_caps_dataset_description(
            "bar", caps_dir, bids_version="1.7.1", caps_version="1.0.0"
        )
