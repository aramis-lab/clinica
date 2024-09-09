import json
from pathlib import Path


def test_build_dataset_table():
    from clinica.iotools.utils.describe import _build_dataset_table  # noqa
    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    table = _build_dataset_table(BIDSDatasetDescription("Test dataset BIDS"))

    assert table.title == "Dataset information"
    assert len(table.columns) == 3
    for i, value in enumerate(("Name", "Type", "BIDS Specifications Version")):
        assert table.columns[i].header == value
    for i, value in enumerate(("bright_cyan", "bright_yellow", "bright_magenta")):
        assert table.columns[i].style == value
    for i, value in enumerate(("Test dataset BIDS", "raw", "1.7.0")):
        assert table.columns[i]._cells[0] == value


def test_build_processing_table_bids():
    from clinica.iotools.utils.describe import _build_processing_table  # noqa
    from clinica.iotools.bids_dataset_description import BIDSDatasetDescription

    assert _build_processing_table(BIDSDatasetDescription("Test dataset BIDS")) is None


def build_test_dataset_caps(folder: Path):
    metadata = {
        "Name": "Test dataset CAPS",
        "BIDSVersion": "1.9.0",
        "CAPSVersion": "1.0.0",
        "DatasetType": "derivative",
        "Processing": [
            {
                "Name": "t1-linear",
                "Date": "2024-08-30T13:41:05.529799",
                "Author": "john.doe",
                "Machine": "machine-name",
                "InputPath": "/Users/john.doe/datasets/adni",
                "ClinicaVersion": "0.9.0",
                "Dependencies": [
                    {
                        "Name": "ants",
                        "VersionConstraint": ">=2.2.0",
                        "InstalledVersion": "2.3.5",
                    }
                ],
            },
            {
                "Name": "pet-linear",
                "Date": "2024-08-30T13:41:05.529799",
                "Author": "sheldon.cooper",
                "Machine": "sheldon-machine",
                "InputPath": "/Users/sheldon.cooper/work/datasets/oasis",
                "ClinicaVersion": "0.9.0",
                "Dependencies": [
                    {
                        "Name": "ants",
                        "VersionConstraint": ">=2.3.0",
                        "InstalledVersion": "2.5.1",
                    },
                    {
                        "Name": "spm",
                        "VersionConstraint": ">=12",
                        "InstalledVersion": "12.2.1",
                    },
                ],
            },
        ],
    }
    with open(folder / "dataset_description.json", "w") as fp:
        json.dump(metadata, fp)


def test_build_processing_table_caps(tmp_path):
    from clinica.iotools.utils.describe import _build_processing_table  # noqa
    from clinica.utils.caps import CAPSDatasetDescription

    (tmp_path / "bids").mkdir()
    build_test_dataset_caps(tmp_path / "bids")
    description = CAPSDatasetDescription.from_file(
        tmp_path / "bids" / "dataset_description.json"
    )
    table = _build_processing_table(description)

    assert table.title == "Processing information"
    assert len(table.columns) == 7
    assert len(table.rows) == 2
    for i, value in enumerate(
        (
            "Name",
            "Date",
            "Author",
            "Machine",
            "InputPath",
            "ClinicaVersion",
            "Dependencies",
        )
    ):
        assert table.columns[i].header == value
    for i, value in enumerate(
        (
            "bright_cyan",
            "bright_yellow",
            "bright_magenta",
            "bright_green",
            "bright_red",
            "bright_blue",
            "bright_green",
        )
    ):
        assert table.columns[i].style == value
    assert table.columns[0]._cells == ["t1-linear", "pet-linear"]
    assert table.columns[1]._cells == [
        "2024-08-30 13:41:05.529799",
        "2024-08-30 13:41:05.529799",
    ]
    assert table.columns[2]._cells == ["john.doe", "sheldon.cooper"]
    assert table.columns[3]._cells == ["machine-name", "sheldon-machine"]
    assert table.columns[4]._cells == [
        "/Users/john.doe/datasets/adni",
        "/Users/sheldon.cooper/work/datasets/oasis",
    ]
    assert table.columns[5]._cells == ["0.9.0", "0.9.0"]
    dependency_table_t1_linear, dependency_table_pet_linear = table.columns[6]._cells
    assert (
        len(dependency_table_t1_linear.columns)
        == len(dependency_table_pet_linear.columns)
        == 3
    )
    assert len(dependency_table_t1_linear.rows) == 1
    assert len(dependency_table_pet_linear.rows) == 2
    for dependency_table in table.columns[6]._cells:
        for i, value in enumerate(("Name", "VersionConstraint", "InstalledVersion")):
            assert dependency_table.columns[i].header == value
    for i, value in enumerate(("ants", ">=2.2.0", "2.3.5")):
        assert dependency_table_t1_linear.columns[i]._cells[0] == value
    for i, value in enumerate(
        (["ants", "spm"], [">=2.3.0", ">=12"], ["2.5.1", "12.2.1"])
    ):
        assert dependency_table_pet_linear.columns[i]._cells == value
