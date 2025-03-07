import json

import pytest

from clinica.dataset import DatasetType


@pytest.mark.parametrize("dataset_type", ["raw", "derivative"])
def test_get_dataset_type(tmp_path, dataset_type: str):
    from clinica.dataset import get_dataset_type

    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump(
            {"Name": "Test", "DatasetType": dataset_type, "BIDSVersion": "1.7.0"}, fp
        )

    assert get_dataset_type(tmp_path) == DatasetType(dataset_type)
