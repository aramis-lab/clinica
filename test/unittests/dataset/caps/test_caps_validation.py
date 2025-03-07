import json
import re

import pytest

from clinica.utils.exceptions import ClinicaCAPSError, ClinicaDatasetError


def test_check_caps_dataset_missing_description_error(tmp_path):
    from clinica.dataset import check_caps_dataset

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided is missing a dataset_description.json file."
        ),
    ):
        check_caps_dataset(tmp_path)


def test_check_caps_dataset_bad_description_error(tmp_path):
    from clinica.dataset import check_caps_dataset

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump(
            {"Name": "Example dataset", "BIDSVersion": "1.0.2", "CAPSVersion": "1.0.0"},
            fp,
        )
    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided has a badly formatted dataset_description.json file."
        ),
    ):
        check_caps_dataset(tmp_path)


def test_check_caps_dataset_bad_description_error_2(tmp_path):
    from clinica.dataset import DatasetType, check_caps_dataset

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump(
            {
                "Name": "Example dataset",
                "DatasetType": DatasetType.DERIVATIVE.value,
                "BIDSVersion": "1.0.2",
                "CAPSVersion": "1.0.0",
            },
            fp,
        )
    assert check_caps_dataset(tmp_path) == tmp_path

    (tmp_path / "sub-01").mkdir()
    with pytest.raises(
        ClinicaCAPSError,
        match="Your CAPS directory contains at least one folder whose name starts with 'sub-'.",
    ):
        check_caps_dataset(tmp_path)
