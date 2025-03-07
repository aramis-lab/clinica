import json
import re

import pytest

from clinica.utils.exceptions import ClinicaBIDSError, ClinicaDatasetError


def test_check_bids_dataset_empty_folder_error(tmp_path):
    from clinica.dataset import check_bids_dataset

    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided is missing a dataset_description.json file."
        ),
    ):
        check_bids_dataset(tmp_path)


def test_check_bids_dataset_empty_description_error(tmp_path):
    from clinica.dataset import check_bids_dataset

    (tmp_path / "dataset_description.json").touch()

    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided has a badly formatted dataset_description.json file."
        ),
    ):
        check_bids_dataset(tmp_path)


def test_check_bids_dataset_wrong_description_error(tmp_path):
    from clinica.dataset import check_bids_dataset

    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump({"foo": "bar"}, fp)

    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided has a badly formatted dataset_description.json file."
        ),
    ):
        check_bids_dataset(tmp_path)


def test_check_bids_dataset_no_subject_error(tmp_path):
    from clinica.dataset import check_bids_dataset

    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump({"Name": "Test", "DatasetType": "raw", "BIDSVersion": "1.7.0"}, fp)
    (tmp_path / "subjects").mkdir()

    with pytest.raises(
        ClinicaBIDSError,
        match=(
            "Your BIDS directory does not contains a single folder whose name starts with 'sub-'. "
            "Check that your folder follow BIDS standard."
        ),
    ):
        check_bids_dataset(tmp_path)
