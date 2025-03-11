import json


def test_get_subjects_from_bids_dataset(tmp_path):
    from clinica.dataset import get_subjects_from_bids_dataset

    (tmp_path / "file").touch()
    (tmp_path / "sub-03").touch()
    (tmp_path / "folder").mkdir()
    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump({"Name": "Test", "DatasetType": "raw", "BIDSVersion": "1.7.0"}, fp)

    for sub in ("sub-01", "sub-02", "sub-16"):
        (tmp_path / sub).mkdir()

    assert set(get_subjects_from_bids_dataset(tmp_path)) == {
        "sub-01",
        "sub-02",
        "sub-16",
    }
